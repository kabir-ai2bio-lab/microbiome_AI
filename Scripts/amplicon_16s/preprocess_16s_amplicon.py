#!/usr/bin/env python3
"""
16S rRNA amplicon preprocessing pipeline (local / WSL2 / Linux version).

Local port of the Colab notebook used for PRJNA496408. It turns public SRA
paired-end 16S runs into ML-ready genus abundance matrices:

    1. metadata   - NCBI RunInfo (+ BioSample attributes), include filter, provisional labels
    2. select     - choose runs for the mode (dry-run / validation / full)
    3. download   - prefetch + fasterq-dump + pigz, several runs in parallel
    4. import     - QIIME 2 manifest import + demux summary
    5. cutadapt   - primer removal (discard untrimmed pairs); skipped if primers: null
    6. figaro     - --figaro: FIGARO; with dada2.truncation: figaro it sets DADA2 trunc/maxEE
    7. qc         - mandatory per-sample primer retention / truncation gate
    8. dada2      - ONE joint DADA2 run across all selected samples
    9. dada2-stats- per-sample retention report + gate
   10. reference  - SILVA 138 NR99 via RESCRIPt (cached, shared across modes)
   11. taxonomy   - VSEARCH consensus classification
   12. genus      - collapse to genus, write count/relative matrices
   13. finalize   - reconcile counts, provenance, step timings, export

Every step is resumable: finished outputs are detected and skipped, and QIIME
artifacts are written to *.partial.qza first and renamed only on success.

All CPU-heavy tools (fasterq-dump, pigz, Cutadapt, DADA2, VSEARCH) use
--threads, which defaults to every CPU visible to the process.

Examples:
    # Check tools, CPUs, memory and disk; print the plan; run nothing
    python preprocess_16s_amplicon.py --config configs/PRJNA496408.yaml --mode dry-run --preflight

    # Small smoke test: 3 samples x 5000 read pairs, end to end
    python preprocess_16s_amplicon.py --config configs/PRJNA496408.yaml --mode dry-run

    # Full 93-sample joint run, also running FIGARO, copying final tables to Windows
    python preprocess_16s_amplicon.py --config configs/PRJNA496408.yaml --mode full \
        --figaro --export-dir /mnt/e/BiomeFM/outputs
"""

import argparse
import gzip
import hashlib
import io
import json
import logging
import os
import random
import re
import shutil
import subprocess
import sys
import time
import urllib.parse
import urllib.request
import zipfile
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from datetime import datetime
from itertools import zip_longest
from pathlib import Path
from typing import Callable, Dict, List, Optional

import numpy as np
import pandas as pd
import yaml

logger = logging.getLogger("amplicon_16s")

STEPS = [
    "metadata", "select", "download", "import", "cutadapt", "figaro", "qc",
    "dada2", "dada2-stats", "reference", "taxonomy", "genus", "finalize",
]
MODES = ("dry-run", "validation", "full")
GENUS_RE = r"(?:^|;\s*)g__([^;]+)"
# QIIME 2 stores paired reads as <sample>_<barcode>_L001_R<mate>_001.fastq.gz
QIIME_FASTQ_RE = re.compile(r"^(?P<sid>.+)_\d+_L\d{3}_R(?P<mate>[12])_001\.fastq\.gz$")


def available_cpus() -> int:
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count() or 1


def available_memory_gb() -> Optional[float]:
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemTotal:"):
                return int(line.split()[1]) / 1024 ** 2
    except OSError:
        pass
    return None


def partial_path(path: Path) -> Path:
    """foo.qza -> foo.partial.qza (QIIME 2 requires the real extension)."""
    return path.with_name(f"{path.stem}.partial{path.suffix}")


def atomic_copy(src: Path, dst: Path) -> None:
    src, dst = Path(src), Path(dst)
    if not src.exists() or src.stat().st_size == 0:
        raise FileNotFoundError(src)
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(str(dst) + ".partial")
    shutil.copy2(src, tmp)
    os.replace(tmp, dst)


def nonempty(path: Path, min_size: int = 1) -> bool:
    return path.exists() and path.stat().st_size >= min_size


def assign_labels(values: pd.Series, lab: dict) -> pd.Series:
    """Provisional labels from one metadata column.

    Either `regex` (first capture group) + optional `map` (e.g. CON -> NC), or `patterns`
    ({group: regex}; the first matching pattern wins). Unmatched values become NaN."""
    values = values.astype("string")
    if lab.get("patterns"):
        out = pd.Series(pd.NA, index=values.index, dtype="object")
        for group, rx in lab["patterns"].items():
            hit = out.isna() & values.str.match(rx, na=False)
            out[hit] = group
        return out
    return values.str.extract(lab["regex"], expand=False).replace(lab.get("map", {}))


def apply_filters(info: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Optional run filters: `include` {column: value or [values]} (exact match), e.g. keep stool and drop
    blood runs, and `include_regex` {column: regex}, e.g. one sequencing batch of a project."""
    for col, allowed in (cfg.get("include") or {}).items():
        allowed = [allowed] if isinstance(allowed, str) else list(allowed)
        before = len(info)
        info = info.loc[info[col].isin(allowed)].copy()
        logger.info("include %s in %s: kept %d of %d runs", col, allowed, len(info), before)
    for col, rx in (cfg.get("include_regex") or {}).items():
        before = len(info)
        info = info.loc[info[col].astype("string").str.match(rx, na=False)].copy()
        logger.info("include %s ~ %s: kept %d of %d runs", col, rx, len(info), before)
    return info


def dataset_name(cfg: dict) -> str:
    """Name for work/output folders: `dataset` (e.g. one sequencing batch of a project) or the BioProject."""
    return cfg.get("dataset", cfg["project"])


IUPAC_RE = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
            "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"}


def primer_regex(primer: str, max_offset: int) -> "re.Pattern":
    """Primer (IUPAC) starting within the first max_offset bases of a read."""
    return re.compile(r"^[ACGTN]{0,%d}" % max_offset + "".join(IUPAC_RE[c] for c in primer.upper()))


def _orient_pair(src_f: str, src_r: str, dst_f: str, dst_r: str, p1: str, p2: str, max_offset: int) -> dict:
    """Mixed-orientation libraries: swap mates when the R2 primer is on read 1 and the R1 primer on read 2,
    so every pair has the forward primer on R1. Pairs with neither pattern are kept for Cutadapt to judge."""
    fwd, rev = primer_regex(p1, max_offset), primer_regex(p2, max_offset)
    n = {"pairs": 0, "swapped": 0}
    tmp_f, tmp_r = dst_f + ".partial", dst_r + ".partial"
    with gzip.open(src_f, "rt") as inf, gzip.open(src_r, "rt") as inr, \
            gzip.open(tmp_f, "wt", compresslevel=1) as outf, gzip.open(tmp_r, "wt", compresslevel=1) as outr:
        for rec_f, rec_r in zip_longest(_iter_fastq(inf, src_f), _iter_fastq(inr, src_r)):
            if rec_f is None or rec_r is None:
                raise RuntimeError(f"Unpaired read counts in {src_f}")
            n["pairs"] += 1
            if not fwd.match(rec_f[1]) and rev.match(rec_f[1]) and fwd.match(rec_r[1]):
                rec_f, rec_r = rec_r, rec_f
                n["swapped"] += 1
            for out, (h, s, p, q) in ((outf, rec_f), (outr, rec_r)):
                out.write(f"{h}{s}\n{p}{q}\n")
    os.replace(tmp_f, dst_f)
    os.replace(tmp_r, dst_r)
    return n


# ---------------------------------------------------------------------------
# Worker functions (module level so ProcessPoolExecutor can pickle them)
# ---------------------------------------------------------------------------

def _iter_fastq(handle, label: str):
    """Yield (header, seq, plus, qual) from a text FASTQ handle, validating records."""
    while True:
        header = handle.readline()
        if not header:
            return
        seq, plus, qual = handle.readline(), handle.readline(), handle.readline()
        if not (seq and plus and qual):
            raise RuntimeError(f"Incomplete FASTQ record in {label}")
        seq, qual = seq.rstrip("\r\n"), qual.rstrip("\r\n")
        if len(seq) != len(qual):
            raise RuntimeError(f"Sequence/quality length mismatch in {label}")
        yield header, seq, plus, qual


def _qza_pair_stats(qza_path: str, member_f: str, member_r: Optional[str],
                    fmin: int, rmin: int) -> dict:
    """Count reads (pairs, or single reads if member_r is None) and truncation eligibility."""
    count = eligible = 0
    shortest_f = shortest_r = 10 ** 9
    with zipfile.ZipFile(qza_path) as zf:
        members = [m for m in (member_f, member_r) if m]
        handles = [gzip.open(zf.open(m), "rt") for m in members]
        try:
            for recs in zip_longest(*(_iter_fastq(h, m) for h, m in zip(handles, members))):
                if any(r is None for r in recs):
                    raise RuntimeError(f"Unpaired read counts in {member_f}")
                lf = len(recs[0][1])
                lr = len(recs[1][1]) if len(recs) > 1 else rmin
                count += 1
                shortest_f, shortest_r = min(shortest_f, lf), min(shortest_r, lr)
                eligible += int(lf >= fmin and lr >= rmin)
        finally:
            for h in handles:
                h.close()
    return {"pairs": count, "eligible": eligible,
            "eligible_fraction": eligible / count if count else 0.0,
            "min_f_len": shortest_f if count else 0,
            "min_r_len": (shortest_r if member_r else 0) if count else 0}


def _subsample(srcs: List[str], dsts: List[str], n_pairs: int) -> int:
    """Copy the first n_pairs reads / read pairs (dry-run only; a smoke test, not a random sample)."""
    written = 0
    tmps = [d + ".partial" for d in dsts]
    ins = [gzip.open(s, "rt") for s in srcs]
    outs = [gzip.open(t, "wt", compresslevel=1) for t in tmps]
    try:
        for recs in zip(*(_iter_fastq(h, s) for h, s in zip(ins, srcs))):
            for out, (h, s, p, q) in zip(outs, recs):
                out.write(f"{h}{s}\n{p}{q}\n")
            written += 1
            if written >= n_pairs:
                break
    finally:
        for f in ins + outs:
            f.close()
    for t, d in zip(tmps, dsts):
        os.replace(t, d)
    return written


def _length_profile(src_f: str, src_r: str) -> dict:
    """1st-percentile read lengths of each mate (FIGARO needs uniform lengths)."""
    lens_f, lens_r = Counter(), Counter()
    with gzip.open(src_f, "rt") as inf, gzip.open(src_r, "rt") as inr:
        for rec_f, rec_r in zip_longest(_iter_fastq(inf, src_f), _iter_fastq(inr, src_r)):
            if rec_f is None or rec_r is None:
                raise RuntimeError(f"Unpaired read counts in {src_f}")
            lens_f[len(rec_f[1])] += 1
            lens_r[len(rec_r[1])] += 1

    def percentile(counter, q):
        target, seen = int(sum(counter.values()) * q), 0
        for length in sorted(counter):
            seen += counter[length]
            if seen > target:
                return length
        return 0
    return {"pairs": sum(lens_f.values()),
            "p01_f": percentile(lens_f, 0.01), "p01_r": percentile(lens_r, 0.01)}


def _write_uniform_pair(src_f: str, src_r: str, dst_f: str, dst_r: str,
                        fcut: int, rcut: int) -> int:
    """Write uncompressed copies trimmed to fcut/rcut, dropping shorter pairs."""
    kept = 0
    with gzip.open(src_f, "rt") as inf, gzip.open(src_r, "rt") as inr, \
            open(dst_f, "w") as outf, open(dst_r, "w") as outr:
        for rec_f, rec_r in zip(_iter_fastq(inf, src_f), _iter_fastq(inr, src_r)):
            if len(rec_f[1]) < fcut or len(rec_r[1]) < rcut:
                continue
            for out, (h, s, p, q), cut in ((outf, rec_f, fcut), (outr, rec_r, rcut)):
                out.write(f"{h}{s[:cut]}\n{p}{q[:cut]}\n")
            kept += 1
    return kept


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

class AmpliconPipeline:
    """Resumable 16S amplicon pipeline for one BioProject and one mode."""

    def __init__(self, config: dict, mode: str, workdir: Path, threads: int,
                 reference_dir: Optional[Path] = None, run_figaro: bool = False,
                 export_dir: Optional[Path] = None):
        self.cfg = config
        self.mode = mode
        self.mode_cfg = config["modes"][mode]
        self.threads = threads
        self.run_figaro = run_figaro
        self.export_dir = export_dir

        self.root = workdir
        self.raw_fastq = self.root / "raw_fastq"        # shared download cache for all modes
        self.sra_dir = self.root / "sra"
        self.tmp = self.root / "tmp"
        self.meta = self.root / "metadata"
        # The SILVA build is shared by every project under the same parent folder.
        self.ref = reference_dir or (self.root.parent / "reference" / f"silva{config['reference']['silva_version']}")
        self.cohort = self.root / "cohorts" / mode
        self.logs = self.cohort / "logs"
        for d in (self.raw_fastq, self.sra_dir, self.tmp, self.meta, self.ref, self.cohort, self.logs):
            d.mkdir(parents=True, exist_ok=True)

        # DADA2 truncation / maxEE: fixed in the config, or chosen by FIGARO (dada2.truncation: figaro).
        self.layout = config.get("layout", "paired").lower()   # paired | single
        if self.layout not in ("paired", "single"):
            raise ValueError("layout must be 'paired' or 'single'")
        self.single = self.layout == "single"
        d = config["dada2"]
        self.trunc_from_figaro = d.get("truncation", "config") == "figaro"
        if self.single and self.trunc_from_figaro:
            raise ValueError("FIGARO needs paired-end reads; single-end configs use dada2.truncation: config")
        self.trunc = (d.get("trunc_f"), d.get("trunc_r"))
        self.max_ee = (d["max_ee_f"], d["max_ee_r"])
        if self.trunc_from_figaro and not run_figaro:
            raise ValueError("dada2.truncation is 'figaro': pass --figaro")

        self.timings: Dict[str, float] = {}
        self.info: Optional[pd.DataFrame] = None
        self.selected: Optional[pd.DataFrame] = None
        self.selected_runs: List[str] = []
        self.stats: Optional[pd.DataFrame] = None
        self.counts: Optional[pd.DataFrame] = None

    # ----- command helpers --------------------------------------------------

    def checked(self, cmd, log: Path, cwd=None, env=None) -> None:
        cmd = [str(c) for c in cmd]
        logger.info("RUN: %s", " ".join(cmd))
        log.parent.mkdir(parents=True, exist_ok=True)
        start = time.time()
        with log.open("w", encoding="utf-8") as fh:
            proc = subprocess.run(cmd, cwd=cwd, env=env, stdout=fh, stderr=subprocess.STDOUT)
        if proc.returncode:
            tail = log.read_text(errors="replace")[-5000:]
            logger.error("Command failed (exit %s). Log tail:\n%s", proc.returncode, tail)
            raise RuntimeError(f"Command failed ({proc.returncode}): {cmd[:4]}. Full log: {log}")
        logger.info("  done in %.1fs (log: %s)", time.time() - start, log.name)

    def qiime(self, *args, log: Path) -> None:
        self.checked(["qiime", *args], log=log)

    def produce(self, outputs: List[Path], build_args: Callable[[Dict[Path, Path]], list],
                log: Path, label: str) -> None:
        """Run one QIIME command writing to *.partial.qza, then rename all outputs.

        build_args receives {final_path: partial_path} and returns the qiime arguments.
        """
        present = [p.exists() for p in outputs]
        if all(present):
            logger.info("%s: outputs exist, skipping", label)
            return
        if any(present):
            raise RuntimeError(f"Partial {label} outputs {[str(p) for p, ok in zip(outputs, present) if ok]}; "
                               "review and delete them before rerunning.")
        tmp = {p: partial_path(p) for p in outputs}
        for t in tmp.values():
            t.unlink(missing_ok=True)
        self.qiime(*build_args(tmp), log=log)
        for final, t in tmp.items():
            os.replace(t, final)

    # ----- step 1: metadata -------------------------------------------------

    def step_metadata(self) -> None:
        project = self.cfg["project"]
        runinfo = self.meta / "SraRunInfo.csv"
        if not runinfo.exists():
            logger.info("Querying NCBI SRA RunInfo for %s", project)
            headers = {"User-Agent": f"Academic {project} reproducibility (public SRA data)"}
            query = urllib.parse.urlencode({"db": "sra", "term": f"{project}[BioProject]",
                                            "retmode": "json", "retmax": "0", "usehistory": "y"})
            req = urllib.request.Request(
                f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{query}", headers=headers)
            with urllib.request.urlopen(req, timeout=90) as r:
                e = json.load(r)["esearchresult"]
            if int(e["count"]) == 0:
                raise RuntimeError("NCBI SRA search returned no matches")
            query = urllib.parse.urlencode({"query_key": e["querykey"], "WebEnv": e["webenv"]})
            req = urllib.request.Request(
                f"https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?{query}", headers=headers)
            with urllib.request.urlopen(req, timeout=240) as r:
                text = r.read().decode("utf-8")
            df = pd.read_csv(io.StringIO(text), dtype=str)
            if not {"Run", "LibraryName", "LibraryLayout"}.issubset(df.columns):
                raise RuntimeError("NCBI RunInfo lacks required columns; inspect the response")
            df.to_csv(runinfo, index=False)

        lab = self.cfg["labels"]
        exp = self.cfg["expected"]
        info = pd.read_csv(runinfo, dtype=str)
        info = info.loc[info["LibraryLayout"].str.upper() == self.layout.upper()].copy()
        if self.cfg.get("biosample_attributes"):
            # Some studies keep diagnosis / sample type only in BioSample attributes.
            info = info.merge(self.biosample_attributes(info["BioSample"].dropna().unique()),
                              on="BioSample", how="left", suffixes=("", "_biosample"))
        info = apply_filters(info, self.cfg)
        if lab.get("file"):
            info = self.merge_label_file(info, lab)
        info["Diagnosis_provisional"] = assign_labels(info[lab["source_column"]], lab)
        unique_col = lab.get("unique_column", lab["source_column"])  # one run per sample
        n_expected = exp.get("runs", exp.get("paired_runs"))
        if len(info) != n_expected or not info["Run"].is_unique or not info[unique_col].is_unique:
            raise RuntimeError(f"Expected {n_expected} distinct {self.layout} runs with distinct "
                               f"{unique_col}; got {len(info)}. Review RunInfo")
        if info["Diagnosis_provisional"].isna().any():
            raise RuntimeError("Some label prefixes do not match the label regex; verify metadata manually")
        actual = info["Diagnosis_provisional"].value_counts().to_dict()
        if actual != exp["groups"]:
            raise RuntimeError(f"Unexpected provisional groups: {actual}; expected {exp['groups']}")
        info[["Run", lab["source_column"], "Diagnosis_provisional"]].to_csv(
            self.meta / "run_labels_provisional.csv", index=False)
        self.info = info
        logger.info("Runs: %d | provisional groups: %s", len(info), actual)
        logger.warning("Labels are inferred from %s; they are NOT independent clinical verification.",
                       lab["source_column"])

    def merge_label_file(self, info: pd.DataFrame, lab: dict) -> pd.DataFrame:
        """Labels supplied by the user (from the original publication) when NCBI has none.

        labels.file: CSV (path relative to the repository) with the column `labels.join_column`
        (e.g. Run or SampleName) and the column `labels.source_column` holding the diagnosis."""
        path = Path(lab["file"])
        if not path.is_absolute():
            path = Path(__file__).resolve().parents[2] / path
        if not path.exists():
            raise FileNotFoundError(
                f"{self.cfg['project']}: NCBI metadata has no diagnosis labels. Put them in {path} "
                f"(columns: {lab['join_column']}, {lab['source_column']}), taken from the original study.")
        extra = pd.read_csv(path, dtype=str)
        key = lab["join_column"]
        if not extra[key].is_unique:
            raise RuntimeError(f"{path}: {key} values must be unique")
        info = info.drop(columns=[lab["source_column"]], errors="ignore")
        # Inner join: the file also decides which runs are used (e.g. 46 of a project's 80 runs).
        return info.merge(extra[[key, lab["source_column"]]], on=key, how="inner", validate="one_to_one")

    def biosample_attributes(self, accessions) -> pd.DataFrame:
        """BioSample attributes (harmonized names) for the given accessions; cached in metadata/."""
        cache = self.meta / "biosample_attributes.csv"
        if cache.exists():
            return pd.read_csv(cache, dtype=str)
        import xml.etree.ElementTree as ET
        headers = {"User-Agent": f"Academic {self.cfg['project']} reproducibility (public SRA data)"}
        eutils = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        rows, accessions = [], sorted(accessions)
        for i in range(0, len(accessions), 150):
            query = urllib.parse.urlencode({"db": "biosample", "term": " OR ".join(accessions[i:i + 150]),
                                            "retmode": "json", "retmax": "500"})
            with urllib.request.urlopen(urllib.request.Request(f"{eutils}/esearch.fcgi?{query}",
                                                               headers=headers), timeout=120) as r:
                ids = json.load(r)["esearchresult"]["idlist"]
            query = urllib.parse.urlencode({"db": "biosample", "id": ",".join(ids), "retmode": "xml"})
            with urllib.request.urlopen(urllib.request.Request(f"{eutils}/efetch.fcgi?{query}",
                                                               headers=headers), timeout=240) as r:
                root = ET.fromstring(r.read())
            for bs in root.iter("BioSample"):
                row = {"BioSample": bs.get("accession")}
                title = bs.find("Description/Title")
                row["bs_title"] = title.text if title is not None else None   # some studies put the diagnosis here
                for a in bs.iter("Attribute"):
                    row[a.get("harmonized_name") or a.get("attribute_name")] = a.text
                rows.append(row)
            time.sleep(0.4)
        df = pd.DataFrame(rows)
        if set(df["BioSample"]) != set(accessions):
            raise RuntimeError(f"Fetched {len(df)} BioSamples, expected {len(accessions)}")
        df.to_csv(cache, index=False)
        return df

    # ----- step 2: select ---------------------------------------------------

    def step_select(self) -> None:
        info = self.info
        n = self.mode_cfg.get("n_samples")
        if n is None:
            runs = sorted(info["Run"])
        else:
            # Stratified: always include the pilot run, then add from whichever
            # group currently has the fewest picks so every group is represented.
            rng = random.Random(self.cfg["selection_seed"])
            groups = sorted(info["Diagnosis_provisional"].unique())
            pools = {}
            for g in groups:
                pools[g] = sorted(info.loc[info["Diagnosis_provisional"] == g, "Run"])
                rng.shuffle(pools[g])
            pilot = self.cfg.get("pilot_run")
            runs = [pilot] if pilot in set(info["Run"]) else []
            picked = Counter(info.set_index("Run").loc[runs, "Diagnosis_provisional"])
            for g in pools:
                pools[g] = [r for r in pools[g] if r not in runs]
            while len(runs) < n:
                candidates = [g for g in groups if pools[g]]
                if not candidates:
                    raise RuntimeError("Cannot choose enough runs")
                g = min(candidates, key=lambda x: (picked[x], groups.index(x)))
                runs.append(pools[g].pop())
                picked[g] += 1

        self.selected_runs = runs
        cols = ["Run", self.cfg["labels"]["source_column"], "Diagnosis_provisional"]
        self.selected = info.set_index("Run").loc[runs].reset_index()
        logger.info("Selected %d runs (%s): %s", len(runs), self.mode,
                    self.selected["Diagnosis_provisional"].value_counts().to_dict())

        # Joint denoising assumes one comparable sequencing protocol.
        for column in ("Platform", "LibraryStrategy", "LibrarySource"):
            if column in self.selected:
                values = sorted(self.selected[column].dropna().astype(str).unique())
                logger.info("%s across cohort: %s", column, values)
                if len(values) > 1:
                    raise RuntimeError(f"Mixed {column}: split by sequencing protocol before joint DADA2")

        # Refuse to reuse outputs produced with a different recipe.
        recipe = {"mode": self.mode, "mode_cfg": self.mode_cfg, "selected_runs": runs,
                  "primers": self.cfg["primers"], "cutadapt": self.cfg["cutadapt"],
                  "dada2": self.cfg["dada2"], "reference": self.cfg["reference"],
                  "vsearch": self.cfg["vsearch"]}
        # Newer options join the recipe only when set, so existing cohorts keep their fingerprint.
        recipe.update({k: self.cfg[k] for k in ("layout", "reorient") if k in self.cfg})
        fingerprint = hashlib.sha256(json.dumps(recipe, sort_keys=True).encode()).hexdigest()
        recipe_file = self.cohort / "analysis_recipe.json"
        if recipe_file.exists():
            old = json.loads(recipe_file.read_text())
            if old["fingerprint"] != fingerprint:
                raise RuntimeError(f"{self.cohort} holds results from a different sample selection or "
                                   "parameters. Move that folder aside before changing the recipe.")
        else:
            recipe_file.write_text(json.dumps({"fingerprint": fingerprint, "recipe": recipe}, indent=2))
        logger.info("Recipe fingerprint: %s", fingerprint[:16])
        self.selected[cols].to_csv(self.cohort / "selected_runs.tsv", sep="\t", index=False)

    # ----- step 3: download -------------------------------------------------

    def mates(self) -> tuple:
        return (1,) if self.single else (1, 2)

    def raw_pair(self, run: str) -> List[Path]:
        """Raw FASTQs of a run: [R1, R2], or [R1] for single-end."""
        return [self.raw_fastq / f"{run}_{m}.fastq.gz" for m in self.mates()]

    def source_pair(self, run: str) -> List[Path]:
        """Raw FASTQs, or the subsampled copies in dry-run mode."""
        if self.mode_cfg.get("subsample_pairs"):
            return [self.cohort / "fastq_subsampled" / f"{run}_{m}.fastq.gz" for m in self.mates()]
        return self.raw_pair(run)

    def input_pair(self, run: str) -> List[Path]:
        """FASTQs fed to QIIME / FIGARO: source_pair, re-oriented first if `reorient: true`."""
        if self.cfg.get("reorient"):
            return [self.cohort / "fastq_oriented" / f"{run}_{m}.fastq.gz" for m in self.mates()]
        return self.source_pair(run)

    def _download_one(self, run: str, threads: int) -> str:
        paths = self.raw_pair(run)
        if all(nonempty(p, 50) for p in paths):
            return f"cached {run}"
        if shutil.disk_usage(self.root).free < 6e9:
            raise RuntimeError("<6 GB free: make space before another fasterq-dump")
        stage, tmp = self.raw_fastq / f".stage_{run}", self.tmp / run
        shutil.rmtree(stage, ignore_errors=True)
        stage.mkdir(parents=True)
        tmp.mkdir(parents=True, exist_ok=True)
        self.checked(["prefetch", run, "--max-size", self.cfg["download"]["max_size"],
                      "--output-directory", self.sra_dir], log=self.logs / f"{run}.prefetch.log")
        sra = self.sra_dir / run / f"{run}.sra"
        self.checked(["fasterq-dump", "--split-files", "-e", threads, "--outdir", stage,
                      "--temp", tmp, sra if sra.exists() else run],
                     log=self.logs / f"{run}.fasterq.log")
        if self.single and not (stage / f"{run}_1.fastq").exists() and (stage / f"{run}.fastq").exists():
            os.replace(stage / f"{run}.fastq", stage / f"{run}_1.fastq")  # fasterq-dump names single reads <run>.fastq
        for m in self.mates():
            raw = stage / f"{run}_{m}.fastq"
            if not nonempty(raw):
                raise RuntimeError(f"Missing mate-{m} FASTQ for {run}")
            self.checked(["pigz", "-p", threads, "-1", "-f", raw], log=self.logs / f"{run}.pigz{m}.log")
        staged = [stage / p.name for p in paths]
        self.checked(["gzip", "-t", *staged], log=self.logs / f"{run}.gzipcheck.log")
        for s, p in zip(staged, paths):
            os.replace(s, p)
        shutil.rmtree(stage, ignore_errors=True)
        shutil.rmtree(tmp, ignore_errors=True)
        if self.cfg["download"].get("delete_sra_after_fastq", True):
            shutil.rmtree(self.sra_dir / run, ignore_errors=True)
        return f"downloaded {run}"

    def step_download(self) -> None:
        runs = self.selected_runs
        jobs = max(1, min(self.cfg["download"]["parallel_jobs"], len(runs)))
        per_job = max(1, self.threads // jobs)
        logger.info("Downloading %d runs: %d parallel jobs x %d threads", len(runs), jobs, per_job)
        done = 0
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            futures = {pool.submit(self._download_one, r, per_job): r for r in runs}
            for fut in as_completed(futures):
                done += 1
                logger.info("[%d/%d] %s", done, len(runs), fut.result())

        sub = self.mode_cfg.get("subsample_pairs")
        if sub:
            (self.cohort / "fastq_subsampled").mkdir(exist_ok=True)
            todo = [r for r in runs if not all(nonempty(p, 50) for p in self.source_pair(r))]
            with ProcessPoolExecutor(max_workers=min(self.threads, max(1, len(todo)))) as pool:
                futures = {pool.submit(_subsample, list(map(str, self.raw_pair(r))),
                                       list(map(str, self.source_pair(r))), sub): r for r in todo}
                for fut in as_completed(futures):
                    logger.info("Subsampled %s to %d pairs", futures[fut], fut.result())
        if self.cfg.get("reorient"):
            if self.single or not self.cfg.get("primers"):
                raise ValueError("reorient needs paired-end reads and both primers")
            p = self.cfg["primers"]
            (self.cohort / "fastq_oriented").mkdir(exist_ok=True)
            todo = [r for r in runs if not all(nonempty(x, 50) for x in self.input_pair(r))]
            swapped = []
            with ProcessPoolExecutor(max_workers=min(self.threads, max(1, len(todo)))) as pool:
                futures = {pool.submit(_orient_pair, *map(str, self.source_pair(r)), *map(str, self.input_pair(r)),
                                       p["r1"], p["r2"], int(p.get("max_offset", 12))): r for r in todo}
                for fut in as_completed(futures):
                    n = fut.result()
                    swapped.append(n["swapped"] / max(n["pairs"], 1))
            if swapped:
                logger.info("Re-oriented %d runs: median %.1f%% of pairs had the primers swapped",
                            len(swapped), 100 * float(np.median(swapped)))
        logger.info("All %d FASTQ pairs ready", len(runs))

    # ----- step 4: import ---------------------------------------------------

    def step_import(self) -> None:
        manifest = self.cohort / "manifest.tsv"
        runs = self.selected_runs
        if self.single:
            cols = {"absolute-filepath": [str(self.input_pair(r)[0].resolve()) for r in runs]}
            qtype, qformat = "SampleData[SequencesWithQuality]", "SingleEndFastqManifestPhred33V2"
        else:
            cols = {"forward-absolute-filepath": [str(self.input_pair(r)[0].resolve()) for r in runs],
                    "reverse-absolute-filepath": [str(self.input_pair(r)[1].resolve()) for r in runs]}
            qtype, qformat = "SampleData[PairedEndSequencesWithQuality]", "PairedEndFastqManifestPhred33V2"
        pd.DataFrame({"sample-id": runs, **cols}).to_csv(manifest, sep="\t", index=False)
        demux = self.cohort / "demux-paired.qza"   # name kept for both layouts (resumability)
        self.produce([demux], lambda o: [
            "tools", "import", "--type", qtype,
            "--input-path", manifest, "--input-format", qformat,
            "--output-path", o[demux]], log=self.logs / "qiime_import.log", label="import")
        summary = self.cohort / "demux-summary.qzv"
        self.produce([summary], lambda o: [
            "demux", "summarize", "--i-data", demux, "--o-visualization", o[summary]],
            log=self.logs / "demux_summarize.log", label="demux summary")

    # ----- step 5: cutadapt -------------------------------------------------

    def step_cutadapt(self) -> None:
        p, c = self.cfg["primers"], self.cfg["cutadapt"]
        demux = self.cohort / "demux-paired.qza"
        trimmed, stats = self.cohort / "demux-primer-trimmed.qza", self.cohort / "cutadapt-stats.qza"
        if not p:
            # Primers are not in the reads (e.g. EMP-style V4, where sequencing starts after
            # the PCR primer): nothing to trim, so the imported reads go straight to QC/DADA2.
            if not trimmed.exists():
                atomic_copy(demux, trimmed)
            logger.info("primers: null -> Cutadapt skipped; reads used as imported")
            return
        if self.single:
            fronts = ["--p-front", p["r1"]] if p.get("r1") else []   # r1: null + r1_3prime = 3' trim only
        else:
            fronts = ["--p-front-f", p["r1"], "--p-front-r", p["r2"]]
        if p.get("r1_3prime"):
            # Reads that span the whole amplicon (e.g. pre-merged pairs) end in the reverse-complemented
            # reverse primer; remove it too (times 2 = front primer and 3' primer in one pass).
            fronts += ["--p-adapter" if self.single else "--p-adapter-f", p["r1_3prime"], "--p-times", 2]
        self.produce([trimmed, stats], lambda o: [
            "cutadapt", "trim-single" if self.single else "trim-paired", "--p-cores", self.threads,
            "--i-demultiplexed-sequences", demux, *fronts,
            "--p-error-rate", c["error_rate"], "--p-overlap", c["min_overlap"],
            "--p-discard-untrimmed",
            "--o-trimmed-sequences", o[trimmed], "--o-stats", o[stats]],
            log=self.logs / "cutadapt.log", label="Cutadapt")
        summary = self.cohort / "trimmed-summary.qzv"
        self.produce([summary], lambda o: [
            "demux", "summarize", "--i-data", trimmed, "--o-visualization", o[summary]],
            log=self.logs / "trimmed_summarize.log", label="trimmed summary")

    # ----- step 6: qc gate --------------------------------------------------

    def _qza_stats(self, qza: Path, fmin: int = 0, rmin: int = 0) -> pd.DataFrame:
        mates: Dict[str, Dict[int, str]] = {}
        with zipfile.ZipFile(qza) as zf:
            for member in zf.namelist():
                m = QIIME_FASTQ_RE.match(Path(member).name)
                if m:
                    mates.setdefault(m["sid"], {})[int(m["mate"])] = member
        if set(mates) != set(self.selected_runs):
            raise RuntimeError(f"{qza.name} sample IDs differ: expected {len(self.selected_runs)}, "
                               f"got {len(mates)}")
        results = {}
        with ProcessPoolExecutor(max_workers=min(self.threads, len(mates))) as pool:
            futures = {pool.submit(_qza_pair_stats, str(qza), mates[r][1], mates[r].get(2), fmin, rmin): r
                       for r in self.selected_runs}
            for fut in as_completed(futures):
                results[futures[fut]] = fut.result()
        return pd.DataFrame.from_dict(results, orient="index").rename_axis("Run").reset_index()

    def step_qc(self) -> None:
        report = self.cohort / "pre_dada2_qc.csv"
        d, q = self.cfg["dada2"], self.cfg["qc"]
        if report.exists():
            qc = pd.read_csv(report)
        else:
            fmin, rmin = self.trunc
            trimmed = self._qza_stats(self.cohort / "demux-primer-trimmed.qza", fmin or 0,
                                      0 if self.single else (rmin or 0))
            raw = self._qza_stats(self.cohort / "demux-paired.qza")
            qc = trimmed.merge(raw[["Run", "pairs"]].rename(columns={"pairs": "raw_pairs"}), on="Run")
            qc["primer_pair_retention"] = qc["pairs"] / qc["raw_pairs"]
            qc["eligible_for_truncation"] = qc["eligible_fraction"]
            qc["qc_pass"] = ((qc["primer_pair_retention"] >= q["min_cutadapt_retained"]) &
                             (qc["eligible_for_truncation"] >= q["min_trunc_eligible"]))
            qc = qc.merge(self.selected[["Run", "Diagnosis_provisional"]], on="Run", validate="one_to_one")
            qc.to_csv(report, index=False)
        logger.info("Pre-DADA2 QC (first rows):\n%s", qc.head(15).to_string())
        failed = qc.loc[~qc["qc_pass"].astype(bool)]
        if len(failed):
            report.unlink()  # recompute after the cause is fixed
            raise RuntimeError(f"{len(failed)} samples fail primer/truncation QC:\n{failed.to_string()}\n"
                               "Investigate primers/orientation/protocol; do not just lower thresholds.")
        logger.info("QC PASS for all %d runs (does not verify primer identity or diagnoses)", len(qc))

    # ----- step 7: figaro (optional) ---------------------------------------

    def step_figaro(self) -> None:
        if not self.run_figaro:
            logger.info("FIGARO skipped (pass --figaro to run it)")
            return
        if self.single:
            logger.info("FIGARO skipped: it models paired-end overlap; single-end truncation comes "
                        "from the config (dada2.trunc_f = %s)", self.trunc[0])
            return
        f = self.cfg["figaro"]
        base = self.cohort / "figaro"
        fig_in, fig_out = base / "uniform_input", base / "report"
        report = fig_out / "trimParameters.json"
        if not report.exists():
            home = Path(os.environ.get("FIGARO_HOME", Path.home() / "opt" / "figaro"))
            py = Path(os.environ.get("FIGARO_PYTHON", Path.home() / "miniforge3" / "envs" / "figaro" / "bin" / "python"))
            script = home / "figaro" / "figaro.py"
            if not script.exists() or not py.exists():
                raise RuntimeError(f"FIGARO not found ({script}, {py}). Run envs/setup_local.sh "
                                   "or set FIGARO_HOME / FIGARO_PYTHON.")
            runs = self.selected_runs
            workers = min(self.threads, len(runs))
            with ProcessPoolExecutor(max_workers=workers) as pool:
                profiles = dict(zip(runs, pool.map(_length_profile,
                                                   *zip(*[map(str, self.input_pair(r)) for r in runs]))))
            # FIGARO needs one read length per direction across the whole cohort.
            fcut = min(p["p01_f"] for p in profiles.values())
            rcut = min(p["p01_r"] for p in profiles.values())
            shutil.rmtree(fig_in, ignore_errors=True)
            fig_in.mkdir(parents=True)
            fig_out.mkdir(parents=True, exist_ok=True)
            dsts = [[str(fig_in / f"{r}_S{i}_L001_R{m}_001.fastq") for m in (1, 2)]
                    for i, r in enumerate(runs, 1)]
            with ProcessPoolExecutor(max_workers=workers) as pool:
                kept = list(pool.map(_write_uniform_pair,
                                     [str(self.input_pair(r)[0]) for r in runs],
                                     [str(self.input_pair(r)[1]) for r in runs],
                                     [d[0] for d in dsts], [d[1] for d in dsts],
                                     [fcut] * len(runs), [rcut] * len(runs)))
            retention = {r: k / profiles[r]["pairs"] for r, k in zip(runs, kept)}
            low = {r: round(v, 3) for r, v in retention.items() if v < 0.95}
            logger.info("FIGARO uniform lengths %d/%d; median retention %.3f", fcut, rcut,
                        float(np.median(list(retention.values()))))
            if low:
                logger.warning("Runs losing >5%% of pairs to FIGARO length normalization: %s", low)
            self.checked([py, script, "-i", fig_in, "-o", fig_out,
                          "-a", f["amplicon_length"], "-f", f["forward_primer_length"],
                          "-r", f["reverse_primer_length"], "-m", f["min_overlap"], "-F", "illumina"],
                         log=base / "figaro.log", cwd=home, env={**os.environ, "MPLBACKEND": "Agg"})
            if not report.exists():
                raise FileNotFoundError(report)
            (base / "normalization.json").write_text(json.dumps(
                {"fcut": fcut, "rcut": rcut, "retention": retention}, indent=2))
            shutil.rmtree(fig_in, ignore_errors=True)  # uncompressed copies are large
        top = pd.json_normalize(json.loads(report.read_text())[:10])
        top.to_csv(base / "figaro_top10.csv", index=False)
        logger.info("Top FIGARO suggestions (RAW read coordinates):\n%s", top.to_string())
        if not self.trunc_from_figaro:
            logger.info("DADA2 uses config trunc %s/%s (post-Cutadapt coordinates); FIGARO is advisory.",
                        *self.trunc)
            return
        # Use FIGARO's best-scoring choice, as in the paper. FIGARO reports positions in raw reads
        # (primer included); DADA2 here runs on primer-trimmed reads, so subtract the primer length
        # (incl. any spacer) given to FIGARO. Its maxExpectedError suggestion is used as-is.
        best = json.loads(report.read_text())[0]
        tf, tr = (int(x) for x in best["trimPosition"])
        ef, er = (int(x) for x in best["maxExpectedError"])
        self.trunc = (tf - int(f["forward_primer_length"]), tr - int(f["reverse_primer_length"]))
        self.max_ee = (ef, er)
        need = self.cfg["expected_insert_upper"] + self.cfg["min_paired_overlap"]
        if sum(self.trunc) < need:
            raise RuntimeError(f"FIGARO truncation {self.trunc} leaves < {self.cfg['min_paired_overlap']} bp "
                               f"overlap for a {self.cfg['expected_insert_upper']} bp insert")
        (base / "dada2_params.json").write_text(json.dumps(
            {"figaro_best": best, "trunc_len_f": self.trunc[0], "trunc_len_r": self.trunc[1],
             "max_ee_f": ef, "max_ee_r": er, "coordinates": "post-Cutadapt (raw minus primer length)"},
            indent=2))
        logger.info("DADA2 will use FIGARO's best: trunc %d/%d, maxEE %d/%d (raw positions %d/%d, "
                    "predicted retention %s%%)", *self.trunc, ef, er, tf, tr, best.get("readRetentionPercent"))

    # ----- step 8: dada2 ----------------------------------------------------

    def step_dada2(self) -> None:
        d = self.cfg["dada2"]
        dd = self.cohort / "dada2"
        dd.mkdir(exist_ok=True)
        outs = {k: dd / f"{k}.qza" for k in ("table", "rep-seqs", "stats", "transitions")}
        n_learn = self.mode_cfg.get("n_reads_learn", d["n_reads_learn"])
        if self.single:   # trunc 0 = keep full length (e.g. reads already merged by the submitters)
            mode = ["denoise-single", "--p-trim-left", d.get("trim_left", 0), "--p-trunc-len", self.trunc[0] or 0,
                    "--p-max-ee", self.max_ee[0]]
        else:
            mode = ["denoise-paired", "--p-trim-left-f", 0, "--p-trim-left-r", 0,
                    "--p-trunc-len-f", self.trunc[0], "--p-trunc-len-r", self.trunc[1],
                    "--p-max-ee-f", self.max_ee[0], "--p-max-ee-r", self.max_ee[1]]
        self.produce(list(outs.values()), lambda o: [
            "dada2", *mode[:1], "--i-demultiplexed-seqs", self.cohort / "demux-primer-trimmed.qza", *mode[1:],
            "--p-n-reads-learn", n_learn, "--p-n-threads", self.threads,
            "--o-table", o[outs["table"]], "--o-representative-sequences", o[outs["rep-seqs"]],
            "--o-denoising-stats", o[outs["stats"]], "--o-base-transition-stats", o[outs["transitions"]]],
            log=self.logs / "dada2_joint.log", label="joint DADA2")

    # ----- step 9: dada2 stats ---------------------------------------------

    def step_dada2_stats(self) -> None:
        export = self.cohort / "dada2" / "exported_stats"
        tsv = export / "stats.tsv"
        if not tsv.exists():
            self.qiime("tools", "export", "--input-path", self.cohort / "dada2" / "stats.qza",
                       "--output-path", export, log=self.logs / "export_dada2_stats.log")
        lines = tsv.read_text(encoding="utf-8").splitlines()
        header = next((i for i, l in enumerate(lines) if "\tinput\t" in l), None)
        if header is None:
            raise ValueError(f"Cannot identify DADA2 stats header in {tsv}")
        body = [l for l in lines[header + 1:] if l.strip() and not l.lstrip().startswith("#")]
        stats = pd.read_csv(io.StringIO("\n".join([lines[header], *body])), sep="\t")
        stats = stats.rename(columns={stats.columns[0]: "sample-id"})
        if "merged" not in stats:          # single-end: nothing to merge
            stats["merged"] = stats["denoised"]
        for col in ("input", "filtered", "denoised", "merged", "non-chimeric"):
            stats[col] = pd.to_numeric(stats[col], errors="raise")
        for name, col in (("retained_pct", "non-chimeric"), ("merged_pct", "merged"), ("filtered_pct", "filtered")):
            stats[name] = (100 * stats[col] / stats["input"].replace(0, np.nan)).round(2)
        if len(stats) != len(self.selected_runs) or set(stats["sample-id"]) != set(self.selected_runs):
            raise RuntimeError("DADA2 stats sample IDs do not match the selected runs")
        if (stats["non-chimeric"] <= 0).any():
            raise RuntimeError("A sample has zero non-chimeric reads; investigate before taxonomy")
        stats.to_csv(self.cohort / "dada2_per_sample_stats.csv", index=False)
        self.stats = stats
        logger.info("Lowest retention:\n%s", stats.nsmallest(min(10, len(stats)), "retained_pct")[
            ["sample-id", "input", "non-chimeric", "merged_pct", "retained_pct"]].to_string())
        logger.info("Input reads %d -> non-chimeric %d | retention median %.2f%% (min %.2f%%, max %.2f%%)",
                    stats["input"].sum(), stats["non-chimeric"].sum(), stats["retained_pct"].median(),
                    stats["retained_pct"].min(), stats["retained_pct"].max())
        low = stats[stats["retained_pct"] < self.cfg["qc"]["min_dada2_retained_pct"]]
        if len(low):
            raise RuntimeError(f"Very low DADA2 retention:\n{low[['sample-id', 'retained_pct']].to_string()}")

    # ----- step 10: reference ----------------------------------------------

    def step_reference(self) -> None:
        r = self.cfg["reference"]
        rna, tax, dna = (self.ref / f"silva-{r['silva_version']}-{k}.qza" for k in ("rna", "taxonomy", "dna"))
        self.produce([rna, tax], lambda o: [
            "rescript", "get-silva-data", "--p-version", r["silva_version"], "--p-target", r["silva_target"],
            "--o-silva-sequences", o[rna], "--o-silva-taxonomy", o[tax]],
            log=self.ref / "silva_download.log", label="SILVA download")
        self.produce([dna], lambda o: [
            "rescript", "reverse-transcribe", "--i-rna-sequences", rna, "--o-dna-sequences", o[dna]],
            log=self.ref / "silva_reverse_transcribe.log", label="SILVA reverse-transcribe")
        self.silva_dna, self.silva_tax = dna, tax

    # ----- step 11: taxonomy -----------------------------------------------

    def step_taxonomy(self) -> None:
        v = self.cfg["vsearch"]
        td = self.cohort / "taxonomy"
        td.mkdir(exist_ok=True)
        tax_qza, hits = td / "taxonomy.qza", td / "vsearch-hits.qza"
        self.produce([tax_qza, hits], lambda o: [
            "feature-classifier", "classify-consensus-vsearch",
            "--i-query", self.cohort / "dada2" / "rep-seqs.qza",
            "--i-reference-reads", self.silva_dna, "--i-reference-taxonomy", self.silva_tax,
            "--p-perc-identity", v["identity"], "--p-query-cov", v["query_coverage"],
            "--p-maxaccepts", v["maxaccepts"], "--p-strand", "both", "--p-threads", self.threads,
            "--o-classification", o[tax_qza], "--o-search-results", o[hits]],
            log=self.logs / "vsearch_taxonomy.log", label="VSEARCH taxonomy")
        tsv = td / "exported" / "taxonomy.tsv"
        if not tsv.exists():
            self.qiime("tools", "export", "--input-path", tax_qza, "--output-path", td / "exported",
                       log=self.logs / "export_taxonomy.log")
        taxonomy = pd.read_csv(tsv, sep="\t")
        taxonomy["Genus"] = taxonomy["Taxon"].str.extract(GENUS_RE)[0]
        taxonomy.to_csv(self.cohort / "asv_taxonomy_consensus.csv", index=False)
        logger.info("ASVs: %d | genus labels: %d", len(taxonomy), taxonomy["Genus"].nunique(dropna=True))
        if "Consensus" in taxonomy:
            low = pd.to_numeric(taxonomy["Consensus"], errors="coerce") < 0.8
            logger.info("Lower-consensus (<0.8) assignments to inspect: %d", int(low.sum()))

    # ----- step 12: genus ---------------------------------------------------

    def step_genus(self) -> None:
        gd = self.cohort / "genus"
        gd.mkdir(exist_ok=True)
        genus_qza = gd / "genus-table.qza"
        self.produce([genus_qza], lambda o: [
            "taxa", "collapse", "--i-table", self.cohort / "dada2" / "table.qza",
            "--i-taxonomy", self.cohort / "taxonomy" / "taxonomy.qza",
            "--p-level", 6, "--o-collapsed-table", o[genus_qza]],
            log=self.logs / "collapse_genus.log", label="genus collapse")
        biom = gd / "biom_export" / "feature-table.biom"
        if not biom.exists():
            self.qiime("tools", "export", "--input-path", genus_qza, "--output-path", biom.parent,
                       log=self.logs / "genus_biom_export.log")
        biom_tsv = gd / "genus-lineages.tsv"
        if not biom_tsv.exists():
            self.checked(["biom", "convert", "-i", biom, "-o", biom_tsv, "--to-tsv"],
                         log=self.logs / "biom_to_tsv.log")

        runs = self.selected_runs
        lineages = pd.read_csv(biom_tsv, sep="\t", skiprows=1)
        lineages["Genus"] = lineages[lineages.columns[0]].astype(str).str.extract(GENUS_RE)[0].fillna("Unassigned")
        missing = set(runs) - set(lineages.columns)
        if missing:
            raise RuntimeError(f"Missing sample columns in genus BIOM output: {missing}")
        for col in runs:
            lineages[col] = pd.to_numeric(lineages[col], errors="raise")
        counts = lineages.groupby("Genus", sort=True)[runs].sum().T
        counts.index.name = "sample-id"
        counts = counts.round().astype("int64")
        if (counts.sum(axis=1) <= 0).any():
            raise RuntimeError("Some sample has zero genus counts")

        relative = counts.div(counts.sum(axis=1), axis=0)  # proportions include Unassigned
        labels = self.selected.set_index("Run")["Diagnosis_provisional"].rename_axis("sample-id")
        counts_lab = counts.join(labels, how="left", validate="one_to_one")
        relative_lab = relative.join(labels, how="left", validate="one_to_one")
        if relative_lab["Diagnosis_provisional"].isna().any():
            raise RuntimeError("Missing inferred labels; do not guess participant diagnoses")

        # Paper-inspired filters applied to THIS project only (not the paper's 15-study pooling).
        g = self.cfg["genus_filter"]
        present = (counts > 0).sum(axis=0)
        corpus_fraction = counts.sum(axis=0) / counts.to_numpy().sum()
        keep = (present >= g["min_samples_present"]) & (corpus_fraction >= g["min_corpus_fraction"])
        tables = {
            "genus_counts.csv": counts,
            "genus_relative_abundance.csv": relative,
            "genus_counts_with_provisional_labels.csv": counts_lab,
            "genus_relative_with_provisional_labels.csv": relative_lab,
            "genus_counts_filtered_single_project.csv": counts.loc[:, keep],
            "genus_relative_filtered_single_project.csv": relative.loc[:, keep],  # not renormalized
            "genus_filter_audit.csv": pd.DataFrame({
                "samples_present": present, "corpus_fraction": corpus_fraction,
                "kept_in_single_project_filter": keep}).rename_axis("Genus"),
        }
        for name, table in tables.items():
            table.to_csv(gd / name)
        self.counts = counts
        logger.info("Samples: %d | named genera: %d | filtered features: %d",
                    counts.shape[0], sum(c != "Unassigned" for c in counts.columns), int(keep.sum()))

    # ----- step 13: finalize ------------------------------------------------

    def check_expected(self) -> None:
        """Compare with the reference numbers in the config (warn, do not fail)."""
        exp = self.cfg.get("expected_results", {}).get(self.mode)
        if not exp:
            return
        nonchim = self.stats.set_index("sample-id")["non-chimeric"].astype(int)
        n_asvs = len(pd.read_csv(self.cohort / "asv_taxonomy_consensus.csv"))
        observed = {"asvs": n_asvs, "total_input": int(self.stats["input"].sum()),
                    "total_non_chimeric": int(nonchim.sum()),
                    "non_chimeric": {k: int(v) for k, v in nonchim.items()}}
        mismatches = {k: {"expected": v, "observed": observed[k]}
                      for k, v in exp.items() if observed.get(k) != v}
        (self.cohort / "expected_results_check.json").write_text(
            json.dumps({"match": not mismatches, "mismatches": mismatches}, indent=2))
        if mismatches:
            logger.warning("Results differ from the reference numbers in the config: %s. Check tool "
                           "versions (QIIME 2 2026.7) and input data before trusting them.", mismatches)
        else:
            logger.info("REPRODUCIBILITY CHECK PASSED: identical to the reference %s results.", self.mode)

    def step_finalize(self) -> None:
        counts, stats = self.counts, self.stats
        by_sample = counts.sum(axis=1).sort_index()
        dada = stats.set_index("sample-id")["non-chimeric"].astype("int64").sort_index()
        if not by_sample.equals(dada):
            diff = pd.concat([by_sample.rename("genus_total"), dada.rename("dada2_nonchimeric")], axis=1)
            raise RuntimeError(f"Genus totals != DADA2 non-chimeric reads:\n"
                               f"{diff[diff.iloc[:, 0] != diff.iloc[:, 1]].to_string()}")
        if set(counts.index) != set(self.selected_runs):
            raise RuntimeError("Final matrix rows do not match selected runs")
        if not np.allclose(counts.div(counts.sum(axis=1), axis=0).sum(axis=1), 1.0):
            raise RuntimeError("Relative abundances do not sum to 1")

        self.check_expected()
        qiime_version = subprocess.run(["qiime", "--version"], capture_output=True, text=True).stdout.strip()
        provenance = {
            "project": self.cfg["project"], "dataset": dataset_name(self.cfg), "mode": self.mode, "mode_settings": self.mode_cfg,
            "samples": self.selected_runs, "joint_dada2_across_samples": True,
            "qiime_version": qiime_version, "cpu_threads": self.threads,
            "primers": self.cfg["primers"], "cutadapt": self.cfg["cutadapt"], "dada2": self.cfg["dada2"],
            "reference": self.cfg["reference"], "vsearch": self.cfg["vsearch"],
            "figaro_run": self.run_figaro,
            "dada2_used": {"trunc_len": list(self.trunc), "max_ee": list(self.max_ee),
                           "source": "FIGARO best" if self.trunc_from_figaro else "config"},
            "clinical_labels": f"INFERRED FROM NCBI {self.cfg['labels']['source_column']}; NOT VERIFIED",
            "fastqc": "skipped (assumed done by the original study); QIIME quality plots available",
            "warnings": self.cfg.get("warnings", []),
            "step_seconds": self.timings,
            "finished": datetime.now().isoformat(timespec="seconds"),
        }
        (self.cohort / "run_provenance.json").write_text(json.dumps(provenance, indent=2))
        logger.info("QC PASSED: genus totals reconcile with DADA2 for all %d samples", len(counts))

        if self.export_dir:
            dest = self.export_dir / dataset_name(self.cfg) / self.mode
            files = [*(self.cohort / "genus").glob("*.csv"),
                     *(self.cohort / n for n in ("selected_runs.tsv", "manifest.tsv", "pre_dada2_qc.csv",
                                                 "dada2_per_sample_stats.csv", "asv_taxonomy_consensus.csv",
                                                 "run_provenance.json", "step_timings.json",
                                                 "expected_results_check.json")),
                     *self.cohort.glob("*.qzv"), *(self.cohort / "figaro").glob("*.csv"),
                     *(self.cohort / "figaro" / "report").glob("*")]
            for f in files:
                if f.is_file() and f.stat().st_size > 0:  # skip e.g. empty FIGARO logs
                    atomic_copy(f, dest / f.relative_to(self.cohort))
            logger.info("Exported final outputs to %s", dest)

    # ----- driver -------------------------------------------------------------

    def run(self, until: str) -> None:
        stop = STEPS.index(until)
        timings_file = self.cohort / "step_timings.json"
        # Keep the longest time seen per step, so a resumed run (where finished steps
        # are skipped in ~0 s) does not overwrite the real cost of the first run.
        previous = json.loads(timings_file.read_text()) if timings_file.exists() else {}
        for name in STEPS[:stop + 1]:
            logger.info("=" * 70)
            logger.info("STEP %s", name)
            start = time.time()
            getattr(self, "step_" + name.replace("-", "_"))()
            self.timings[name] = max(round(time.time() - start, 1), previous.get(name, 0.0))
            timings_file.write_text(json.dumps(self.timings, indent=2))
        logger.info("=" * 70)
        logger.info("Step timings (s): %s", self.timings)
        logger.info("Total: %.1f min", sum(self.timings.values()) / 60)


# ---------------------------------------------------------------------------
# Preflight
# ---------------------------------------------------------------------------

def preflight(args, cfg: dict) -> bool:
    ok = True
    logger.info("Platform: %s | Python %s", sys.platform, sys.version.split()[0])
    if not sys.platform.startswith(("linux", "darwin")):
        logger.error("QIIME 2 needs Linux or macOS. On Windows run this inside WSL (see docs/16S_PREPROCESSING.md).")
        ok = False
    tools = ["qiime", "biom", "prefetch", "fasterq-dump", "pigz", "gzip"]
    for t in tools:
        path = shutil.which(t)
        logger.info("  %-13s %s", t, path or "MISSING")
        ok &= path is not None
    if shutil.which("qiime"):
        out = subprocess.run(["qiime", "info"], capture_output=True, text=True).stdout
        for plugin in ("cutadapt", "dada2", "feature-classifier", "rescript", "taxa", "demux"):
            found = re.search(rf"^{plugin}:", out, re.M) is not None
            logger.info("  plugin %-19s %s", plugin, "ok" if found else "MISSING")
            ok &= found
    if args.figaro:
        home = Path(os.environ.get("FIGARO_HOME", Path.home() / "opt" / "figaro"))
        found = (home / "figaro" / "figaro.py").exists()
        logger.info("  FIGARO         %s", home if found else "MISSING")
        ok &= found
    mem = available_memory_gb()
    args.workdir.mkdir(parents=True, exist_ok=True)
    free = shutil.disk_usage(args.workdir).free / 1e9
    logger.info("CPUs visible: %d | threads used: %d | RAM: %s GB | free disk at workdir: %.0f GB",
                available_cpus(), args.threads, f"{mem:.1f}" if mem else "?", free)
    if mem and mem < 12:
        logger.warning("Less than 12 GB RAM visible; VSEARCH vs full SILVA may run out of memory. "
                       "Raise 'memory' in %UserProfile%\\.wslconfig.")
    if str(args.workdir.resolve()).startswith("/mnt/"):
        logger.warning("workdir is on a Windows drive (/mnt/...). This works but I/O is much slower; "
                       "prefer a Linux path like ~/biomefm_runs and use --export-dir for Windows copies.")
    mc = cfg["modes"][args.mode]
    logger.info("Plan: mode=%s | samples=%s | subsample=%s | steps up to '%s' | FIGARO=%s",
                args.mode, mc.get("n_samples", "all"), mc.get("subsample_pairs", "none"),
                args.until, args.figaro)
    logger.info("Preflight %s", "PASSED" if ok else "FAILED")
    return ok


def load_config(path: Path) -> dict:
    cfg = yaml.safe_load(path.read_text())
    d = cfg["dada2"]
    if cfg.get("layout", "paired") == "paired" and d.get("truncation", "config") == "config" and \
            d["trunc_f"] + d["trunc_r"] < cfg["expected_insert_upper"] + cfg["min_paired_overlap"]:
        raise ValueError("Truncation lengths do not leave the planned read-pair overlap")
    return cfg


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Local 16S amplicon preprocessing (QIIME 2 + DADA2 + VSEARCH/SILVA)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("Examples:")[1])
    parser.add_argument("--config", type=Path, required=True, help="Project YAML (e.g. configs/PRJNA496408.yaml)")
    parser.add_argument("--mode", choices=MODES, default="dry-run",
                        help="dry-run: 3 samples x 5000 pairs; validation: 10 samples; full: all")
    parser.add_argument("--workdir", type=Path, default=None,
                        help="Working directory (default: ~/biomefm_runs/<project>)")
    parser.add_argument("--reference-dir", type=Path, default=None,
                        help="Where the SILVA artifacts are cached (default: <workdir>/../reference/silvaXXX, "
                             "shared by all projects)")
    parser.add_argument("--threads", type=int, default=available_cpus(),
                        help=f"CPU threads for all tools (default: all visible = {available_cpus()})")
    parser.add_argument("--until", choices=STEPS, default="finalize", help="Stop after this step")
    parser.add_argument("--figaro", action="store_true", help="Run FIGARO (sets DADA2 truncation when dada2.truncation is 'figaro', "
                             "otherwise an advisory report)")
    parser.add_argument("--export-dir", type=Path, default=None,
                        help="Copy final tables/reports here, e.g. /mnt/e/BiomeFM/outputs")
    parser.add_argument("--preflight", action="store_true", help="Check tools/resources, print plan, run nothing")
    args = parser.parse_args()

    cfg = load_config(args.config)
    args.workdir = (args.workdir or Path.home() / "biomefm_runs" / dataset_name(cfg)).expanduser()
    args.workdir.mkdir(parents=True, exist_ok=True)

    log_dir = args.workdir / "cohorts" / args.mode / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)-7s %(message)s", datefmt="%H:%M:%S",
        handlers=[logging.StreamHandler(sys.stdout),
                  logging.FileHandler(log_dir / f"pipeline_{datetime.now():%Y%m%d_%H%M%S}.log")])

    if args.preflight:
        return 0 if preflight(args, cfg) else 1
    if not preflight(args, cfg):
        logger.error("Fix the preflight problems above first (see docs/16S_PREPROCESSING.md)")
        return 1
    pipeline = AmpliconPipeline(cfg, args.mode, args.workdir, args.threads,
                                reference_dir=args.reference_dir, run_figaro=args.figaro,
                                export_dir=args.export_dir)
    pipeline.run(args.until)
    return 0


if __name__ == "__main__":
    sys.exit(main())
