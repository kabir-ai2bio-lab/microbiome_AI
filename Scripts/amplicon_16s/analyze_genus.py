#!/usr/bin/env python3
"""
Automatic downstream analysis of a 16S genus matrix produced by
preprocess_16s_amplicon.py. Writes one self-contained HTML report plus CSVs:

    1. QC          - reads per sample, DADA2 retention by group
    2. Composition - top genera per group
    3. Alpha       - Shannon / Simpson / observed genera (rarefied when counts exist),
                     Kruskal-Wallis + pairwise Mann-Whitney (BH-FDR)
    4. Beta        - Bray-Curtis PCoA, PERMANOVA overall and pairwise
    5. Differential abundance - CLR + Kruskal-Wallis / Mann-Whitney (BH-FDR)
    6. ML baseline - Random forest and L2 logistic regression, repeated stratified
                     CV ROC-AUC, label-permutation test, RF feature importance
    7. Reproduction check (optional) - compare against a reference matrix (e.g. Colab)

Steps that need more samples than available (e.g. a 3-sample dry run) are skipped
with a note instead of failing. Model fitting uses every CPU (n_jobs=-1).

Examples:
    python analyze_genus.py --cohort-dir ~/biomefm_runs/PRJNA496408/cohorts/full
    python analyze_genus.py --genus-csv genus_relative_with_provisional_labels.csv --out-dir report/
    python analyze_genus.py --cohort-dir .../cohorts/full --compare-to colab/genus_relative_with_provisional_labels.csv
"""

import argparse
import base64
import io
import json
import logging
import sys
import time
from datetime import datetime
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import (RepeatedStratifiedKFold, StratifiedKFold, cross_val_predict,
                                     cross_validate, permutation_test_score)
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

logger = logging.getLogger("analyze_genus")

LABEL = "Diagnosis_provisional"
GROUP_ORDER = ["NC", "MCI", "AD"]
COLORS = {"NC": "#4C78A8", "MCI": "#E8A33D", "AD": "#C8553D"}
NOT_A_GENUS = {"Unassigned", "uncultured"}


def bh_fdr(p) -> np.ndarray:
    """Benjamini-Hochberg q-values; NaN p-values stay NaN."""
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    n = ok.sum()
    if n:
        order = np.argsort(p[ok])
        ranked = p[ok][order] * n / np.arange(1, n + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        vals = np.empty(n)
        vals[order] = np.minimum(ranked, 1.0)
        q[ok] = vals
    return q


def fig_to_base64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


class Report:
    """Collects sections (text, tables, figures) and renders one HTML file."""

    def __init__(self, title: str):
        self.title = title
        self.parts: List[str] = []
        self.summary: Dict = {}

    def h(self, text: str):
        self.parts.append(f"<h2>{text}</h2>")

    def p(self, text: str):
        self.parts.append(f"<p>{text}</p>")

    def note(self, text: str):
        self.parts.append(f'<p class="note">{text}</p>')

    def table(self, df: pd.DataFrame, float_format="{:.4g}".format, index=True):
        self.parts.append(df.to_html(float_format=float_format, index=index, border=0, classes="t"))

    def figure(self, fig, caption: str = ""):
        self.parts.append(f'<figure><img src="data:image/png;base64,{fig_to_base64(fig)}">'
                          f"<figcaption>{caption}</figcaption></figure>")

    def write(self, path: Path):
        css = """
        :root{--bg:#fff;--fg:#1f2328;--muted:#59636e;--line:#d1d9e0;--note:#fff8e6}
        @media (prefers-color-scheme:dark){:root{--bg:#0d1117;--fg:#e6edf3;--muted:#9198a1;--line:#30363d;--note:#2b2410}}
        body{background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,sans-serif;max-width:1100px;margin:auto;padding:16px}
        h1{font-size:1.6em}h2{border-bottom:1px solid var(--line);padding-bottom:4px;margin-top:2em}
        table.t{border-collapse:collapse;font-size:13px;margin:8px 0;display:block;overflow-x:auto}
        table.t th,table.t td{border-bottom:1px solid var(--line);padding:4px 8px;text-align:right}
        img{max-width:100%;background:#fff;border-radius:4px}figcaption{color:var(--muted);font-size:13px}
        .note{background:var(--note);padding:8px 12px;border-radius:6px}"""
        html = (f"<!doctype html><html><head><meta charset='utf-8'>"
                f"<meta name='viewport' content='width=device-width,initial-scale=1'>"
                f"<title>{self.title}</title><style>{css}</style></head><body>"
                f"<h1>{self.title}</h1><p>Generated {datetime.now():%Y-%m-%d %H:%M}</p>"
                + "\n".join(self.parts) + "</body></html>")
        path.write_text(html, encoding="utf-8")


class GenusAnalysis:

    def __init__(self, rel: pd.DataFrame, labels: pd.Series, counts: Optional[pd.DataFrame],
                 dada2_stats: Optional[pd.DataFrame], out_dir: Path, args):
        self.rel, self.labels, self.counts, self.dada2_stats = rel, labels, counts, dada2_stats
        self.out, self.args = out_dir, args
        self.rng = np.random.default_rng(args.seed)
        self.groups = [g for g in GROUP_ORDER if g in set(labels)] + \
                      sorted(set(labels) - set(GROUP_ORDER))
        self.genera = [c for c in rel.columns if c not in NOT_A_GENUS]
        self.report = Report(args.title)
        self.min_group = labels.value_counts().min()

    def color(self, g):
        return COLORS.get(g, "#777777")

    # ----- 1. QC ------------------------------------------------------------

    def section_qc(self):
        r = self.report
        r.h("1. Cohort and read QC")
        n = self.labels.value_counts().reindex(self.groups)
        r.p(f"{len(self.rel)} samples, {self.rel.shape[1]} genus-level categories "
            f"({len(self.genera)} named genera). Group sizes: "
            + ", ".join(f"{g} = {n[g]}" for g in self.groups) + ".")
        r.note("Labels are <b>provisional</b>: inferred from NCBI LibraryName prefixes and "
               "not clinically verified.")
        r.summary["n_samples"] = int(len(self.rel))
        r.summary["group_sizes"] = {g: int(n[g]) for g in self.groups}
        unassigned = self.rel["Unassigned"] if "Unassigned" in self.rel else pd.Series(0.0, index=self.rel.index)
        r.summary["median_unassigned_fraction"] = float(unassigned.median())
        r.p(f"Median fraction of reads without a genus assignment: {unassigned.median():.2%}.")

        if self.dada2_stats is not None:
            s = self.dada2_stats.set_index("sample-id").reindex(self.rel.index)
            s["group"] = self.labels
            tab = s.groupby("group")[["input", "non-chimeric", "retained_pct"]].median().reindex(self.groups)
            r.table(tab.rename(columns=lambda c: f"median {c}"))
            fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
            for i, (col, title) in enumerate((("non-chimeric", "Final reads per sample"),
                                              ("retained_pct", "DADA2 retention (%)"))):
                data = [s.loc[s["group"] == g, col].dropna() for g in self.groups]
                bp = ax[i].boxplot(data, patch_artist=True)
                ax[i].set_xticks(range(1, len(self.groups) + 1), self.groups)
                for patch, g in zip(bp["boxes"], self.groups):
                    patch.set_facecolor(self.color(g))
                ax[i].set_title(title)
            r.figure(fig, "Sequencing depth and DADA2 retention by group. Depth differences "
                          "between groups can confound diversity comparisons.")
            if s[["non-chimeric"]].notna().all().all() and len(self.groups) > 1 and self.min_group >= 2:
                p = stats.kruskal(*[s.loc[s["group"] == g, "non-chimeric"] for g in self.groups]).pvalue
                r.p(f"Kruskal-Wallis on final read depth across groups: p = {p:.3g}.")
                r.summary["depth_kruskal_p"] = float(p)

    # ----- 2. composition ---------------------------------------------------

    def section_composition(self):
        r = self.report
        r.h("2. Composition")
        k = self.args.top_n
        means = self.rel.groupby(self.labels).mean().reindex(self.groups)
        top = self.rel[self.genera].mean().nlargest(k).index.tolist()
        plot = means[top].copy()
        plot["Other"] = 1 - plot.sum(axis=1)
        fig, ax = plt.subplots(figsize=(9, 4.5))
        bottom = np.zeros(len(plot))
        cmap = plt.get_cmap("tab20")
        for i, col in enumerate(plot.columns):
            ax.bar(plot.index, plot[col], bottom=bottom, label=col,
                   color="#bbbbbb" if col == "Other" else cmap(i % 20))
            bottom += plot[col].to_numpy()
        ax.set_ylabel("Mean relative abundance")
        ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
        r.figure(fig, f"Mean relative abundance of the {k} most abundant genera per group.")
        tab = means[top].T
        tab["overall"] = self.rel[top].mean()
        tab["prevalence"] = (self.rel[top] > 0).mean()
        r.table(tab)
        tab.to_csv(self.out / "top_genera_by_group.csv")

    # ----- 3. alpha diversity -----------------------------------------------

    def alpha_table(self) -> pd.DataFrame:
        if self.counts is not None:
            base = self.counts.loc[self.rel.index, [c for c in self.genera if c in self.counts]]
            depth = int(base.sum(axis=1).min())
            rows = [self.rng.multivariate_hypergeometric(row.astype(np.int64), depth)
                    for row in base.to_numpy()]
            base = pd.DataFrame(rows, index=base.index, columns=base.columns)
            self.rarefied_depth = depth
        else:
            base = self.rel[self.genera]
            self.rarefied_depth = None
        p = base.div(base.sum(axis=1), axis=0)
        logp = np.log(p.where(p > 0))
        return pd.DataFrame({
            "shannon": -(p * logp).sum(axis=1),
            "simpson": 1 - (p ** 2).sum(axis=1),
            "observed_genera": (base > 0).sum(axis=1),
        })

    def group_tests(self, values: pd.Series) -> dict:
        data = {g: values[self.labels == g].dropna() for g in self.groups}
        out = {"kruskal_p": np.nan}
        try:
            out["kruskal_p"] = stats.kruskal(*data.values()).pvalue
        except ValueError:
            pass
        pairs = list(combinations(self.groups, 2))
        ps = []
        for a, b in pairs:
            try:
                ps.append(stats.mannwhitneyu(data[a], data[b], alternative="two-sided").pvalue)
            except ValueError:
                ps.append(np.nan)
        for (a, b), pv, qv in zip(pairs, ps, bh_fdr(ps)):
            out[f"{a}_vs_{b}_p"], out[f"{a}_vs_{b}_q"] = pv, qv
        return out

    def section_alpha(self):
        r = self.report
        r.h("3. Alpha diversity")
        alpha = self.alpha_table()
        alpha.join(self.labels).to_csv(self.out / "alpha_diversity.csv")
        if self.rarefied_depth:
            r.p(f"Computed on genus counts rarefied to {self.rarefied_depth:,} reads per sample "
                f"(seed {self.args.seed}); 'Unassigned'/'uncultured' excluded.")
        else:
            r.p("Computed on relative abundances (no counts available), so richness depends on depth.")
        if self.min_group < 2 or len(self.groups) < 2:
            r.note("Too few samples per group for statistical tests; values are reported only.")
            r.table(alpha.join(self.labels))
            return
        tests = pd.DataFrame({m: self.group_tests(alpha[m]) for m in alpha.columns}).T
        r.table(tests)
        tests.to_csv(self.out / "alpha_tests.csv")
        r.summary["alpha_kruskal_p"] = tests["kruskal_p"].to_dict()
        fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
        for ax, m in zip(axes, alpha.columns):
            data = [alpha.loc[self.labels == g, m] for g in self.groups]
            bp = ax.boxplot(data, patch_artist=True, showfliers=False)
            for patch, g in zip(bp["boxes"], self.groups):
                patch.set_facecolor(self.color(g))
                patch.set_alpha(0.6)
            for i, d in enumerate(data, 1):
                ax.scatter(np.full(len(d), i) + self.rng.uniform(-0.12, 0.12, len(d)), d, s=10, c="k", alpha=0.6)
            ax.set_xticks(range(1, len(self.groups) + 1), self.groups)
            ax.set_title(f"{m}  (KW p={tests.loc[m, 'kruskal_p']:.3g})")
        r.figure(fig, "Alpha diversity by group (Kruskal-Wallis p in titles).")

    # ----- 4. beta diversity ------------------------------------------------

    @staticmethod
    def pcoa(d: np.ndarray):
        n = d.shape[0]
        j = np.eye(n) - np.ones((n, n)) / n
        b = -0.5 * j @ (d ** 2) @ j
        w, v = np.linalg.eigh(b)
        order = np.argsort(w)[::-1]
        w, v = w[order], v[:, order]
        pos = w > 1e-10
        return v[:, pos] * np.sqrt(w[pos]), w[pos] / w[pos].sum()

    def permanova(self, d: np.ndarray, groups: np.ndarray, n_perm: int):
        n = len(groups)
        d2 = d ** 2
        iu = np.triu_indices(n, 1)
        sst = d2[iu].sum() / n

        def ssw_of(g):
            total = 0.0
            for lvl in np.unique(g):
                idx = np.flatnonzero(g == lvl)
                sub = d2[np.ix_(idx, idx)]
                total += sub[np.triu_indices(len(idx), 1)].sum() / len(idx)
            return total

        a = len(np.unique(groups))
        ssw = ssw_of(groups)
        f = ((sst - ssw) / (a - 1)) / (ssw / (n - a))
        perm_f = np.empty(n_perm)
        for i in range(n_perm):
            sw = ssw_of(self.rng.permutation(groups))
            perm_f[i] = ((sst - sw) / (a - 1)) / (sw / (n - a))
        return {"pseudo_F": f, "R2": (sst - ssw) / sst,
                "p": (1 + (perm_f >= f).sum()) / (n_perm + 1), "n": n}

    def section_beta(self):
        r = self.report
        r.h("4. Beta diversity (Bray-Curtis)")
        x = self.rel.to_numpy()
        d = squareform(pdist(x, metric="braycurtis"))
        coords, expl = self.pcoa(d)
        pd.DataFrame(coords[:, :5], index=self.rel.index,
                     columns=[f"PC{i + 1}" for i in range(min(5, coords.shape[1]))]).join(self.labels) \
            .to_csv(self.out / "pcoa_braycurtis.csv")
        fig, ax = plt.subplots(figsize=(6, 5))
        for g in self.groups:
            m = (self.labels == g).to_numpy()
            ax.scatter(coords[m, 0], coords[m, 1] if coords.shape[1] > 1 else np.zeros(m.sum()),
                       label=f"{g} (n={m.sum()})", c=self.color(g), s=28, alpha=0.8, edgecolors="none")
        ax.set_xlabel(f"PC1 ({expl[0]:.1%})")
        if len(expl) > 1:
            ax.set_ylabel(f"PC2 ({expl[1]:.1%})")
        ax.legend()
        if self.min_group < 2 or len(self.groups) < 2:
            r.figure(fig, "Bray-Curtis PCoA.")
            r.note("Too few samples per group for PERMANOVA.")
            return
        g = self.labels.to_numpy()
        rows = {"all groups": self.permanova(d, g, self.args.n_permutations)}
        for a, b in combinations(self.groups, 2):
            m = np.isin(g, [a, b])
            rows[f"{a} vs {b}"] = self.permanova(d[np.ix_(m, m)], g[m], self.args.n_permutations)
        res = pd.DataFrame(rows).T
        res["q (pairwise BH)"] = np.nan
        pair_idx = res.index[1:]
        res.loc[pair_idx, "q (pairwise BH)"] = bh_fdr(res.loc[pair_idx, "p"].astype(float))
        r.figure(fig, f"Bray-Curtis PCoA. PERMANOVA (all groups): R² = {res.iloc[0]['R2']:.3f}, "
                      f"p = {res.iloc[0]['p']:.3g}.")
        r.table(res)
        res.to_csv(self.out / "permanova_braycurtis.csv")
        r.summary["permanova"] = {k: {"R2": float(v["R2"]), "p": float(v["p"])} for k, v in rows.items()}

    # ----- 5. differential abundance ----------------------------------------

    def clr_features(self) -> pd.DataFrame:
        keep = [c for c in self.genera if (self.rel[c] > 0).mean() >= self.args.min_prevalence]
        x = self.rel[keep]
        pseudo = x[x > 0].min().min() / 2
        logx = np.log(x + pseudo)
        return logx.sub(logx.mean(axis=1), axis=0)

    def section_differential(self):
        r = self.report
        r.h("5. Differential abundance (CLR, nonparametric)")
        clr = self.clr_features()
        self.clr = clr
        r.p(f"{clr.shape[1]} genera present in ≥{self.args.min_prevalence:.0%} of samples; CLR with "
            "half-minimum pseudocount. Kruskal-Wallis across groups, pairwise Mann-Whitney; BH-FDR "
            "within each test.")
        if self.min_group < 3 or len(self.groups) < 2:
            r.note("Too few samples per group for differential abundance testing; skipped.")
            return
        rows = {}
        for gname in clr.columns:
            t = self.group_tests(clr[gname])
            for g in self.groups:
                t[f"mean_rel_{g}"] = self.rel.loc[self.labels == g, gname].mean()
            for a, b in combinations(self.groups, 2):
                t[f"{a}_vs_{b}_median_clr_diff"] = (clr.loc[self.labels == b, gname].median()
                                                    - clr.loc[self.labels == a, gname].median())
            rows[gname] = t
        da = pd.DataFrame(rows).T
        da["kruskal_q"] = bh_fdr(da["kruskal_p"])
        for a, b in combinations(self.groups, 2):
            da[f"{a}_vs_{b}_q"] = bh_fdr(da[f"{a}_vs_{b}_p"])
        da = da.sort_values("kruskal_p")
        da.to_csv(self.out / "differential_abundance.csv")
        n_sig = int((da["kruskal_q"] < 0.05).sum())
        r.summary["da_kruskal_q_lt_0.05"] = n_sig
        r.summary["da_top10"] = da.index[:10].tolist()
        r.p(f"Genera with Kruskal-Wallis q &lt; 0.05: <b>{n_sig}</b> of {len(da)}.")
        cols = ["kruskal_p", "kruskal_q"] + [f"mean_rel_{g}" for g in self.groups] + \
               [c for c in da.columns if c.endswith("_q") and "_vs_" in c]
        r.table(da[cols].head(20))

        top = da.index[:8]
        fig, axes = plt.subplots(2, 4, figsize=(13, 6))
        for ax, gname in zip(axes.flat, top):
            data = [clr.loc[self.labels == g, gname] for g in self.groups]
            bp = ax.boxplot(data, patch_artist=True, showfliers=False)
            for patch, g in zip(bp["boxes"], self.groups):
                patch.set_facecolor(self.color(g))
                patch.set_alpha(0.6)
            for i, d in enumerate(data, 1):
                ax.scatter(np.full(len(d), i) + self.rng.uniform(-0.12, 0.12, len(d)), d, s=6, c="k", alpha=0.5)
            ax.set_xticks(range(1, len(self.groups) + 1), self.groups)
            ax.set_title(f"{gname[:28]}\nq={da.loc[gname, 'kruskal_q']:.2g}", fontsize=9)
        for ax in list(axes.flat)[len(top):]:
            ax.axis("off")
        fig.tight_layout()
        r.figure(fig, "CLR abundance of the 8 genera with the smallest Kruskal-Wallis p-values.")

        if {"NC", "AD"} <= set(self.groups):
            x = da["NC_vs_AD_median_clr_diff"].astype(float)
            y = -np.log10(da["NC_vs_AD_p"].astype(float))
            sig = da["NC_vs_AD_q"].astype(float) < 0.05
            fig, ax = plt.subplots(figsize=(6.5, 4.5))
            ax.scatter(x[~sig], y[~sig], s=12, c="#999999")
            ax.scatter(x[sig], y[sig], s=16, c=COLORS["AD"])
            for gname in y.nlargest(8).index:
                ax.annotate(gname[:22], (x[gname], y[gname]), fontsize=7)
            ax.set_xlabel("Median CLR difference (AD − NC)")
            ax.set_ylabel("−log10 p (Mann-Whitney)")
            r.figure(fig, "AD vs NC volcano plot; red = q &lt; 0.05.")

    # ----- 6. machine learning baseline -------------------------------------

    def section_ml(self):
        r = self.report
        r.h("6. Machine-learning baseline")
        a = self.args
        tasks = {}
        for x, y in (("AD", "NC"), ("MCI", "NC"), ("AD", "MCI")):
            if {x, y} <= set(self.groups):
                tasks[f"{x}_vs_{y}"] = [x, y]
        if len(self.groups) > 2:
            tasks["_vs_".join(reversed(self.groups))] = self.groups
        if self.min_group < a.cv_folds or not tasks:
            r.note(f"Needs ≥{a.cv_folds} samples per group for {a.cv_folds}-fold CV; skipped.")
            return
        X_all = self.clr
        models = {
            "RandomForest": RandomForestClassifier(n_estimators=500, class_weight="balanced",
                                                   min_samples_leaf=2, random_state=a.seed, n_jobs=1),
            "LogisticL2": make_pipeline(StandardScaler(), LogisticRegression(
                C=0.1, max_iter=5000, class_weight="balanced")),
        }
        r.p(f"Features: CLR of the {X_all.shape[1]} prevalence-filtered genera (filter is label-free). "
            f"{a.cv_folds}-fold stratified CV × {a.cv_repeats} repeats; {a.n_ml_permutations} label "
            "permutations give the chance-level p-value. Results are exploratory baselines, "
            "not a validated classifier.")
        rows, roc_data = [], {}
        for task, groups in tasks.items():
            m = self.labels.isin(groups).to_numpy()
            X, y = X_all.to_numpy()[m], self.labels.to_numpy()[m]
            binary = len(groups) == 2
            scoring = "roc_auc" if binary else "roc_auc_ovr"
            # positive class = disease (first in the task name)
            y_enc = (y == groups[0]).astype(int) if binary else y
            cv = RepeatedStratifiedKFold(n_splits=a.cv_folds, n_repeats=a.cv_repeats, random_state=a.seed)
            for name, model in models.items():
                t0 = time.time()
                res = cross_validate(model, X, y_enc, cv=cv, n_jobs=-1,
                                     scoring={"auc": scoring, "bal_acc": "balanced_accuracy"})
                _, _, perm_p = permutation_test_score(
                    model, X, y_enc, scoring=scoring, n_permutations=a.n_ml_permutations, n_jobs=-1,
                    cv=StratifiedKFold(a.cv_folds, shuffle=True, random_state=a.seed), random_state=a.seed)
                rows.append({"task": task, "model": name, "n": int(m.sum()),
                             "AUC_mean": res["test_auc"].mean(), "AUC_sd": res["test_auc"].std(),
                             "bal_acc_mean": res["test_bal_acc"].mean(), "perm_p": perm_p,
                             "seconds": round(time.time() - t0, 1)})
                logger.info("ML %-22s %-12s AUC %.3f ± %.3f (perm p=%.3g)", task, name,
                            res["test_auc"].mean(), res["test_auc"].std(), perm_p)
                if binary:
                    proba = cross_val_predict(model, X, y_enc, method="predict_proba", n_jobs=-1,
                                              cv=StratifiedKFold(a.cv_folds, shuffle=True, random_state=a.seed))[:, 1]
                    roc_data[(task, name)] = roc_curve(y_enc, proba)
        ml = pd.DataFrame(rows)
        ml.to_csv(self.out / "ml_baseline.csv", index=False)
        r.table(ml, index=False)
        r.summary["ml"] = ml[["task", "model", "AUC_mean", "AUC_sd", "perm_p"]].to_dict("records")

        binary_tasks = [t for t in tasks if len(tasks[t]) == 2]
        if binary_tasks:
            fig, axes = plt.subplots(1, len(binary_tasks), figsize=(4.2 * len(binary_tasks), 4))
            axes = np.atleast_1d(axes)
            for ax, task in zip(axes, binary_tasks):
                for name, ls in zip(models, ("-", "--")):
                    fpr, tpr, _ = roc_data[(task, name)]
                    ax.plot(fpr, tpr, ls, label=f"{name} (AUC {auc(fpr, tpr):.2f})")
                ax.plot([0, 1], [0, 1], ":", c="#999999")
                ax.set_title(task)
                ax.set_xlabel("FPR")
                ax.set_ylabel("TPR")
                ax.legend(fontsize=8, loc="lower right")
            fig.tight_layout()
            r.figure(fig, "ROC curves from one 5-fold cross-validated prediction per sample.")

        if "AD_vs_NC" in tasks:
            mask = self.labels.isin(["AD", "NC"]).to_numpy()
            rf = clone(models["RandomForest"]).set_params(n_jobs=-1)
            rf.fit(X_all.to_numpy()[mask], (self.labels.to_numpy()[mask] == "AD").astype(int))
            imp = pd.Series(rf.feature_importances_, index=X_all.columns).nlargest(20)
            imp.to_csv(self.out / "rf_importance_AD_vs_NC.csv", header=["importance"])
            fig, ax = plt.subplots(figsize=(7, 5))
            ax.barh(imp.index[::-1], imp.values[::-1], color=COLORS["AD"])
            ax.set_xlabel("Impurity importance")
            r.figure(fig, "Top 20 genera by random-forest importance (AD vs NC, fit on all samples; "
                          "descriptive, not a selection procedure).")

    # ----- 7. reproduction check --------------------------------------------

    def section_compare(self, ref_path: Path):
        r = self.report
        r.h("7. Reproduction check against reference matrix")
        ref = pd.read_csv(ref_path, index_col=0)
        ref_labels = ref.pop(LABEL) if LABEL in ref else None
        ref = ref.astype(float)
        common = self.rel.index.intersection(ref.index)
        cols = self.rel.columns.union(ref.columns)
        a = self.rel.reindex(index=common, columns=cols, fill_value=0.0)
        b = ref.reindex(index=common, columns=cols, fill_value=0.0)
        per_sample = pd.Series({s: np.corrcoef(a.loc[s], b.loc[s])[0, 1] for s in common})
        diff = (a - b).abs()
        out = {
            "reference": str(ref_path), "samples_here": len(self.rel), "samples_reference": len(ref),
            "samples_common": len(common),
            "genera_only_here": sorted(set(self.rel.columns) - set(ref.columns)),
            "genera_only_reference": sorted(set(ref.columns) - set(self.rel.columns)),
            "per_sample_pearson_min": float(per_sample.min()) if len(common) else None,
            "per_sample_pearson_median": float(per_sample.median()) if len(common) else None,
            "max_abs_diff": float(diff.to_numpy().max()) if len(common) else None,
            "identical": bool(len(common) and np.allclose(a, b, atol=1e-9)),
        }
        if ref_labels is not None and len(common):
            out["labels_agree"] = bool((ref_labels.reindex(common) == self.labels.reindex(common)).all())
        r.summary["reproduction"] = {k: v for k, v in out.items() if not isinstance(v, list)}
        verdict = ("identical to the reference" if out["identical"] else
                   f"per-sample Pearson r median {out['per_sample_pearson_median']:.4f} "
                   f"(min {out['per_sample_pearson_min']:.4f}), max |Δ| = {out['max_abs_diff']:.3g}")
        r.p(f"{len(common)} common samples; local matrix is {verdict}.")
        r.p(f"Genera only here: {len(out['genera_only_here'])}; only in reference: "
            f"{len(out['genera_only_reference'])}.")
        (self.out / "reproduction_check.json").write_text(json.dumps(out, indent=2))
        logger.info("Reproduction check: %s", verdict)

    def run(self, compare_to: Optional[Path]):
        for name in ("qc", "composition", "alpha", "beta", "differential", "ml"):
            t0 = time.time()
            logger.info("Section: %s", name)
            getattr(self, "section_" + name)()
            logger.info("  %.1fs", time.time() - t0)
        if compare_to:
            self.section_compare(compare_to)
        self.report.h("Caveats")
        self.report.note("Single BioProject (PRJNA496408); provisional labels; genus taxonomy at a "
                         "permissive 90% VSEARCH identity; no cross-study batch correction. Treat all "
                         "p-values and AUCs as exploratory.")
        html = self.out / "analysis_report.html"
        self.report.write(html)
        (self.out / "analysis_summary.json").write_text(json.dumps(self.report.summary, indent=2, default=float))
        logger.info("Report: %s", html)


def main() -> int:
    parser = argparse.ArgumentParser(description="Automatic analysis of a 16S genus matrix",
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=__doc__.split("Examples:")[1])
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--cohort-dir", type=Path, help="cohorts/<mode> folder from the preprocessing pipeline")
    src.add_argument("--genus-csv", type=Path, help="genus_relative_with_provisional_labels.csv")
    parser.add_argument("--out-dir", type=Path, help="default: <cohort-dir>/analysis")
    parser.add_argument("--compare-to", type=Path, help="reference relative-abundance CSV (e.g. Colab output)")
    parser.add_argument("--title", default="PRJNA496408 16S genus analysis")
    parser.add_argument("--top-n", type=int, default=15)
    parser.add_argument("--min-prevalence", type=float, default=0.10)
    parser.add_argument("--n-permutations", type=int, default=999, help="PERMANOVA permutations")
    parser.add_argument("--cv-folds", type=int, default=5)
    parser.add_argument("--cv-repeats", type=int, default=10)
    parser.add_argument("--n-ml-permutations", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20261011)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)-7s %(message)s",
                        datefmt="%H:%M:%S", stream=sys.stdout)

    counts = dada2_stats = None
    if args.cohort_dir:
        gdir = args.cohort_dir / "genus"
        rel_path = gdir / "genus_relative_with_provisional_labels.csv"
        if (gdir / "genus_counts.csv").exists():
            counts = pd.read_csv(gdir / "genus_counts.csv", index_col=0)
        if (args.cohort_dir / "dada2_per_sample_stats.csv").exists():
            dada2_stats = pd.read_csv(args.cohort_dir / "dada2_per_sample_stats.csv")
        out = args.out_dir or args.cohort_dir / "analysis"
    else:
        rel_path = args.genus_csv
        out = args.out_dir or rel_path.parent / "analysis"
    out.mkdir(parents=True, exist_ok=True)

    rel = pd.read_csv(rel_path, index_col=0)
    if LABEL not in rel:
        raise SystemExit(f"{rel_path} has no {LABEL} column")
    labels = rel.pop(LABEL)
    rel = rel.astype(float)
    if not np.allclose(rel.sum(axis=1), 1.0, atol=1e-6):
        logger.warning("Rows do not sum to 1 (filtered matrix?); continuing without renormalizing")

    GenusAnalysis(rel, labels, counts, dada2_stats, out, args).run(args.compare_to)
    return 0


if __name__ == "__main__":
    sys.exit(main())
