# Windows helper: create an Ubuntu 24.04 WSL distro named "BiomeFM" and install the
# 16S pipeline inside it. Run from the repository folder in PowerShell:
#   powershell -ExecutionPolicy Bypass -File envs\windows\setup_wsl.ps1
#
# Uses WSL2 when CPU virtualization is enabled in the BIOS/UEFI, otherwise WSL1
# (no virtualization needed; CPU-bound steps run at near-native speed).
# If the WSL Windows feature is missing, it is enabled (admin prompt) and you are
# asked to reboot once and run this script again.

param(
    [string]$Distro  = "BiomeFM",
    [string]$InstallDir = "$env:LOCALAPPDATA\BiomeFM\wsl"
)
$ErrorActionPreference = "Stop"
$env:WSL_UTF8 = "1"
$repo = (Resolve-Path "$PSScriptRoot\..\..").Path

$virt = (Get-CimInstance Win32_Processor | Select-Object -First 1).VirtualizationFirmwareEnabled
$ver  = if ($virt) { 2 } else { 1 }
Write-Host "CPU virtualization enabled in firmware: $virt -> WSL$ver"

$status = (wsl --status 2>&1 | Out-String)
if ($status -match "not supported|enable the .*optional component" -or -not (Get-Command wsl.exe -ErrorAction SilentlyContinue)) {
    Write-Host "Enabling the Windows Subsystem for Linux features (admin prompt)..."
    $cmd = "dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart; " +
           "dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart"
    Start-Process powershell.exe -Verb RunAs -Wait -ArgumentList "-NoProfile", "-Command", $cmd
    Write-Host "`nReboot Windows, then run this script again." -ForegroundColor Yellow
    exit 0
}

$distros = (wsl --list --quiet) -split "`r?`n" | ForEach-Object { $_.Trim() } | Where-Object { $_ }
if ($distros -notcontains $Distro) {
    $tar = "$InstallDir\ubuntu-24.04-rootfs.tar.gz"
    New-Item -ItemType Directory -Force $InstallDir | Out-Null
    if (-not (Test-Path $tar)) {
        Write-Host "Downloading Ubuntu 24.04 WSL image (~340 MB)..."
        $ProgressPreference = "SilentlyContinue"
        Invoke-WebRequest "https://cloud-images.ubuntu.com/wsl/releases/24.04/current/ubuntu-noble-wsl-amd64-24.04lts.rootfs.tar.gz" `
            -OutFile $tar -UseBasicParsing
    }
    Write-Host "Importing as WSL distro '$Distro' (WSL$ver)..."
    wsl --import $Distro "$InstallDir\$Distro" $tar --version $ver
    if ($LASTEXITCODE) { throw "wsl --import failed ($LASTEXITCODE)" }
}

# Copy the repo into the Linux filesystem (much faster than working on /mnt/c).
$wslRepo = (wsl -d $Distro -u root -- wslpath -a "$repo").Trim()
wsl -d $Distro -u root -- bash -c "rm -rf /root/microbiome_AI && cp -r '$wslRepo' /root/microbiome_AI && find /root/microbiome_AI -name '*.sh' -exec sed -i 's/\r$//' {} + && cd /root/microbiome_AI && bash envs/setup_local.sh"
if ($LASTEXITCODE) { throw "setup_local.sh failed ($LASTEXITCODE)" }

Write-Host "`nSetup complete. To run the pipeline:" -ForegroundColor Green
Write-Host "  wsl -d $Distro -u root"
Write-Host "  cd ~/microbiome_AI && bash Scripts/amplicon_16s/run_end_to_end.sh --dry-run-only"
