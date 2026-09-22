[CmdletBinding()]
param(
    [switch]$ForceRepair,
    [switch]$Quiet
)

$ErrorActionPreference = "Stop"
$CodexRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$RuntimeRoot = Join-Path $CodexRoot "runtime\geng-beva"
$VenvRoot = Join-Path $RuntimeRoot ".venv"
$PythonExe = Join-Path $VenvRoot "Scripts\python.exe"
$Requirements = Join-Path $CodexRoot "runtime-requirements.txt"
$Verifier = Join-Path $CodexRoot "tools\verify_runtime.py"
$Manifest = Join-Path $RuntimeRoot "runtime_manifest.json"
$BootstrapDir = Join-Path $RuntimeRoot "bootstrap"
$UvDir = Join-Path $BootstrapDir "uv"
$UvExe = Join-Path $UvDir "uv.exe"
$PythonVersion = (Get-Content (Join-Path $CodexRoot ".python-version") -Raw).Trim()
$UvVersion = "0.10.0"

New-Item -ItemType Directory -Force -Path $RuntimeRoot | Out-Null
New-Item -ItemType Directory -Force -Path $UvDir | Out-Null

function Write-Status([string]$Message) {
    if (-not $Quiet) { Write-Host "[Geng-BEVA runtime] $Message" }
}

function Get-Uv {
    $pathUv = Get-Command uv.exe -ErrorAction SilentlyContinue
    if (-not $pathUv) { $pathUv = Get-Command uv -ErrorAction SilentlyContinue }
    if ($pathUv) {
        Write-Status "Using uv from PATH: $($pathUv.Source)"
        return $pathUv.Source
    }
    if (Test-Path $UvExe) {
        Write-Status "Reusing workspace uv: $UvExe"
        return $UvExe
    }

    $arch = $env:PROCESSOR_ARCHITECTURE
    if ($arch -eq "ARM64") {
        $asset = "uv-aarch64-pc-windows-msvc.zip"
    } else {
        $asset = "uv-x86_64-pc-windows-msvc.zip"
    }
    $url = "https://github.com/astral-sh/uv/releases/download/$UvVersion/$asset"
    $zip = Join-Path $BootstrapDir $asset
    Write-Status "uv is missing; downloading fixed bootstrap uv $UvVersion to the workspace runtime."
    try {
        Invoke-WebRequest -Uri $url -OutFile $zip -UseBasicParsing
        $extract = Join-Path $BootstrapDir "uv-extract"
        if (Test-Path $extract) { Remove-Item -Recurse -Force $extract }
        Expand-Archive -Path $zip -DestinationPath $extract -Force
        $found = Get-ChildItem -Path $extract -Recurse -Filter uv.exe | Select-Object -First 1
        if (-not $found) { throw "uv.exe was not present in downloaded archive." }
        Copy-Item $found.FullName $UvExe -Force
        Remove-Item $zip -Force -ErrorAction SilentlyContinue
        Remove-Item $extract -Recurse -Force -ErrorAction SilentlyContinue
    } catch {
        throw "Could not bootstrap uv. Network or sandbox policy may block downloads/writes. $($_.Exception.Message)"
    }
    return $UvExe
}

function Test-Runtime {
    if (-not (Test-Path $PythonExe)) { return $false }
    & $PythonExe $Verifier *> $null
    return ($LASTEXITCODE -eq 0)
}

if ((-not $ForceRepair) -and (Test-Runtime)) {
    Write-Status "Verified fixed runtime; reusing it."
    Write-Output $PythonExe
    exit 0
}

$uv = Get-Uv
Write-Status "Preparing managed Python $PythonVersion."
& $uv python install $PythonVersion
if ($LASTEXITCODE -ne 0) { throw "uv could not install/locate Python $PythonVersion." }

if (-not (Test-Path $PythonExe)) {
    Write-Status "Creating isolated runtime at $VenvRoot."
    & $uv venv --python $PythonVersion $VenvRoot
    if ($LASTEXITCODE -ne 0) { throw "uv could not create the Geng-BEVA virtual environment." }
}

Write-Status "Synchronizing fixed Geng-BEVA dependencies."
& $uv pip sync --python $PythonExe $Requirements
if ($LASTEXITCODE -ne 0) { throw "uv could not synchronize Geng-BEVA runtime dependencies." }

$verifyJson = & $PythonExe $Verifier
if ($LASTEXITCODE -ne 0) {
    throw "Geng-BEVA runtime verification failed after repair.`n$verifyJson"
}

$reqHash = (Get-FileHash $Requirements -Algorithm SHA256).Hash
$manifestObject = [ordered]@{
    status = "verified"
    python_version = $PythonVersion
    python_executable = $PythonExe
    requirements_sha256 = $reqHash
    uv_version = $UvVersion
    updated_at = (Get-Date).ToString("o")
    verification = ($verifyJson | ConvertFrom-Json)
}
$manifestObject | ConvertTo-Json -Depth 8 | Set-Content -Encoding UTF8 $Manifest
Write-Status "Runtime is ready and will be reused on future runs."
Write-Output $PythonExe
