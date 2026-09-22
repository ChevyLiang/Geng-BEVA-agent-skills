[CmdletBinding()]
param(
    [Parameter(Mandatory=$true, Position=0)]
    [ValidateSet("population-structure-analysis","functional-gene-diagnostic","igci-evaluation","transgene-detection")]
    [string]$Skill,

    [Parameter(ValueFromRemainingArguments=$true, Position=1)]
    [string[]]$SkillArgs
)

$ErrorActionPreference = "Stop"
$CodexRoot = $PSScriptRoot
$Bootstrap = Join-Path $CodexRoot "tools\bootstrap_runtime.ps1"
& $Bootstrap -Quiet | Out-Null
if ($LASTEXITCODE -ne 0) { throw "Geng-BEVA runtime bootstrap failed." }

$PythonExe = Join-Path $CodexRoot "runtime\geng-beva\.venv\Scripts\python.exe"
if (-not (Test-Path $PythonExe)) { throw "Verified runtime Python was not found at $PythonExe" }

switch ($Skill) {
    "igci-evaluation" {
        $Entry = Join-Path $CodexRoot "skills\igci-evaluation\scripts\igci_evaluation.py"
    }
    default {
        $Entry = Join-Path $CodexRoot "skills\$Skill\main.py"
    }
}

$env:GENGBEVA_FIXED_RUNTIME = "1"
& $PythonExe $Entry @SkillArgs
exit $LASTEXITCODE
