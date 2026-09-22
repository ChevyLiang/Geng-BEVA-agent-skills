# Geng-BEVA managed runtime

This directory is a persistent workspace-local runtime cache.

On Windows, `.codex/run-geng-beva.ps1` invokes `tools/bootstrap_runtime.ps1` before every skill run. The bootstrap verifies the existing runtime and reuses it when healthy. If missing or damaged, it repairs the runtime using a workspace-local `uv` bootstrap and managed Python 3.11.

The generated `.venv/`, `bootstrap/`, and `runtime_manifest.json` are runtime artifacts and should not be committed to GitHub.

PLINK is intentionally excluded from this runtime. Population-structure analysis continues to use the fixed PLINK executable supplied separately under `skills/population-structure-analysis/tools/plink/windows-x64/plink.exe` (or an explicitly configured executable).
