# Geng-BEVA fixed runtime

Geng-BEVA treats its scientific Python stack as a required execution environment, not as optional visualization extras.

## Windows Codex entry point

Use `.codex/run-geng-beva.ps1` as the canonical entry point. It always verifies the workspace runtime before analysis.

Examples:

```powershell
powershell -ExecutionPolicy Bypass -File .codex\run-geng-beva.ps1 population-structure-analysis pca data.vcf --out-dir outputs\pca
powershell -ExecutionPolicy Bypass -File .codex\run-geng-beva.ps1 population-structure-analysis ibs data.vcf --out-dir outputs\ibs
powershell -ExecutionPolicy Bypass -File .codex\run-geng-beva.ps1 igci-evaluation --vcf data.vcf --sample YF47 --out-dir outputs\igci
```

## Bootstrap behavior

The first run:

1. reuses `uv` from PATH when available;
2. otherwise downloads a fixed Windows `uv` bootstrap into `.codex/runtime/geng-beva/bootstrap/`;
3. uses `uv` to obtain managed Python 3.11 when necessary;
4. creates `.codex/runtime/geng-beva/.venv/`;
5. synchronizes the exact packages in `runtime-requirements.txt`;
6. verifies imports and exact versions;
7. records runtime provenance in `runtime_manifest.json`.

Later runs verify and reuse the same `.venv`. If it is incomplete or damaged, bootstrap repairs it before analysis.

## Required packages

The fixed runtime includes NumPy, pandas, Matplotlib, Plotly, SciPy, and scikit-learn. Plotly is mandatory because interactive PCA and IBS/IBD HTML outputs are part of the standard population-structure contract.

## PLINK exception

PLINK is deliberately not installed or upgraded by runtime bootstrap. Use the fixed executable under `skills/population-structure-analysis/tools/plink/windows-x64/plink.exe`, or explicitly configure `--plink-cmd` / `PLINK_BIN`.

## Failure semantics

A task that promises an interactive PCA or IBS/IBD visualization is not complete unless the required HTML file is created. Missing Plotly must trigger runtime repair and rerun; it must not be silently downgraded to a CSV-only result. Only report runtime failure after bootstrap/repair itself is blocked by network, sandbox, or filesystem policy.
