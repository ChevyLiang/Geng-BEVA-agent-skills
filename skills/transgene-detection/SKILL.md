---
name: transgene-detection
description: Evaluate sample-level evidence from transgene assay records encoded in standardized Geng-BEVA VCF files, including known BT, 35S, CP4_EPSPS and related element aliases, using sample GT/DP and configurable heuristic depth tiers. Use for Geng-BEVA transgene-assay evidence review, not for de novo insertion discovery or genome-wide transgene calling from sequencing reads.
---

# Transgene assay evidence

Use `main.py` to evaluate project-specific transgene assay records and, optionally, export a canonical chromosome-1-12 VCF. Do not describe this workflow as de novo transgene insertion discovery.

Under the Geng-BEVA VCF convention used by this repository, non-chromosome-1-12 **data rows** represent assay-element records. Header metadata is never evidence. Known aliases are grouped into families; unknown noncanonical record names retain their original identity. Do not generalize this convention to unrelated VCFs.

## Evidence rules

Read sample `FORMAT/GT` and `FORMAT/DP`; do not substitute site-level `INFO/DP` for missing sample depth. Require a fully called diploid GT; missing or partial-missing GT gives `no_data` regardless of DP.

The default evidence tiers are heuristic and configurable:

| Sample DP | Default `EvidenceLevel` |
| --- | --- |
| >100 | `high_support` |
| 50-100 | `moderate_support` |
| 10-<50 | `low_support` |
| 0-<10 | `insufficient_support` |
| Missing/invalid DP or missing GT | `no_data` |

Treat these cutoffs as operational defaults, not validated biological detection thresholds. DP below 10 is insufficient evidence, not proof of absence. `HasTransgeneSupport` means at least one fully called assay record meets the configured low-support threshold; it does not establish an experimentally confirmed insertion.

Do not infer insertion position, copy number, sensitivity, specificity, inheritance, or phenotype from this workflow.

Use the exact input supplied by the user and preserve source files. Example files are for demonstration only and must not replace the user's biological input.

Read:
- `references/evidence_rules.md` for interpretation.
- `references/usage.md` for commands and configurable thresholds.
- `references/contract.md` for the Python/manifest contract.

## Fixed runtime and self-healing

On Windows Codex, use the workspace entry point `.codex/run-geng-beva.ps1`; do not treat the active system/Conda Python as the authoritative environment. Before analysis, the entry point verifies `.codex/runtime/geng-beva/.venv`. If the runtime is missing or damaged, it bootstraps/repairs the managed Python 3.11 environment and synchronizes the pinned packages in `.codex/runtime-requirements.txt`. Later runs reuse the same verified runtime.

NumPy, pandas, Matplotlib, Plotly, SciPy, scikit-learn, and seaborn are required runtime components, not optional extras. Do not skip an analysis or promised visualization because one is missing: repair the runtime first and rerun. Only report a dependency failure after bootstrap/repair itself is blocked by network, sandbox, or filesystem policy.

PLINK is the exception: never install or upgrade PLINK automatically. Population-structure analysis resolves the explicitly configured/bundled PLINK executable separately.
