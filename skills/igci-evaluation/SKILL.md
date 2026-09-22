---
name: igci-evaluation
description: Calculate the Indica Genomic Contribution Index (IGCI) from standardized Geng-BEVA rice VCF data and an indica/japonica diagnostic-allele panel, producing site scores, chromosome/genome summaries, sliding-window tables, and single-sample introgression profiles. Use when evaluating indica allele contribution in Geng-BEVA materials.
---

# IGCI evaluation

Use the bundled production `assets/IJ_panel.csv` by default. Its schema is `chr,POS,Geng,Xian`, where `Geng` is the japonica diagnostic allele and `Xian` is the indica diagnostic allele. Accept an explicit `--panel` path when the user wants to override the bundled panel; also accept the legacy `chr,POS,Ref,Alt` schema for backward compatibility, with `Ref` interpreted as japonica and `Alt` as indica.

For each called, allele-aligned diploid site, calculate indica allele dosage / 2:
- japonica homozygote: 0;
- heterozygote: 0.5;
- indica homozygote: 1.

Treat the panel as the biological authority and the VCF REF/ALT fields only as encoding. For every diagnostic site, first align VCF REF/ALT to the panel-defined `Geng`/`Xian` alleles, then convert the sample GT to Xian allele dosage. Handle exact and reversed REF/ALT orientation explicitly. Never assume VCF ALT is indica. Exclude missing, partially missing, multiallelic, or unalignable records from the sample-specific denominator.

The per-site output must retain `GENG_ALLELE`, `XIAN_ALLELE`, `VCF_REF`, `VCF_ALT`, `VCF_GT`, `ALIGNMENT`, `ANCESTRY_GENOTYPE`, `XIAN_ALLELE_DOSAGE`, and `SCORE` so that every IGCI contribution is auditable.

Report chromosome and genome values as:

`IGCI = SUPPORT_SUM / INFORMATIVE_SITES`

Pool all informative sites for genome-wide IGCI; do not average chromosome IGCI values without weighting. Interpret IGCI as a diagnostic-allele contribution index, not as a calibrated ancestry probability or an exact physical fraction of introgressed bases.

No informative sites means missing IGCI, never zero. Empty sliding windows remain missing and are displayed as gray. Unless the user requests otherwise, use a 500-kb window and 100-kb step. When the user asks where indica contribution occurs in one accession, run a single-sample analysis and generate the chromosome-level introgression profile from these windows. Describe it as a diagnostic-allele contribution profile, not an exact ancestry-boundary map.

Use the exact input path supplied by the user, regardless of where the VCF is stored. Do not require production VCFs to be copied into the repository. The workflow targets standardized, diploid, biallelic Geng-BEVA VCFs and does not require PLINK. The optional `examples/` directory is for format demonstration only.

Do not ask the user to supply a panel unless they explicitly want a different one or the bundled production panel is unavailable.

Read:
- `references/igci_reference.md` for definitions.
- `references/usage.md` for CLI usage and panel schema.
- `references/contract.md` for the Python `run(...)` and manifest contract.

## Fixed runtime and self-healing

On Windows Codex, use the workspace entry point `.codex/run-geng-beva.ps1`; do not treat the active system/Conda Python as the authoritative environment. Before analysis, the entry point verifies `.codex/runtime/geng-beva/.venv`. If the runtime is missing or damaged, it bootstraps/repairs the managed Python 3.11 environment and synchronizes the pinned packages in `.codex/runtime-requirements.txt`. Later runs reuse the same verified runtime.

NumPy, pandas, Matplotlib, Plotly, SciPy, scikit-learn, and seaborn are required runtime components, not optional extras. Do not skip an analysis or promised visualization because one is missing: repair the runtime first and rerun. Only report a dependency failure after bootstrap/repair itself is blocked by network, sandbox, or filesystem policy.

PLINK is the exception: never install or upgrade PLINK automatically. Population-structure analysis resolves the explicitly configured/bundled PLINK executable separately.
