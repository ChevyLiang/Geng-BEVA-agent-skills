---
name: population-structure-analysis
description: Analyze population structure and genetic relationships in standardized Geng-BEVA rice VCF data using PCA, PLINK IBS/IBD, pairwise comparisons, genotype summaries, and Weir-Cockerham FST with Manhattan-style genome scans. Use for Geng-BEVA genetic structure, similarity, differentiation, or sample-comparison tasks; use transgene-detection for assay-element evidence.
---

# Population structure analysis

On Windows Codex, invoke this skill through `.codex/run-geng-beva.ps1 population-structure-analysis ...` so the fixed runtime is verified before analysis. Direct Python entry points remain available for development after the runtime is ready. For genomic analyses, first work from canonical rice chromosomes 1-12. This preprocessing step removes project-specific assay records but does not classify transgene evidence.

Use the sibling `transgene-detection` skill when the task is about BT, 35S, CP4_EPSPS, or other transgene assay elements.

## Analysis routing

- `pca`: PLINK PCA coordinates plus **required** `pca_interactive.html`. Hovering a point shows the sample name (and group/PC coordinates when available). A static PNG is secondary. Do not report PCA as complete if the interactive HTML is missing; repair the fixed runtime and rerun.
- `ibs`: PLINK `DST` allele-sharing coefficient matrix plus **required** `ibs_heatmap_interactive.html`. Hovering a cell shows Sample A, Sample B, and the IBS value. Do not report IBS as complete if this HTML is missing; repair the fixed runtime and rerun.
- `ibd`: PLINK `PI_HAT` estimated IBD-proportion matrix plus **required** `ibd_heatmap_interactive.html`.
- `pairwise-compare`: PCA distance, IBS, IBD, and differing called-genotype loci for two samples.
- `pairwise-matrix`: PCA distances and IBS matrix for selected samples.
- `fst`: per-site Weir-Cockerham FST (theta) for two requested groups, with a Manhattan-style genome scan. Per-site FST is the default. Only calculate sliding-window FST when the user requests it; when requested without parameters, use a 500-kb window and 100-kb step. User-provided window/step values override these defaults.
- `genotype-frequency`: genotype-count/frequency summaries.
- `genotype-summary`: per-sample call, heterozygosity, homozygous-alt, missing, and alternate-allele dosage summaries.
- `sample-names`: sample IDs from the supplied VCF header.

Keep IBS and IBD distinct; never calculate IBD as `1 - IBS`. Exclude missing genotypes from pairwise genotype-difference counts. Retain negative finite Weir-Cockerham FST estimates rather than silently truncating them.

`AltAlleleDosageSum = HeterozygousCount + 2 * HomoAltCount` is a genotype-summary statistic and is **not nucleotide diversity (pi)**. The CLI alias `diversity` is retained only for backward compatibility and should not be described as nucleotide diversity.

Do not silently add LD pruning, MAF filtering, or other population-genetic preprocessing unless requested by the user or required by a documented command.

## Runtime

Use the exact input supplied by the user. The workflow targets standardized diploid, biallelic Geng-BEVA VCFs. Under the project convention, non-chromosome-1-12 data rows may be assay records; do not generalize this rule to unrelated VCFs. Files under `examples/` are format examples only.

PLINK 1.9 is required for PCA/IBS/IBD. On Windows, prefer the skill-local executable at `tools/plink/windows-x64/plink.exe`. Resolution order is: explicit `--plink-cmd`, the skill-local bundled PLINK executable, `PLINK_BIN`, then `PATH`. Do **not** auto-install PLINK, write to the user profile, or create a `~/.geng-beva` tool cache. If no usable executable is found, report the expected skill-local path and stop.

Treat the VCF path supplied in the user prompt as the biological input. Do not require VCF files to live inside the repository. `examples/example_geng_beva.vcf` is only the recommended optional demo location.

Read:
- `references/usage.md` for CLI commands.
- `references/plink.md` for executable setup and IBS/IBD fields.
- `references/contract.md` for the Python/manifest contract.

## Fixed runtime and self-healing

On Windows Codex, use the workspace entry point `.codex/run-geng-beva.ps1`; do not treat the active system/Conda Python as the authoritative environment. Before analysis, the entry point verifies `.codex/runtime/geng-beva/.venv`. If the runtime is missing or damaged, it bootstraps/repairs the managed Python 3.11 environment and synchronizes the pinned packages in `.codex/runtime-requirements.txt`. Later runs reuse the same verified runtime.

NumPy, pandas, Matplotlib, Plotly, SciPy, scikit-learn, and seaborn are required runtime components, not optional extras. Do not skip an analysis or promised visualization because one is missing: repair the runtime first and rerun. Only report a dependency failure after bootstrap/repair itself is blocked by network, sandbox, or filesystem policy.

PLINK is the exception: never install or upgrade PLINK automatically. Population-structure analysis resolves the explicitly configured/bundled PLINK executable separately.
