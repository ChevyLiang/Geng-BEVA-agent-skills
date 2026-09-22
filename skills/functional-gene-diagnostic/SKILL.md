---
name: functional-gene-diagnostic
description: Annotate Geng-BEVA rice VCF genotypes against the bundled RiceNavi-derived functional-locus subset, correct VCF REF/ALT into RiceNavi Nipponbare/alternative allele coordinates, identify default favorable alleles (RiceNavi Positive and earlier-heading rules), summarize samples/populations/traits, and identify breeding targets. Use for Geng-BEVA functional-marker interpretation, not GWAS, genomic selection, or causal discovery.
---

# Functional gene diagnostic

Use `assets/RiceNavi_subset.csv` as the default annotation table. It is a RiceNavi-derived subset and must retain RiceNavi attribution.

## Critical allele semantics

Never interpret VCF `REF` as Nipponbare or VCF `ALT` as non-Nipponbare by position alone. In user VCFs, REF/ALT can reflect the encoding chosen during variant construction and may be reversed relative to RiceNavi biological alleles.

The bundled CSV uses these explicit biological columns:

- `Nipponbare_Allele`: the RiceNavi Nipponbare/reference allele.
- `Alternative_Allele`: the non-Nipponbare alternative allele.
- `Alternative_Allele_Function`: phenotype/function attributed to the Alternative allele.
- `Alternative_Allele_Direction`: normalized direction token for the Alternative allele.
- `Alternative_Allele_Effect`: RiceNavi effect label for the Alternative allele (`Positive`, `Negative`, `Uncertainty`).

Before interpreting a genotype, align `VCFRefAllele/VCFAltAllele` to `Nipponbare_Allele/Alternative_Allele`. Only then assign `RiceNaviGenotypeState` (`NIPPONBARE_HOM`, `ALTERNATIVE_HOM`, `HET`, `MISSING`, `AMBIGUOUS`).

## Default favorable-allele rules

Unless the user explicitly supplies a different breeding preference, apply these project defaults in this order:

1. **Earlier heading is favorable.** For heading/flowering loci, if the Alternative allele promotes/advances heading, the Alternative allele is favorable. If the Alternative allele delays heading, the Nipponbare allele is favorable.
2. **RiceNavi Positive defines a binary favorable pair.** For other traits, `Alternative_Allele_Effect=Positive` means Alternative is favorable and Nipponbare is unfavorable.
3. **RiceNavi Negative defines the opposite binary pair.** `Alternative_Allele_Effect=Negative` means Nipponbare is favorable and Alternative is unfavorable.
4. Non-heading `Uncertainty` remains unresolved unless an explicit user preference resolves it.
5. An explicit user breeding-preference record overrides the defaults above.

For every resolvable locus, the two homozygous states are complementary: one is favorable and the other is unfavorable. Do not infer dominance for heterozygotes; report them as a separate `HETEROZYGOUS` class.

## Workflow

1. Resolve requested VCF sample IDs.
2. Stream the VCF directly in Python and retain only bundled functional loci; do not launch PLINK for the default functional-diagnostic path.
3. Align each VCF REF/ALT pair to the RiceNavi Nipponbare/Alternative alleles.
4. Convert each sample call into `RiceNaviGenotypeState` and `AlternativeAlleleDosage`.
5. Apply the default favorable rules above, unless an explicit preference overrides them.
6. Preserve every locus-trait annotation; do not collapse multiple trait records at one locus.
7. Classify resolvable called genotypes into five biological classes: `NIPPONBARE_FAVORABLE`, `NIPPONBARE_UNFAVORABLE`, `ALTERNATIVE_FAVORABLE`, `ALTERNATIVE_UNFAVORABLE`, or `HETEROZYGOUS`. Keep missing/unalignable/true-Uncertainty records as `UNRESOLVED` QC records.
8. Report all annotations plus separate class-specific outputs.

For single-sample mode, produce only two primary result CSV files:

- `favorable_genotypes.csv`: favorable homozygous genotypes plus heterozygous calls. Heterozygotes must be labeled `Heterozygous` and must not be described as homozygous favorable.
- `unfavorable_genotypes.csv`: unfavorable homozygous genotypes only.

Both files use the same compact eight-column schema:

- `GeneName`: RiceNavi gene name.
- `Position`: chromosome and physical position (`ChrN:POS`).
- `NipponbareGenotype`: homozygous Nipponbare genotype using literal biological alleles.
- `AlternativeGenotype`: homozygous non-Nipponbare genotype using literal biological alleles.
- `SampleGenotype`: the sample genotype after VCF REF/ALT has been re-oriented to RiceNavi biological alleles.
- `GenotypeBackground`: `Nipponbare`, `Alternative`, or `Heterozygous`, describing which RiceNavi biological background the current call belongs to.
- `PhenotypicFunction`: trait plus the annotated function of the non-Nipponbare allele.
- `GenotypeEffect`: `Favorable`, `Unfavorable`, or `Heterozygous`.

Do not expose VCF REF/ALT, 0/1 genotype coding, alignment status, dosage, rule-source fields, or other internal diagnostic columns in these primary result files. Unresolved records are omitted from the two CSVs and counted in the manifest warning/metrics.

## Interpretation safeguards

- `FavorableAlleleGroup` is in RiceNavi biological coordinates: `NIPPONBARE`, `ALTERNATIVE`, or `UNKNOWN`.
- `FavorableVCFAlleleGroup` is in the current VCF encoding: `REF`, `ALT`, or `UNKNOWN`.
- `FavorableAllele` is the literal nucleotide/indel sequence.
- `FavorableRule` records why the allele was selected (`RiceNavi_positive`, `RiceNavi_negative`, `early_heading_alternative`, `early_heading_nipponbare`, `explicit_preference`, or `unresolved`).
- `UnfavorableAlleleGroup` is always the opposite homozygous biological group when favorability is resolved.
- `AdvantageClass` is the final five-class genotype interpretation; `UNRESOLVED` is reserved for QC-only cases.
- Missing or unalignable calls are never classified as unfavorable.
- Do not call a site favorable solely because the VCF genotype is ALT/ALT.
- Do not treat `Alternative_Allele_Effect=Positive` as a generic increase/decrease direction; the phenotype direction is read separately from `Alternative_Allele_Function` and `Alternative_Allele_Direction`.

Read `references/dataflow.md` for field definitions and `references/breeding_preferences.md` before preference-based overrides.

## Runtime

Use the exact input supplied by the user. Files under `examples/` are only for format demonstration and must not replace the user's biological input unless the user explicitly requests a demo.

The workflow targets standardized diploid, biallelic Geng-BEVA VCF data. Chromosomes 1-12 are canonical rice genomic records under this project convention.

Resolve relative paths from the caller's working directory. Use a separate output directory per run and preserve source inputs. Functional-gene diagnosis is zero-external-process by default and does not require PLINK. `--plink-cmd` is retained only for backward compatibility.

On Windows/Codex, invoke the existing `main.py` rather than generating ad-hoc analysis scripts. If `CreateProcessWithLogonW failed: 267` occurs, use an ASCII-only process working directory (for example `%TEMP%\geng-beva-run`) and pass the Python interpreter, skill path, biological input, and output paths as absolute paths.

Read:
- `references/usage.md` for CLI modes.
- `references/dataflow.md` for allele-coordinate semantics and outputs.
- `references/breeding_preferences.md` for override precedence.
- `references/contract.md` for the Python `run(...)` and manifest contract.

## Fixed runtime and self-healing

On Windows Codex, use the workspace entry point `.codex/run-geng-beva.ps1`; do not treat the active system/Conda Python as the authoritative environment. Before analysis, the entry point verifies `.codex/runtime/geng-beva/.venv`. If the runtime is missing or damaged, it bootstraps/repairs the managed Python 3.11 environment and synchronizes the pinned packages in `.codex/runtime-requirements.txt`. Later runs reuse the same verified runtime.

NumPy, pandas, Matplotlib, Plotly, SciPy, scikit-learn, and seaborn are required runtime components, not optional extras. Do not skip an analysis or promised visualization because one is missing: repair the runtime first and rerun. Only report a dependency failure after bootstrap/repair itself is blocked by network, sandbox, or filesystem policy.

PLINK is the exception: never install or upgrade PLINK automatically. Population-structure analysis resolves the explicitly configured/bundled PLINK executable separately.
