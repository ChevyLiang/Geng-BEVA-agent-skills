# Geng-BEVA Skills

Codex skills and command-line workflows for analyzing standardized **Geng-BEVA rice genotype data**.

This repository packages four focused skills for researchers who use Geng-BEVA data with Codex: functional-gene diagnostics, indica genomic contribution evaluation, population-structure analysis, and transgene-assay evidence evaluation. Each skill follows the Codex skill layout (`SKILL.md`, `agents/`, and optional `scripts/`, `references/`, and `assets/`), while the suite shares a common managed runtime and launcher.

> **Scope.** These workflows are designed for the Geng-BEVA data conventions documented in each skill. They should not be applied blindly to unrelated VCF files.

## Skills

| Skill | Purpose | PLINK 1.9 |
| --- | --- | --- |
| `functional-gene-diagnostic` | Annotate functional loci using a RiceNavi-derived QTN subset and summarize breeding-relevant alleles | No external tool required |
| `igci-evaluation` | Calculate the Indica Genomic Contribution Index (IGCI) from a diagnostic indica/japonica allele panel | No |
| `population-structure-analysis` | PCA, IBS/IBD, pairwise comparison, genotype summaries, and Weir-Cockerham FST | Required for PCA/IBS/IBD |
| `transgene-detection` | Evaluate transgene-assay records using sample GT/DP evidence tiers | No |

## Install for Codex

Install/deploy the **entire `.codex` suite as one unit**. The four skills share `tools/`, `runtime-requirements.txt`, and the managed runtime bootstrap; copying only an individual skill folder is not a complete installation.

For a workspace-local deployment, the expected structure is:

```text
<workspace>/
└── .codex/
    ├── run-geng-beva.ps1
    ├── runtime-requirements.txt
    ├── tools/
    ├── runtime/
    └── skills/
        ├── functional-gene-diagnostic/
        ├── igci-evaluation/
        ├── population-structure-analysis/
        └── transgene-detection/
```

Restart/reload Codex after replacing the suite. Skills can still be selected by their normal names, but on Windows their execution must use the shared fixed-runtime entry point described below.

## Fixed Windows runtime

For Windows Codex, Geng-BEVA uses its own persistent managed runtime rather than assuming the user already has Python or scientific packages installed. Use the canonical entry point:

```powershell
powershell -ExecutionPolicy Bypass -File .codex\run-geng-beva.ps1 population-structure-analysis pca input.vcf --out-dir outputs\pca
```

On first use, the bootstrap reuses `uv` if available or downloads a fixed workspace-local copy, obtains managed Python 3.11 when needed, creates `.codex/runtime/geng-beva/.venv`, synchronizes the exact packages in `runtime-requirements.txt`, verifies them, and records a runtime manifest. Subsequent runs verify and reuse the same environment; damaged or incomplete runtimes are repaired before analysis.

The required stack includes NumPy, pandas, Matplotlib, Plotly, SciPy, scikit-learn, and seaborn. Plotly is not optional: interactive PCA and IBS/IBD HTML files are standard outputs. Missing dependencies must trigger runtime repair rather than silent output degradation. See `RUNTIME.md`.

### PLINK 1.9

PCA, IBS, and IBD use genomic PLINK 1.9. For the primary Windows Codex workflow, place your own `plink.exe` at:

```text
skills/population-structure-analysis/tools/plink/windows-x64/plink.exe
```

The binary is intentionally **not bundled**. The resolver checks, in order: `--plink-cmd`, the skill-local Windows executable, `PLINK_BIN`, and then `PATH`. It validates that the executable is genomic PLINK 1.9. The skill does **not** auto-install PLINK and does not write a tool cache under the Windows user profile.

You can alternatively point to an existing installation:

```powershell
$env:PLINK_BIN = "C:\\path\\to\\plink.exe"
```

### Windows/Codex path robustness

Functional-gene diagnostics use a pure-Python streaming VCF backend and do not require PLINK. If a Windows Codex sandbox reports `CreateProcessWithLogonW failed: 267`, run the existing skill entry point from an ASCII-only working directory such as `%TEMP%\\geng-beva-run` and pass the Python interpreter, skill path, VCF path, and output directory as absolute paths. Do not replace the workflow with an ad-hoc heredoc/piped Python script.

The transgene convention is project-specific and must not be generalized to arbitrary VCF files.

### IGCI panel

IGCI ships with the production diagnostic panel at `skills/igci-evaluation/assets/IJ_panel.csv`. Its preferred schema is:

```text
chr,POS,Geng,Xian
```

`Geng` is the confirmed japonica (geng) diagnostic allele and `Xian` is the confirmed indica (xian) diagnostic allele. The bundled panel contains 16,182 diagnostic sites across rice chromosomes 1-12. If `--panel` is omitted, the workflow uses this bundled panel automatically. The analysis handles VCF REF/ALT reversal relative to these diagnostic alleles and excludes missing or unalignable genotypes from each sample's denominator.

For backward compatibility, custom panels using `chr,POS,Ref,Alt` are also accepted, with `Ref` interpreted as japonica and `Alt` as indica.

## Method notes

### Functional-gene allele correction and default favorability

The functional-gene workflow does **not** assume that VCF `REF`/`ALT` correspond to Nipponbare/non-Nipponbare alleles. It first aligns literal VCF alleles to the bundled RiceNavi biological alleles (`Nipponbare_Allele` and `Alternative_Allele`) and only then interprets the sample genotype. This prevents frequency-oriented or otherwise reversed VCF encoding from changing the biological conclusion.

Unless the user explicitly overrides the breeding objective, resolvable homozygous functional genotypes are classified as complementary favorable/unfavorable pairs. `Alternative_Allele_Effect=Positive` means Alternative favorable / Nipponbare unfavorable; `Negative` means Nipponbare favorable / Alternative unfavorable. Earlier heading is favorable by default, so heading `promote` means Alternative favorable and heading `delay` means Nipponbare favorable. Heterozygotes are reported separately. Non-heading `Uncertainty` remains unresolved unless the user supplies an explicit preference.

Single-sample analysis writes two compact breeding-facing result files: `favorable_genotypes.csv` (favorable homozygotes plus heterozygous calls) and `unfavorable_genotypes.csv` (unfavorable homozygotes). Both contain `GeneName`, `Position`, `NipponbareGenotype`, `AlternativeGenotype`, `SampleGenotype`, `GenotypeBackground`, `PhenotypicFunction`, and `GenotypeEffect`. Heterozygotes are explicitly labeled `Heterozygous` and are not treated as homozygous favorable.

### IGCI

For each informative diploid site, indica allele dosage is normalized to 0, 0.5, or 1. For a sample:

```text
IGCI = sum(indica allele dosage / 2) / number of informative sites
```

Genome-wide IGCI pools informative sites across chromosomes. IGCI is a diagnostic-allele contribution index; it is not a calibrated ancestry probability and is not identical to the physical fraction of introgressed sequence.

For one requested sample, the default 500-kb/100-kb sliding-window analysis also generates a chromosome-level **indica introgression profile**. This visualizes the local concentration of indica diagnostic alleles across chromosomes 1-12; it should not be interpreted as an exact sequence-ancestry boundary map.

### Population differentiation

The `fst` workflow calculates per-site **Weir-Cockerham FST (theta)** by default and writes a Manhattan-style whole-genome scatter plot. Sliding-window FST is optional and is only calculated when requested. If the user requests sliding windows without specifying parameters, the defaults are a 500-kb window and a 100-kb step; user-provided values take precedence. Window FST is aggregated from Weir-Cockerham variance components (`sum(a) / sum(a+b+c)`). Negative finite estimates are retained rather than silently truncated to zero. The estimator follows Weir BS & Cockerham CC (1984), *Evolution* 38:1358-1370.

### Individual genotype summary

The former `DiversityScore` label has been replaced by the explicit `AltAlleleDosageSum` (`heterozygous count + 2 x homozygous-alt count`). This is a genotype-summary statistic, **not nucleotide diversity (pi)**. The legacy CLI alias `diversity` remains available for compatibility.

### Transgene assay evidence

The transgene workflow evaluates explicitly encoded assay records. Its default DP cutoffs (10, 50, and 100) are **heuristic evidence tiers**, not validated biological detection thresholds. Low depth is not proof of absence, and the workflow does not infer insertion position, copy number, sensitivity, specificity, or inheritance.

## RiceNavi-derived annotation subset

`skills/functional-gene-diagnostic/assets/RiceNavi_subset.csv` is a selected subset derived from the RiceNavi resource developed by Xuehui Huang and colleagues. It is included for functional-locus annotation and must not be presented as an original Geng-BEVA database.

Please cite the original RiceNavi study when this resource contributes to an analysis:

> Wei X, Qiu J, Yong K, et al. A quantitative genomics map of rice provides genetic insights and guides breeding. *Nature Genetics*. 2021;53:243-253. doi:10.1038/s41588-020-00769-9.

See `THIRD_PARTY_NOTICES.md` for attribution and scope.

## Examples

A minimal example area is provided at `examples/`. It intentionally contains only format guidance and one small group-file example.

- `examples/README.md` explains when each example input is used and the recommended format.
- `examples/example_groups.csv` demonstrates the recommended sample-to-group mapping for FST/population comparisons.
- Add your own small publishable VCF at `examples/example_geng_beva.vcf` when you want a repository demo. The VCF is not bundled here.

Normal analyses are not restricted to the `examples/` directory: use the exact VCF path supplied by the user.

Production reference resources are kept with their skills, not in `examples/`:

- `skills/functional-gene-diagnostic/assets/RiceNavi_subset.csv`
- `skills/igci-evaluation/assets/IJ_panel.csv`

No test fixtures, fake VCFs, duplicate reference panels, or required preference CSVs are distributed in the public suite.

## Citation

If you use this software, cite the repository using `CITATION.cff`. When the Geng-BEVA manuscript becomes available, add the manuscript citation/DOI alongside the software citation.

For analyses using the RiceNavi-derived subset, also cite Wei et al. (2021) as listed above.

## Authors

**Chevy Liang & Jian Sun**  
Rice Research Institute, Shenyang Agricultural University  
Contact: <chevyliang0905@163.com>

## License

Code in this repository is released under the MIT License. Third-party data/resources retain their original attribution and applicable terms; see `THIRD_PARTY_NOTICES.md`.
