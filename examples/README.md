# Example input files

This directory contains only the minimal example material needed to understand the recommended input formats for Geng-BEVA skills.

The files here are **examples**, not production reference resources. The two production reference datasets bundled with the skills are:

- `skills/functional-gene-diagnostic/assets/RiceNavi_subset.csv` — functional-gene annotation resource.
- `skills/igci-evaluation/assets/IJ_panel.csv` — japonica/indica diagnostic allele panel used by IGCI.

## `example_geng_beva.vcf`

**Not bundled in the repository. Add your own small, publishable VCF at this path when you want a repository demo:**

```text
.codex/examples/example_geng_beva.vcf
```

### Recommended use

A small real Geng-BEVA VCF can be used to demonstrate:

- PCA and IBS/IBD in `population-structure-analysis`;
- per-site or requested sliding-window FST when suitable sample groups are defined;
- IGCI and single-sample indica introgression profiles;
- functional-gene diagnostics;
- transgene-assay evidence evaluation when the corresponding project-specific assay records are present.

### Recommended format

Use a standard diploid VCF with sample genotype fields. `GT` is required for genotype-based analyses; `GT:DP` is recommended when depth-aware analyses are needed.

Minimal structure:

```text
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample01\tSample02
1\t100000\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1
```

VCF `REF` and `ALT` describe only how alleles are encoded in that VCF. They must not be interpreted automatically as Nipponbare/non-Nipponbare or japonica/indica alleles. The functional-gene and IGCI workflows reorient alleles against their bundled biological reference resources before interpretation.

The example VCF may be located elsewhere. `examples/example_geng_beva.vcf` is only the recommended repository demo location; normal analyses should use the exact path supplied by the user.

## `example_groups.csv`

### Used by

- FST and other explicit population comparisons when group membership is supplied as a file.

### When to use

Use a group file when sample membership is not stated directly in the user prompt or when a reusable grouping definition is convenient.

### Recommended format

```csv
Sample,Group
Sample01,GroupA
Sample02,GroupA
Sample03,GroupB
Sample04,GroupB
```

Rules:

- `Sample` values must exactly match sample IDs in the VCF header.
- `Group` contains the population/group label used for comparison.
- Replace the example sample names and groups with the actual samples in your VCF.

## Files intentionally not included

No test fixtures, synthetic VCFs, fake reference panels, preference CSVs, or duplicate copies of production resources are included. User-specific breeding objectives should normally be stated in the prompt rather than stored as a required CSV.
