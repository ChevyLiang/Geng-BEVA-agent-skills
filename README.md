# Codex Rice Analysis Skills

## A modular analysis framework for rice functional genomics, introgression, and population structure

This repository provides a set of Codex-based analytical skills for **Shennongxin-1 Toolkit (Geng-BEVA)**, enabling reproducible and scalable analysis of:

- Functional genetic variation (QTN and trait-associated loci)
- Indica–Japonica genomic introgression patterns
- Population genetic structure and relatedness

All modules are designed around standardized VCF inputs and curated biological reference resources, allowing consistent application across diverse rice germplasm datasets.

---

## Abstract

High-throughput sequencing has enabled large-scale characterization of rice genetic diversity; however, translating genomic variation into interpretable breeding-relevant knowledge remains challenging. Here, we present a modular analytical framework composed of three interoperable Codex skills that integrate functional annotation, introgression inference, and population structure analysis. The system leverages curated prior knowledge (Rice_Navi and IJ reference panels) together with user-provided VCF data to generate interpretable summaries for breeding decision support.

---

## System Overview

| Module | Function | Biological Focus | Output |
|---|---|---|---|
| functional-gene-diagnostic | Functional variant interpretation | QTN / trait-associated loci | Trait-level functional profiles |
| igci-evaluation | Introgression quantification | Indica–Japonica differentiation | Genome- and chromosome-level ancestry proportions |
| population_structure_agent | Genetic structure inference | Population stratification | PCA, IBS/IBD, clustering |

---

## Reference Data Resources

Two curated datasets are included and must remain unchanged:

### Rice_Navi.csv
A curated functional annotation resource containing:
- experimentally supported or high-confidence QTN
- trait-associated loci

### IJ_panel.csv
A fixed diagnostic panel for indica–japonica differentiation:
- Ref allele: japonica state  
- Alt allele: indica state  

These resources serve as biological priors for downstream inference.

---

## Input Data Format

All modules operate on standard VCF (Variant Call Format) files.

A minimal test dataset is provided:

.codex/test_data/test.vcf

Users should replace this file with their own dataset for analysis.

---

## Installation and Execution

After cloning the repository:

```bash
cd <project_root_above_.codex>
codex
```

Inside Codex:

```bash
$ functional-gene-diagnostic
$ igci-evaluation
$ population_structure_agent
```

Each command triggers an independent analysis pipeline operating on the active VCF file.

---

## Data Usage Guidelines

### Single VCF constraint
Only one VCF file should exist in `.codex/test_data/` to ensure deterministic behavior.

### Sample naming convention (strongly recommended)

Recommended:
- YF47
- CJCx37
- IR64_01
- Japonica_12

Not recommended:
- Non-English characters
- Spaces or mixed encoding
- Chinese sample names

---

## Modules

### 1. Functional Gene Diagnostic
Identifies functional alleles and QTN distributions across samples or populations.

Outputs:
- functional variant matrices
- trait/category summaries
- gene-level annotation tables

---

### 2. IGCI Evaluation
Estimates indica–japonica introgression using diagnostic alleles:

- Ref → japonica
- Alt → indica

Outputs:
- sample introgression scores
- chromosome-level profiles
- genome-wide summaries
- heatmaps (PNG/PDF)

---

### 3. Population Structure Analysis
Performs:
- PCA
- IBS/IBD similarity
- clustering

Outputs:
- PCA coordinates
- similarity matrices
- population clusters

---

## Outputs

All results are saved in:

outputs/

Key fields:
- SampleID
- Trait
- Category
- Gene
- LocusID

---

## Reproducibility

- Deterministic given identical VCF + reference data
- Functional inference depends on Rice_Navi.csv
- Introgression inference depends on IJ_panel.csv
- Population structure derived from genotype matrix

---


