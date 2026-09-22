
## Windows Codex runtime

Prefer `.codex/run-geng-beva.ps1` as the execution entry point. It verifies or repairs the fixed Geng-BEVA Python 3.11 runtime before dispatching the skill. Do not consider a requested output complete when a required Python package is missing; repair the runtime and rerun. PLINK remains separately managed and is never auto-installed.

# Usage

Input VCF and panel files may be located anywhere accessible to the runtime. The paths below are examples only.

```bash
python skills/igci-evaluation/scripts/igci_evaluation.py --vcf /path/to/input.vcf --out-dir outputs/igci
python skills/igci-evaluation/scripts/igci_evaluation.py --vcf /path/to/input.vcf --sample SAMPLE_A --out-dir outputs/igci_A

# Override the bundled production panel only when needed
python skills/igci-evaluation/scripts/igci_evaluation.py --vcf /path/to/input.vcf --panel /path/to/custom_IJ_panel.csv --out-dir outputs/igci
```

Default regional analysis uses `--window 500000 --step 100000`. These parameters can be changed explicitly, and `--no-plot` suppresses figures.

When `--sample SAMPLE_A` is supplied, the workflow produces a **single-sample indica introgression profile** across rice chromosomes 1-12. Each 500-kb window is colored by its mean indica-allele contribution score and windows advance every 100 kb. Gray windows have no informative markers. The PNG and PDF profile are intended to show where indica-derived diagnostic alleles are concentrated along the genome; they should not be described as exact physical ancestry boundaries.

Bundled production panel schema:

```text
chr,POS,Geng,Xian
```

`Geng` is the japonica diagnostic allele and `Xian` is the indica diagnostic allele. These panel alleles define ancestry; VCF REF/ALT only define how that particular VCF encodes the alleles. The workflow always aligns the VCF alleles to Geng/Xian before computing dosage, including when REF/ALT are reversed. The bundled panel is `assets/IJ_panel.csv` and contains 16,182 diagnostic sites across chromosomes 1-12. If `--panel` is omitted, this bundled panel is used automatically. The legacy `chr,POS,Ref,Alt` schema is accepted only for backward compatibility. The bundled production panel is always the default reference panel unless the user explicitly supplies another panel.

Per-sample outputs include site scores, chromosome summary, genome summary, sliding windows, PNG, and PDF. Cohort outputs include all-sample score and summary tables. The manifest records explicit output paths and the window/step used.
