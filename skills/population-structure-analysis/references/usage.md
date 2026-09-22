
## Windows Codex runtime

Prefer `.codex/run-geng-beva.ps1` as the execution entry point. It verifies or repairs the fixed Geng-BEVA Python 3.11 runtime before dispatching the skill. Do not consider a requested output complete when a required Python package is missing; repair the runtime and rerun. PLINK remains separately managed and is never auto-installed.

# Usage

Input VCF files may be located anywhere accessible to the runtime. The paths below are examples only; the user does not need to copy biological data into this repository.

```bash
python skills/population-structure-analysis/main.py sample-names /path/to/input.vcf --out-dir outputs/names
python skills/population-structure-analysis/main.py pca /path/to/input.vcf --group-file /path/to/groups.csv --sample-id SAMPLE_A --out-dir outputs/pca
python skills/population-structure-analysis/main.py ibs /path/to/input.vcf --out-dir outputs/ibs
python skills/population-structure-analysis/main.py ibd /path/to/input.vcf --out-dir outputs/ibd
python skills/population-structure-analysis/main.py pairwise-compare /path/to/input.vcf --sample-a SAMPLE_A --sample-b SAMPLE_B --out-dir outputs/pair
python skills/population-structure-analysis/main.py pairwise-matrix /path/to/input.vcf --samples SAMPLE_A SAMPLE_B SAMPLE_C --out-dir outputs/matrix
python skills/population-structure-analysis/main.py fst /path/to/input.vcf --group-a SAMPLE_A SAMPLE_B --group-b SAMPLE_C SAMPLE_D --out-dir outputs/fst
python skills/population-structure-analysis/main.py fst /path/to/input.vcf --group-a SAMPLE_A SAMPLE_B --group-b SAMPLE_C SAMPLE_D --sliding-window --out-dir outputs/fst_window
python skills/population-structure-analysis/main.py fst /path/to/input.vcf --group-a SAMPLE_A SAMPLE_B --group-b SAMPLE_C SAMPLE_D --sliding-window --window 1000000 --step 200000 --out-dir outputs/fst_custom_window
python skills/population-structure-analysis/main.py genotype-summary /path/to/input.vcf --out-dir outputs/genotype_summary
```

The FST workflow writes per-site Weir-Cockerham estimates and a Manhattan-style genome scan by default. It does not calculate sliding windows unless requested. Use `--sliding-window` to enable windowed FST; if `--window`/`--step` are omitted, the defaults are 500 kb and 100 kb. Supplying `--window` or `--step` also implies sliding-window mode. User-provided values take precedence. Use `--min-sites` to require a minimum number of valid variants in each window.

Additional commands: `genotype-frequency` and `nontransgene-vcf`.

For backward compatibility, `diversity` aliases `genotype-summary`. The legacy `transgene-events` and `transgene-filter` entry points remain compatibility bridges; new transgene analyses should use the sibling `transgene-detection` skill directly.

Group file schema: `SampleID,Group`. Use exact VCF sample IDs.
## Interactive-output dependency behavior

PCA and IBS/IBD interactive HTML outputs require Plotly. If Plotly is absent, the skill must invoke the shared dependency resolver and install Plotly into the persistent Geng-BEVA runtime automatically. Do not tell the user to install Plotly manually unless automatic recovery fails. The installed package is reused on subsequent runs.

