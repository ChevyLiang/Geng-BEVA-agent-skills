
## Windows Codex runtime

Prefer `.codex/run-geng-beva.ps1` as the execution entry point. It verifies or repairs the fixed Geng-BEVA Python 3.11 runtime before dispatching the skill. Do not consider a requested output complete when a required Python package is missing; repair the runtime and rerun. PLINK remains separately managed and is never auto-installed.

# Usage

```bash
python skills/transgene-detection/main.py data/input.vcf --out-dir outputs/transgene
python skills/transgene-detection/main.py data/input.vcf --sample-id SAMPLE_A --export-clean --out-dir outputs/transgene_A
```

Default heuristic depth thresholds can be adjusted:

```bash
python skills/transgene-detection/main.py data/input.vcf \
  --low-dp 10 --moderate-dp 50 --high-dp 100 \
  --out-dir outputs/transgene
```

Interpretation with the defaults is: low support at DP >= 10, moderate support at DP >= 50, and high support at DP > 100. The high threshold is implemented as a strict `>` cutoff to preserve the existing convention.

Primary output: `transgene_sample_summary.csv`. Supporting output: `transgene_locus_evidence.csv`. With `--export-clean`, `canonical.nontransgene.vcf` contains chromosome-1-12 data records and retained headers.
