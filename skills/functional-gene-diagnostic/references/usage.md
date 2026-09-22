
## Windows Codex runtime

Prefer `.codex/run-geng-beva.ps1` as the execution entry point. It verifies or repairs the fixed Geng-BEVA Python 3.11 runtime before dispatching the skill. Do not consider a requested output complete when a required Python package is missing; repair the runtime and rerun. PLINK remains separately managed and is never auto-installed.

# Usage

From the repository root:

```bash
python skills/functional-gene-diagnostic/main.py samples data/input.vcf --out-dir outputs/sample_ids
python skills/functional-gene-diagnostic/main.py single data/input.vcf --sample-id SAMPLE_A --out-dir outputs/functional
python skills/functional-gene-diagnostic/main.py population data/input.vcf --out-dir outputs/pop_functional
python skills/functional-gene-diagnostic/main.py trait data/input.vcf --trait "heading date" --out-dir outputs/heading
python skills/functional-gene-diagnostic/main.py targets data/input.vcf --out-dir outputs/targets
```

The default functional annotation is `assets/RiceNavi_subset.csv`. It uses explicit RiceNavi biological allele names (`Nipponbare_Allele`, `Alternative_Allele`) and the workflow corrects VCF REF/ALT orientation before genotype interpretation.

Without any preference file, the workflow uses binary homozygous rules wherever RiceNavi is resolvable: `Positive` => Alternative favorable / Nipponbare unfavorable; `Negative` => Nipponbare favorable / Alternative unfavorable; heading `promote` => Alternative favorable; heading `delay` => Nipponbare favorable. Heterozygotes are always reported separately. Use `--breeding-preferences` only to override these defaults.

Single-sample mode writes:

- `favorable_genotypes.csv`
- optional heatmap

All analysis modes except `samples` can also accept a precomputed functional long-table CSV. Use `--functional-qtn-csv` to override the bundled RiceNavi-derived table. Legacy annotation CSV column names are accepted for backward compatibility but are normalized to the new schema.

VCF functional-locus extraction is performed directly in Python; PLINK is not required. The legacy `--plink-cmd` option remains accepted for compatibility.

On Windows, absolute paths are recommended. If a Codex sandbox reports `CreateProcessWithLogonW failed: 267`, invoke the existing skill entry point from an ASCII-only working directory (for example `%TEMP%\geng-beva-run`) while passing the VCF and output directory as absolute paths.

Single-sample primary outputs are `favorable_genotypes.csv` (favorable homozygotes plus heterozygotes) and `unfavorable_genotypes.csv` (unfavorable homozygotes). Both use only: `GeneName`, `Position`, `NipponbareGenotype`, `AlternativeGenotype`, `SampleGenotype`, `GenotypeBackground`, `PhenotypicFunction`, and `GenotypeEffect`.
