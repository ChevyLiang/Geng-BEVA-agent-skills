# IGCI definitions

## Diagnostic panel

Preferred production columns: `chr,POS,Geng,Xian`.

- `POS` is 1-based.
- `Geng` is the confirmed japonica (geng) diagnostic allele.
- `Xian` is the confirmed indica (xian) diagnostic allele.

For backward compatibility only, `chr,POS,Ref,Alt` is accepted, with panel `Ref` interpreted as japonica and panel `Alt` as indica. These legacy column names describe the panel definition and must never be conflated with the REF/ALT fields of the input VCF. New panels should use `Geng/Xian`.

When the VCF REF allele equals the panel japonica allele and VCF ALT equals the panel indica allele: `0/0 -> 0`, heterozygote `-> 0.5`, `1/1 -> 1`.
When the VCF REF/ALT orientation is reversed relative to the panel: `0/0 -> 1`, heterozygote `-> 0.5`, `1/1 -> 0`.
Phasing separators do not change dosage. Missing or partial-missing calls and incompatible allele pairs are excluded.

## Output statistics

- `GENG_ALLELE` / `XIAN_ALLELE`: panel-defined biological ancestry alleles.
- `VCF_REF` / `VCF_ALT` / `VCF_GT`: local VCF encoding and raw genotype.
- `ALIGNMENT`: `exact` or `swapped` after comparing VCF alleles with the panel.
- `ANCESTRY_GENOTYPE`: `GENG_HOM`, `GENG_XIAN_HET`, or `XIAN_HOM`.
- `XIAN_ALLELE_DOSAGE`: corrected indica allele dosage (0, 1, or 2).
- `SCORE`: `XIAN_ALLELE_DOSAGE / 2` (0, 0.5, or 1).
- `SUPPORT_SUM`: sum of informative site scores.
- `INFORMATIVE_SITES`: number of eligible, called, allele-aligned sites for the sample.
- `IGCI`: `SUPPORT_SUM / INFORMATIVE_SITES`.

Genome-wide IGCI pools all informative sites rather than taking an unweighted mean of chromosome values.

## Sliding windows

Window coordinates use `START` inclusive and `END` exclusive. `COUNT` is the number of informative sites in the window. Empty windows have `COUNT = 0`, missing `SCORE`, and a gray display in the heatmap.

`assets/IJ_panel.csv` is the bundled production diagnostic panel and uses the preferred `Geng/Xian` schema. It is the default panel unless the user explicitly supplies another panel.
