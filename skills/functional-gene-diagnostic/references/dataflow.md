# Dataflow and field meanings

## Allele coordinate systems

The functional workflow deliberately keeps two coordinate systems separate:

1. **VCF encoding**: `VCFRefAllele`, `VCFAltAllele`, and the raw `VCFGenotype`.
2. **RiceNavi biological encoding**: `NipponbareAllele`, `AlternativeAllele`, and `RiceNaviGenotypeState`.

VCF REF/ALT is not assumed to equal Nipponbare/Alternative. The workflow first aligns literal allele sequences. `AlignmentStatus=exact` means VCF REF=Nipponbare and VCF ALT=Alternative; `swapped` means VCF REF=Alternative and VCF ALT=Nipponbare. Unresolved alignment remains ambiguous and is not used for favorable-genotype calls.

`RiceNaviGenotypeState` values:

- `NIPPONBARE_HOM`
- `ALTERNATIVE_HOM`
- `HET`
- `MISSING`
- `AMBIGUOUS`

`AlternativeAlleleDosage` is 0, 1, or 2 after allele correction. It is not the raw VCF ALT dosage when alignment is swapped.

## Bundled RiceNavi-derived CSV

Preferred columns:

- `CHROM`, `POS`, `ID`
- `Alternative_Allele_Effect`
- `Nipponbare_Allele`
- `Alternative_Allele`
- `Alternative_Allele_Function`
- `GENEname`, `Category`, `Trait`
- `Alternative_Allele_Direction`

Legacy custom tables containing `Effect,Unfunctional,Functional,ALT_FUNCTION,alt_function` are accepted for backward compatibility and normalized internally, but new data should use the preferred names above.

## Favorable-allele defaults

`FavorableRule` records the decision source:

- `RiceNavi_positive`: Alternative favorable, Nipponbare unfavorable.
- `RiceNavi_negative`: Nipponbare favorable, Alternative unfavorable.
- `early_heading_alternative`: Alternative promotes/advances heading, so Alternative favorable and Nipponbare unfavorable.
- `early_heading_nipponbare`: Alternative delays heading, so Nipponbare favorable and Alternative unfavorable.
- `explicit_preference`: user preference overrides project defaults.
- `unresolved`: no default or explicit favorable allele can be assigned safely.

`FavorableAlleleGroup` is `NIPPONBARE`, `ALTERNATIVE`, or `UNKNOWN` in RiceNavi coordinates.
`FavorableVCFAlleleGroup` is `REF`, `ALT`, or `UNKNOWN` in the current VCF encoding.
`CarriesFavorableAllele` is yes/no, blank when unresolved or missing.
`FavorableGenotypeStatus` distinguishes `homozygous_favorable`, `homozygous_unfavorable`, `heterozygous`, and `unresolved`.

`AdvantageClass` provides the requested five biological genotype classes:

- `NIPPONBARE_FAVORABLE`
- `NIPPONBARE_UNFAVORABLE`
- `ALTERNATIVE_FAVORABLE`
- `ALTERNATIVE_UNFAVORABLE`
- `HETEROZYGOUS`

`UNRESOLVED` is not a sixth biological class; it is a QC status used only when the site is missing/unalignable or when RiceNavi effect is genuinely unresolved (for example a non-heading `Uncertainty` record without an explicit preference).

## Phenotype direction

`EffectDirection` describes the phenotype associated with the RiceNavi Alternative allele. Values may include increase, decrease, gain, loss, delay, promote, prevent, neutral, or unknown. It is separate from favorability.

The default earlier-heading rule interprets promoting/advancing heading as favorable and delaying heading as unfavorable for the Alternative allele. Conditional descriptions such as `under LD` are preserved in the raw annotation.

## Outputs

Single-sample mode writes the complete aligned table plus:

- `favorable_genotypes.csv`

Heterozygotes remain a separate class; no dominance is assumed.

Single-sample primary outputs are `favorable_genotypes.csv` (favorable homozygotes plus heterozygotes) and `unfavorable_genotypes.csv` (unfavorable homozygotes). Both use only: `GeneName`, `Position`, `NipponbareGenotype`, `AlternativeGenotype`, `SampleGenotype`, `GenotypeBackground`, `PhenotypicFunction`, and `GenotypeEffect`.
