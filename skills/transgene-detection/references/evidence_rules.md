# Evidence and scope

The Geng-BEVA standardized VCF convention used here contains canonical rice chromosome records plus separately encoded transgene-assay records. Noncanonical **data records** are therefore interpreted as assay elements by this workflow. This is a project-specific data contract, not a universal rule about VCF contig names.

Known families include 35S (including CaMV/FMV aliases), BT, CP4_EPSPS, caoganlin, Nos, tNOS, and ubiquitin. Unknown element names are retained rather than discarded.

The evidence unit is one sample at one explicit assay record. `FORMAT` determines the GT and DP fields. Because no independent presence/absence allele map is supplied, evidence tiers depend on a fully called genotype plus sample depth rather than on allele 0 versus allele 1.

Default thresholds are `low >= 10`, `moderate >= 50`, and `high > 100` reads. These are heuristic operational tiers and may be changed with CLI/Python arguments. They are not validated sensitivity/specificity thresholds.

`HasTransgeneSupport` means that at least one called assay record meets the configured low-support depth cutoff. It is not equivalent to a confirmed transgene insertion. No claim is made about insertion position, copy number, inheritance, or phenotypic expression.
