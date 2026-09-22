# IGCI production reference resource

`IJ_panel.csv` is the bundled production japonica/indica diagnostic allele panel used by the IGCI workflow.

Preferred schema:

```text
chr,POS,Geng,Xian
```

- `Geng`: japonica (geng) diagnostic allele.
- `Xian`: indica (xian) diagnostic allele.

These panel alleles define biological ancestry. Input VCF `REF`/`ALT` values only define VCF encoding and are reoriented to `Geng`/`Xian` before indica dosage and IGCI are calculated.

If `--panel` is omitted, the workflow uses this bundled file automatically. A user-supplied panel may be provided explicitly when required.
