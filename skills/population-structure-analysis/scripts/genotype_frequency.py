from __future__ import annotations
from core.runtime import contract

from typing import Any, Dict

import pandas as pd

from core.manifest import build_manifest
from core.vcf_utils import vcf_to_variant_table
from tools.plink_wrapper import ensure_out_dir


@contract
def genotype_frequency(input_path, out_dir="outputs", **kwargs) -> Dict[str, Any]:
    out_path = ensure_out_dir(out_dir)
    table = vcf_to_variant_table(input_path)
    sample_cols = [c for c in table.columns if c not in {"CHROM", "POS", "LocusID", "REF", "ALT"}]
    long_rows = []
    for _, row in table.iterrows():
        for sample in sample_cols:
            long_rows.append(
                {
                    "LocusID": row["LocusID"],
                    "SampleID": sample,
                    "Genotype": row[sample],
                }
            )
    long_df = pd.DataFrame(long_rows)
    freq = (
        long_df.groupby(["LocusID", "Genotype"], dropna=False)
        .size()
        .reset_index(name="Count")
    )
    result_csv = out_path / "genotype_frequency.csv"
    freq.to_csv(result_csv, index=False, encoding="utf-8-sig")
    return build_manifest("genotype_frequency", input_path, out_path, result_csv=result_csv, extra={"variant_count": int(table.shape[0])})


def run(input_path, out_dir="outputs", **kwargs):
    return genotype_frequency(input_path,out_dir=out_dir,**kwargs)


if __name__ == "__main__":
    raise SystemExit("Import and call genotype_frequency().")

