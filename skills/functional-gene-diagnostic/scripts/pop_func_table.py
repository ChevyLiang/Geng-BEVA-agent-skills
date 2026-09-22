from __future__ import annotations
from _runtime import contract

from pathlib import Path
from typing import Any, Dict

import pandas as pd

from functional_workflow import (
    annotate_long_table,
    ensure_out_dir,
    maybe_prepare_functional_table,
)


@contract
def pop_func_table(input_path, out_dir="outputs", **kwargs) -> Dict[str, Any]:
    out_path = ensure_out_dir(out_dir)
    qtn_path = kwargs.get("functional_qtn_csv") or kwargs.get("qtn_csv")
    plink_cmd = kwargs.get("plink_cmd")

    long_df, prep_info = maybe_prepare_functional_table(
        input_path,
        out_path,
        functional_qtn_csv=qtn_path,
        plink_cmd=plink_cmd,
    )
    qtn_csv = prep_info["functional_qtn_csv"]
    annotated = annotate_long_table(long_df, qtn_csv, kwargs.get("breeding_preferences"), kwargs.get("environment"))

    result = (
        annotated.groupby(
            [
                "Gene",
                "LocusID",
                "Trait",
                "Category",
                "Alternative_Allele_Function",
                "Alternative_Allele_Direction",
                "DirectionSource",
                "AlignmentStatus",
                "FavorableAlleleGroup",
                "EffectDirection",
                "RiceNaviGenotypeState",
                "CarriesFavorableAllele",
            ],
            dropna=False,
        )
        .size()
        .reset_index(name="Count")
        .sort_values(["Trait", "Category", "Gene", "LocusID", "EffectDirection"])
    )

    result_csv = out_path / "pop_func_table.csv"
    result.to_csv(result_csv, index=False, encoding="utf-8-sig")

    summary = (
        annotated.groupby(["SampleID", "DirectionSource", "EffectDirection", "AlignmentStatus"], dropna=False)
        .size()
        .reset_index(name="Count")
        .sort_values(["SampleID", "EffectDirection"])
    )
    summary_csv = out_path / "pop_func_table_summary.csv"
    summary.to_csv(summary_csv, index=False, encoding="utf-8-sig")

    return {
        "module": "pop_func_table",
        "input_path": str(input_path),
        "out_dir": str(out_path),
        "status": "success",
        "outputs": {
            "result_csv": str(result_csv),
            "summary_csv": str(summary_csv),
            "heatmap_png": None,
        },
        "prep": prep_info,
    }


def run(input_path, out_dir="outputs", **kwargs):
    return pop_func_table(input_path, out_dir=out_dir, **kwargs)


if __name__ == "__main__":
    raise SystemExit("Import and call pop_func_table().")

