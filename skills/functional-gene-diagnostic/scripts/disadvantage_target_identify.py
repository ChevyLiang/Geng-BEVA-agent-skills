from __future__ import annotations
from _runtime import contract

from typing import Any, Dict

import pandas as pd

from functional_workflow import (
    ARROW_DOWN,
    annotate_long_table,
    ensure_out_dir,
    maybe_prepare_functional_table,
)


@contract
def disadvantage_target_identify(input_path, out_dir="outputs", **kwargs) -> Dict[str, Any]:
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

    disadvantage = annotated.loc[annotated["CarriesFavorableAllele"] == "no"].copy()

    ranking = (
        disadvantage.groupby(["Gene", "LocusID", "Trait", "Category", "Alternative_Allele_Function", "Alternative_Allele_Direction", "DirectionSource", "AlignmentStatus", "FavorableAlleleGroup"], dropna=False)
        .agg(
            DisadvantageCount=("EffectDirection", "size"),
            SampleCount=("SampleID", "nunique"),
            MissingPreferredAlleleCount=("CarriesFavorableAllele", lambda s: int((s == "no").sum())),
        )
        .reset_index()
        .sort_values(
            ["DisadvantageCount", "MissingPreferredAlleleCount", "SampleCount", "Gene", "LocusID"],
            ascending=[False, False, False, True, True],
        )
    )

    result_csv = out_path / "disadvantage_target_identify.csv"
    ranking.to_csv(result_csv, index=False, encoding="utf-8-sig")

    summary = (
        annotated.groupby(["Trait", "Category", "Alternative_Allele_Function", "Alternative_Allele_Direction", "DirectionSource", "AlignmentStatus", "FavorableAlleleGroup", "EffectDirection"], dropna=False)
        .size()
        .reset_index(name="Count")
        .sort_values(["Trait", "Category", "EffectDirection"])
    )
    summary_csv = out_path / "disadvantage_target_identify_summary.csv"
    summary.to_csv(summary_csv, index=False, encoding="utf-8-sig")

    return {
        "module": "disadvantage_target_identify",
        "input_path": str(input_path),
        "out_dir": str(out_path),
        "status": "success",
        "outputs": {
            "result_csv": str(result_csv),
            "summary_csv": str(summary_csv),
            "heatmap_png": None,
        },
        "prep": prep_info,
        "warnings": [],
    }


def run(input_path, out_dir="outputs", **kwargs):
    return disadvantage_target_identify(input_path, out_dir=out_dir, **kwargs)


if __name__ == "__main__":
    raise SystemExit("Import and call disadvantage_target_identify().")

