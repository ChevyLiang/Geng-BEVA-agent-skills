from __future__ import annotations
from _runtime import contract

from typing import Any, Dict

from functional_workflow import (
    annotate_long_table,
    ensure_out_dir,
    filter_annotated_table,
    maybe_prepare_functional_table,
)


@contract
def trait_category_func_extract(input_path, out_dir="outputs", **kwargs) -> Dict[str, Any]:
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
    filtered = filter_annotated_table(
        annotated,
        sample_ids=kwargs.get("sample_ids") or kwargs.get("sample_id"),
        traits=kwargs.get("traits") or kwargs.get("trait"),
        categories=kwargs.get("categories") or kwargs.get("category"),
        query_logic=kwargs.get("query_logic", "and"),
    )

    result_csv = out_path / "trait_category_func_extract.csv"
    filtered.to_csv(result_csv, index=False, encoding="utf-8-sig")

    summary_csv = out_path / "trait_category_func_extract_summary.csv"
    summary = (
        filtered.groupby(["SampleID", "Trait", "Category", "Alternative_Allele_Function", "Alternative_Allele_Direction", "DirectionSource", "AlignmentStatus", "FavorableAlleleGroup", "EffectDirection"], dropna=False)
        .size()
        .reset_index(name="Count")
        .sort_values(["SampleID", "Trait", "Category", "EffectDirection"])
    )
    summary.to_csv(summary_csv, index=False, encoding="utf-8-sig")

    return {
        "module": "trait_category_func_extract",
        "input_path": str(input_path),
        "out_dir": str(out_path),
        "status": "success",
        "outputs": {
            "result_csv": str(result_csv),
            "summary_csv": str(summary_csv),
            "heatmap_png": None,
        },
        "prep": prep_info,
        "filters": {
            "sample_ids": kwargs.get("sample_ids") or kwargs.get("sample_id"),
            "traits": kwargs.get("traits") or kwargs.get("trait"),
            "categories": kwargs.get("categories") or kwargs.get("category"),
            "query_logic": kwargs.get("query_logic", "and"),
        },
    }


def run(input_path, out_dir="outputs", **kwargs):
    return trait_category_func_extract(input_path, out_dir=out_dir, **kwargs)


if __name__ == "__main__":
    raise SystemExit("Import and call trait_category_func_extract().")

