from __future__ import annotations
import argparse
import json
from pathlib import Path
import pandas as pd
from _runtime import contract, build_manifest
from functional_workflow import annotate_long_table, maybe_prepare_functional_table, filter_annotated_table, ensure_out_dir


def _biological_genotype(row: pd.Series) -> str:
    nip = str(row.get("NipponbareAllele", "") or "")
    alt = str(row.get("AlternativeAllele", "") or "")
    state = str(row.get("RiceNaviGenotypeState", "") or "")
    if state == "NIPPONBARE_HOM":
        return f"{nip}/{nip}"
    if state == "ALTERNATIVE_HOM":
        return f"{alt}/{alt}"
    if state == "HET":
        return f"{nip}/{alt}"
    return ""


def _compact_output(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=[
            "GeneName", "Position", "NipponbareGenotype", "AlternativeGenotype",
            "SampleGenotype", "GenotypeBackground", "PhenotypicFunction", "GenotypeEffect",
        ])
    out = pd.DataFrame(index=df.index)
    gene_col = next((c for c in ("Gene", "GENEname", "GeneName") if c in df.columns), None)
    out["GeneName"] = df[gene_col].fillna("").astype(str).str.strip() if gene_col else ""
    out["Position"] = df.apply(lambda r: f"Chr{str(r.get('CHROM','')).replace('chr','').replace('Chr','')}:{r.get('POS','')}", axis=1)
    out["NipponbareGenotype"] = df["NipponbareAllele"].astype(str).map(lambda a: f"{a}/{a}" if a else "")
    out["AlternativeGenotype"] = df["AlternativeAllele"].astype(str).map(lambda a: f"{a}/{a}" if a else "")
    out["SampleGenotype"] = df.apply(_biological_genotype, axis=1)
    background_map = {
        "NIPPONBARE_HOM": "Nipponbare",
        "ALTERNATIVE_HOM": "Alternative",
        "HET": "Heterozygous",
    }
    out["GenotypeBackground"] = df["RiceNaviGenotypeState"].map(background_map).fillna("Unresolved")
    trait = df["Trait"].fillna("").astype(str).str.strip()
    func = df["Alternative_Allele_Function"].fillna("").astype(str).str.strip()
    out["PhenotypicFunction"] = [
        f"{t}: {f}" if t and f else (f or t) for t, f in zip(trait, func)
    ]
    class_to_effect = {
        "NIPPONBARE_FAVORABLE": "Favorable",
        "ALTERNATIVE_FAVORABLE": "Favorable",
        "NIPPONBARE_UNFAVORABLE": "Unfavorable",
        "ALTERNATIVE_UNFAVORABLE": "Unfavorable",
        "HETEROZYGOUS": "Heterozygous",
    }
    out["GenotypeEffect"] = df["AdvantageClass"].map(class_to_effect).fillna("Unresolved")
    return out.reset_index(drop=True)


@contract
def single_sample_func_table(input_path, out_dir="outputs", **kwargs):
    out_path = ensure_out_dir(out_dir)
    long_df, prep = maybe_prepare_functional_table(
        input_path,
        out_path,
        functional_qtn_csv=kwargs.get("functional_qtn_csv") or kwargs.get("qtn_csv"),
        plink_cmd=kwargs.get("plink_cmd"),
    )
    result = annotate_long_table(
        long_df,
        prep["functional_qtn_csv"],
        kwargs.get("breeding_preferences"),
        kwargs.get("environment"),
    )
    result = filter_annotated_table(result, sample_ids=kwargs.get("sample_id") or kwargs.get("sample_ids"))

    # The long table is an internal extraction artifact for VCF inputs, not a user-facing result.
    internal_long = prep.pop("functional_long_csv", None)
    if prep.get("extraction_backend") == "python_stream" and internal_long:
        try:
            Path(internal_long).unlink(missing_ok=True)
        except OSError:
            pass

    favorable_csv = out_path / "favorable_genotypes.csv"
    unfavorable_csv = out_path / "unfavorable_genotypes.csv"

    favorable_classes = {"NIPPONBARE_FAVORABLE", "ALTERNATIVE_FAVORABLE", "HETEROZYGOUS"}
    unfavorable_classes = {"NIPPONBARE_UNFAVORABLE", "ALTERNATIVE_UNFAVORABLE"}

    favorable = _compact_output(result.loc[result["AdvantageClass"].isin(favorable_classes)].copy())
    unfavorable = _compact_output(result.loc[result["AdvantageClass"].isin(unfavorable_classes)].copy())
    favorable.to_csv(favorable_csv, index=False, encoding="utf-8-sig")
    unfavorable.to_csv(unfavorable_csv, index=False, encoding="utf-8-sig")

    unresolved_count = int((result["AdvantageClass"] == "UNRESOLVED").sum()) if not result.empty else 0
    m = build_manifest(
        "single_sample_func_table",
        input_path,
        out_path,
        extra={"prep": prep},
    )
    m["outputs"].update({
        "favorable_genotypes_csv": str(favorable_csv),
        "unfavorable_genotypes_csv": str(unfavorable_csv),
    })
    m["warnings"] = ([f"{unresolved_count} records were unresolved and were omitted from the two breeding-summary CSV files."] if unresolved_count else [])
    m["metrics"] = {
        "annotated_rows": len(result),
        "favorable_or_heterozygous_rows": len(favorable),
        "unfavorable_rows": len(unfavorable),
        "unresolved_rows": unresolved_count,
        "samples": int(result.SampleID.nunique()) if not result.empty else 0,
    }
    return m


def run(input_path, out_dir="outputs", **kwargs):
    return single_sample_func_table(input_path, out_dir=out_dir, **kwargs)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Functional-locus annotation summarized into favorable and unfavorable genotype files")
    p.add_argument("input_path")
    p.add_argument("--out-dir", default="outputs")
    p.add_argument("--sample-id")
    p.add_argument("--functional-qtn-csv")
    p.add_argument("--plink-cmd")
    p.add_argument("--breeding-preferences")
    p.add_argument("--environment")
    print(json.dumps(run(**vars(p.parse_args())), ensure_ascii=False, indent=2))
