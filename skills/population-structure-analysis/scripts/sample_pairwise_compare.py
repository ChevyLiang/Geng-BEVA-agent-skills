from __future__ import annotations
from core.runtime import contract

from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd

from core.manifest import build_manifest
from core.vcf_utils import is_missing_genotype, is_standard_chromosome, normalize_chrom
from tools.plink_wrapper import ensure_out_dir, run_plink_ibs, run_plink_pca


def _read_eigenvec(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", header=None)
    cols = ["FID", "IID"] + [f"PC{i}" for i in range(1, df.shape[1] - 1)]
    df.columns = cols
    return df


def _pairwise_distance_from_pca(eigenvec: pd.DataFrame, sample_a: str, sample_b: str) -> pd.DataFrame:
    a = eigenvec.loc[eigenvec["IID"].astype(str) == str(sample_a)]
    b = eigenvec.loc[eigenvec["IID"].astype(str) == str(sample_b)]
    if a.empty or b.empty:
        raise ValueError(f"Samples not found in PCA results: {sample_a}, {sample_b}")
    pc_cols = [c for c in eigenvec.columns if c.startswith("PC")]
    va = a.iloc[0][pc_cols].astype(float).to_numpy()
    vb = b.iloc[0][pc_cols].astype(float).to_numpy()
    diff = va - vb
    dist = float(np.sqrt(np.sum(diff ** 2)))
    return pd.DataFrame(
        [
            {
                "SampleA": sample_a,
                "SampleB": sample_b,
                "PC_Distance": dist,
                **{f"Delta_{pc}": float(da - db) for pc, da, db in zip(pc_cols, va, vb)},
            }
        ]
    )


def _read_genome(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+")
    if df.empty:
        raise ValueError(f"Empty PLINK genome file: {path}")
    return df


def _pairwise_metric(genome: pd.DataFrame, sample_a: str, sample_b: str, value_col: str) -> float:
    mask = (
        (genome["IID1"].astype(str) == str(sample_a)) & (genome["IID2"].astype(str) == str(sample_b))
    ) | (
        (genome["IID1"].astype(str) == str(sample_b)) & (genome["IID2"].astype(str) == str(sample_a))
    )
    hit = genome.loc[mask]
    if hit.empty:
        raise ValueError(f"Pair not found in genome file: {sample_a}, {sample_b}")
    return float(hit.iloc[0][value_col])


def _pairwise_vcf_differences(vcf_path: str | Path, sample_a: str, sample_b: str) -> pd.DataFrame:
    vcf_path = Path(vcf_path)
    rows: List[Dict[str, Any]] = []
    ia = ib = None
    with vcf_path.open("r", encoding="utf-8-sig") as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                parts = line.rstrip("\n").split("\t")
                samples = parts[9:]
                sample_index = {s: i for i, s in enumerate(samples)}
                if str(sample_a) not in sample_index or str(sample_b) not in sample_index:
                    raise ValueError(f"Samples not found in VCF: {sample_a}, {sample_b}")
                ia = 9 + sample_index[str(sample_a)]
                ib = 9 + sample_index[str(sample_b)]
                continue
            if line.startswith("#") or ia is None or ib is None:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            chrom, pos, vid, ref, alt = parts[:5]
            if not is_standard_chromosome(chrom):
                continue
            locus_id = vid if vid != "." else f"{normalize_chrom(chrom)}_{pos}"
            gt_a = parts[ia].split(":", 1)[0].replace("|", "/")
            gt_b = parts[ib].split(":", 1)[0].replace("|", "/")
            if is_missing_genotype(gt_a) or is_missing_genotype(gt_b):
                continue
            if gt_a == gt_b:
                continue
            rows.append(
                {
                    "CHROM": chrom,
                    "POS": int(pos),
                    "LocusID": locus_id,
                    "REF": ref,
                    "ALT": alt,
                    "SampleA_GT": gt_a,
                    "SampleB_GT": gt_b,
                }
            )
    return pd.DataFrame.from_records(rows)


def _difference_heatmap_data(diff_df: pd.DataFrame) -> pd.DataFrame:
    if diff_df.empty:
        return pd.DataFrame()
    table = diff_df.copy()
    table["ChromOrder"] = pd.to_numeric(table["CHROM"].astype(str).map(normalize_chrom), errors="coerce")
    table = table.loc[table["ChromOrder"].between(1, 12, inclusive="both")].copy()
    table["BinMb"] = (table["POS"].astype(int) // 1_000_000).astype(int)
    heatmap = (
        table.groupby(["ChromOrder", "BinMb"])
        .size()
        .reset_index(name="DiffCount")
        .pivot(index="ChromOrder", columns="BinMb", values="DiffCount")
        .fillna(0)
    )
    heatmap = heatmap.reindex(index=list(range(1, 13)), fill_value=0)
    heatmap.index.name = "Chromosome"
    return heatmap


def _plot_blocked_difference_heatmap(diff_df: pd.DataFrame, heatmap_png: Path, ibs_value: float) -> None:
    import matplotlib.pyplot as plt
    import seaborn as sns

    heatmap = _difference_heatmap_data(diff_df)
    if heatmap.empty:
        raise ValueError("No difference data available for heatmap")

    fig_width = max(10, heatmap.shape[1] * 0.5)
    fig_height = max(5, heatmap.shape[0] * 0.42)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    sns.heatmap(
        heatmap,
        cmap="magma",
        annot=False,
        linewidths=0.4,
        linecolor="#f0f0f0",
        cbar_kws={"label": "Difference count"},
        ax=ax,
    )
    ax.set_title(f"Pairwise Difference Density | chr1-12 only, missing excluded | IBS={ibs_value:.4f}")
    ax.set_xlabel("Physical position bin (Mb)")
    ax.set_ylabel("Chromosome")
    ax.set_yticklabels([str(i) for i in range(1, 13)], rotation=0)
    max_col = int(heatmap.shape[1])
    if max_col > 0:
        ax.set_xticks(np.arange(0.5, max_col + 0.5, 1.0))
        ax.set_xticklabels([str(i) for i in heatmap.columns], rotation=45, ha="right")
    for y in range(1, 12):
        ax.axhline(y, color="white", lw=1.2, alpha=0.9)
    plt.tight_layout()
    fig.savefig(heatmap_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _heatmap_summary(diff_df: pd.DataFrame, sample_a: str, sample_b: str, ibs_value: float) -> pd.DataFrame:
    if diff_df.empty:
        return pd.DataFrame(
            [
                {
                    "SampleA": sample_a,
                    "SampleB": sample_b,
                    "DiffSites": 0,
                    "IBS": ibs_value,
                    "ChromosomesWithDiff": 0,
                    "TopChromosome": "",
                    "TopChromosomeDiffCount": 0,
                }
            ]
        )
    counts = diff_df["CHROM"].astype(str).value_counts()
    return pd.DataFrame(
        [
            {
                "SampleA": sample_a,
                "SampleB": sample_b,
                "DiffSites": int(diff_df.shape[0]),
                "IBS": ibs_value,
                "ChromosomesWithDiff": int(diff_df["CHROM"].astype(str).nunique()),
                "TopChromosome": str(counts.idxmax()),
                "TopChromosomeDiffCount": int(counts.max()),
            }
        ]
    )


@contract
def sample_pairwise_compare(input_path, out_dir="outputs", **kwargs) -> Dict[str, Any]:
    out_path = ensure_out_dir(out_dir)
    sample_a = kwargs.get("sample_a") or kwargs.get("sample1")
    sample_b = kwargs.get("sample_b") or kwargs.get("sample2")
    if not sample_a or not sample_b:
        raise ValueError("sample_a and sample_b are required")
    plink_cmd = kwargs.get("plink_cmd")

    pca_prefix = out_path / f"pairwise_{sample_a}_{sample_b}_pca"
    pca_outputs = run_plink_pca(
        input_path,
        pca_prefix,
        plink_cmd=plink_cmd,
    )
    eigenvec = _read_eigenvec(Path(pca_outputs["eigenvec"]))
    pca_distance = _pairwise_distance_from_pca(eigenvec, sample_a, sample_b)

    genome_prefix = out_path / f"pairwise_{sample_a}_{sample_b}_ibs"
    genome_outputs = run_plink_ibs(
        input_path,
        genome_prefix,
        plink_cmd=plink_cmd,
    )
    genome = _read_genome(Path(genome_outputs["genome"]))
    ibs_value_col = "DST"
    ibs_value = _pairwise_metric(genome, sample_a, sample_b, ibs_value_col)

    result = pca_distance.copy()
    result["IBS"] = ibs_value
    ibd_value = _pairwise_metric(genome, sample_a, sample_b, "PI_HAT")
    result["IBD"] = ibd_value

    diff_df = _pairwise_vcf_differences(input_path, sample_a, sample_b)
    diff_csv = out_path / f"pairwise_{sample_a}_{sample_b}_variant_differences.csv"
    diff_df.to_csv(diff_csv, index=False, encoding="utf-8-sig")

    result_csv = out_path / f"pairwise_{sample_a}_{sample_b}_compare.csv"
    result.to_csv(result_csv, index=False, encoding="utf-8-sig")

    heatmap_png = None
    planned_heatmap = out_path / f"pairwise_{sample_a}_{sample_b}_heatmap.png"
    heatmap_summary_csv = None
    try:
        diff_heatmap = _difference_heatmap_data(diff_df)
        if not diff_heatmap.empty:
            _plot_blocked_difference_heatmap(diff_df, planned_heatmap, ibs_value)
            heatmap_png = planned_heatmap
            heatmap_summary_csv = out_path / f"pairwise_{sample_a}_{sample_b}_heatmap_summary.csv"
            _heatmap_summary(diff_df, sample_a, sample_b, ibs_value).to_csv(
                heatmap_summary_csv, index=False, encoding="utf-8-sig"
            )
    except Exception:
        heatmap_png = None

    return build_manifest(
        "sample_pairwise_compare",
        input_path,
        out_path,
        result_csv=result_csv,
        heatmap_png=heatmap_png,
        extra={
            "sample_a": sample_a,
            "sample_b": sample_b,
            "variant_differences_csv": str(diff_csv),
            "pca_outputs": pca_outputs,
            "genome_outputs": genome_outputs,
            "ibs_similarity": ibs_value,
            "ibd_similarity": float(ibd_value) if np.isfinite(ibd_value) else None,
            "heatmap_summary_csv": str(heatmap_summary_csv) if heatmap_summary_csv else None,
        },
    )


def run(input_path, out_dir="outputs", **kwargs):
    return sample_pairwise_compare(input_path,out_dir=out_dir,**kwargs)


if __name__ == "__main__":
    raise SystemExit("Import and call sample_pairwise_compare().")

