from __future__ import annotations
from core.runtime import contract

from pathlib import Path
from typing import Any, Dict, Sequence

import numpy as np
import pandas as pd

from core.manifest import build_manifest
from tools.plink_wrapper import ensure_out_dir, run_plink_ibs, run_plink_pca


def _read_eigenvec(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", header=None)
    cols = ["FID", "IID"] + [f"PC{i}" for i in range(1, df.shape[1] - 1)]
    df.columns = cols
    return df


def _pca_distances(eigenvec: pd.DataFrame, samples: Sequence[str]) -> pd.DataFrame:
    pc_cols = [c for c in eigenvec.columns if c.startswith("PC")]
    rows = []
    for i, s1 in enumerate(samples):
        for s2 in samples[i + 1 :]:
            a = eigenvec.loc[eigenvec["IID"].astype(str) == str(s1)].iloc[0][pc_cols].astype(float).to_numpy()
            b = eigenvec.loc[eigenvec["IID"].astype(str) == str(s2)].iloc[0][pc_cols].astype(float).to_numpy()
            dist = float(np.sqrt(np.sum((a - b) ** 2)))
            rows.append({"SampleA": s1, "SampleB": s2, "PC_Distance": dist})
    return pd.DataFrame(rows)


def _read_genome(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=r"\s+")


@contract
def pairwise_matrix(input_path, out_dir="outputs", **kwargs) -> Dict[str, Any]:
    out_path = ensure_out_dir(out_dir)
    samples = kwargs.get("samples")
    if not samples or len(samples) < 2:
        raise ValueError("At least two samples are required")
    plink_cmd = kwargs.get("plink_cmd")

    pca_prefix = out_path / "pairwise_matrix_pca"
    pca_outputs = run_plink_pca(input_path, pca_prefix, plink_cmd=plink_cmd)
    eigenvec = _read_eigenvec(Path(pca_outputs["eigenvec"]))
    pca_df = _pca_distances(eigenvec, samples)

    ibs_prefix = out_path / "pairwise_matrix_ibs"
    ibs_outputs = run_plink_ibs(input_path, ibs_prefix, plink_cmd=plink_cmd)
    genome = _read_genome(Path(ibs_outputs["genome"]))
    value_col = "DST"
    if value_col not in genome.columns:
        raise ValueError(f"PLINK genome output missing {value_col}")

    ibs_matrix = pd.DataFrame(index=samples, columns=samples, dtype=float)
    for s in samples:
        ibs_matrix.loc[s, s] = 1.0
    for _, row in genome.iterrows():
        i = str(row["IID1"])
        j = str(row["IID2"])
        if i in samples and j in samples:
            ibs_matrix.loc[i, j] = float(row[value_col])
            ibs_matrix.loc[j, i] = float(row[value_col])

    pca_csv = out_path / "pairwise_matrix_pca_distances.csv"
    ibs_csv = out_path / "pairwise_matrix_ibs.csv"
    pca_df.to_csv(pca_csv, index=False, encoding="utf-8-sig")
    ibs_matrix.to_csv(ibs_csv, encoding="utf-8-sig")

    heatmap_png = out_path / "pairwise_matrix_heatmap.png"
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns

        plt.figure(figsize=(max(6, len(samples) * 0.35), max(5, len(samples) * 0.35)))
        sns.heatmap(ibs_matrix.astype(float), cmap="viridis", annot=False)
        plt.title("Pairwise IBS Matrix")
        plt.tight_layout()
        plt.savefig(heatmap_png, dpi=300, bbox_inches="tight")
        plt.close()
    except Exception:
        heatmap_png = None

    return build_manifest(
        "pairwise_matrix",
        input_path,
        out_path,
        result_csv=pca_csv,
        summary_csv=ibs_csv,
        heatmap_png=heatmap_png,
        extra={"samples": list(samples), "pca_outputs": pca_outputs, "ibs_outputs": ibs_outputs},
    )


def run(input_path, out_dir="outputs", **kwargs):
    return pairwise_matrix(input_path,out_dir=out_dir,**kwargs)


if __name__ == "__main__":
    raise SystemExit("Import and call pairwise_matrix().")

