from __future__ import annotations
from core.runtime import contract

from pathlib import Path
import sys

_CODEX_ROOT = next(p for p in Path(__file__).resolve().parents if p.name == ".codex")
if str(_CODEX_ROOT / "tools") not in sys.path:
    sys.path.insert(0, str(_CODEX_ROOT / "tools"))
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from core.manifest import build_manifest
from tools.plink_wrapper import ensure_out_dir, run_plink_pca


def _read_eigenvec(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", header=None)
    if df.shape[1] < 3:
        raise ValueError(f"Unexpected PCA eigenvec format: {path}")
    cols = ["FID", "IID"] + [f"PC{i}" for i in range(1, df.shape[1] - 1)]
    df.columns = cols
    return df


def _read_group_file(group_file: str | Path) -> pd.DataFrame:
    path = Path(group_file)
    if not path.exists():
        raise FileNotFoundError(f"Group file not found: {path}")
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    df = pd.read_csv(path, sep=sep)
    if "SampleID" not in df.columns:
        cols_lower = {c.lower(): c for c in df.columns}
        if "iid" in cols_lower:
            df = df.rename(columns={cols_lower["iid"]: "SampleID"})
        elif df.shape[1] >= 1:
            df = df.rename(columns={df.columns[0]: "SampleID"})
    if "Group" not in df.columns:
        cols_lower = {c.lower(): c for c in df.columns}
        if "population" in cols_lower:
            df = df.rename(columns={cols_lower["population"]: "Group"})
        elif df.shape[1] >= 2:
            df = df.rename(columns={df.columns[1]: "Group"})
        else:
            df["Group"] = "Group1"
    return df[["SampleID", "Group"]].copy()


def _pc_columns(df: pd.DataFrame) -> List[str]:
    return [c for c in df.columns if c.startswith("PC")]


def _distance(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.sum((a - b) ** 2)))


def _nearest_farthest_summary(eigenvec: pd.DataFrame, sample_id: str) -> pd.DataFrame:
    matches = eigenvec.loc[eigenvec["IID"].astype(str) == str(sample_id)]
    if matches.empty:
        raise ValueError(f"Sample not found in PCA results: {sample_id}")
    pc_cols = _pc_columns(eigenvec)
    if not pc_cols:
        raise ValueError("No PCA coordinates available")
    base = matches.iloc[0][pc_cols].astype(float).to_numpy()
    rows: List[Tuple[str, float]] = []
    for _, row in eigenvec.iterrows():
        other = str(row["IID"])
        if other == str(sample_id):
            continue
        vec = row[pc_cols].astype(float).to_numpy()
        rows.append((other, _distance(base, vec)))
    if not rows:
        raise ValueError("Need at least two samples for nearest/farthest summary")
    rows_sorted = sorted(rows, key=lambda x: x[1])
    nearest = rows_sorted[0]
    farthest = rows_sorted[-1]
    return pd.DataFrame(
        [
            {
                "SampleID": sample_id,
                "NearestSample": nearest[0],
                "NearestDistance": nearest[1],
                "FarthestSample": farthest[0],
                "FarthestDistance": farthest[1],
                "DistanceMetric": "euclidean_on_all_pcs",
            }
        ]
    )


def _prepare_plot_df(eigenvec: pd.DataFrame, group_file: str | Path | None) -> pd.DataFrame:
    plot_df = eigenvec.copy()
    if group_file:
        groups = _read_group_file(group_file)
        plot_df = plot_df.merge(groups, how="left", left_on="IID", right_on="SampleID")
        plot_df["Group"] = plot_df["Group"].fillna("Ungrouped")
    else:
        plot_df["Group"] = "All"
    return plot_df


def _write_interactive_pca(plot_df: pd.DataFrame, out_html: Path) -> Path:
    from dependency_resolver import ensure_python_package
    ensure_python_package("plotly", start=__file__)
    import plotly.express as px

    hover_data = {"IID": True, "PC1": ":.5f", "PC2": ":.5f", "Group": True}
    fig = px.scatter(
        plot_df,
        x="PC1",
        y="PC2",
        color="Group",
        hover_name="IID",
        hover_data=hover_data,
        title="PCA",
    )
    fig.update_traces(marker={"size": 9, "opacity": 0.85})
    fig.update_layout(
        xaxis_title="PC1",
        yaxis_title="PC2",
        legend_title_text="Group",
        template="plotly_white",
    )
    fig.write_html(str(out_html), include_plotlyjs=True, full_html=True)
    return out_html


def _write_static_pca(plot_df: pd.DataFrame, out_png: Path) -> Optional[Path]:
    try:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(9, 7))
        for group, sub in plot_df.groupby("Group", dropna=False):
            ax.scatter(sub["PC1"], sub["PC2"], s=55, alpha=0.8, label=str(group))
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        ax.set_title("PCA")
        if plot_df["Group"].nunique(dropna=False) > 1:
            ax.legend(title="Group", frameon=False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close(fig)
        return out_png
    except Exception:
        return None


@contract
def pca_analysis(input_path, out_dir="outputs", **kwargs) -> Dict[str, Any]:
    out_path = ensure_out_dir(out_dir)
    plink_cmd = kwargs.get("plink_cmd")
    group_file = kwargs.get("group_file") or kwargs.get("groups")
    focus_sample = kwargs.get("sample_id") or kwargs.get("sample") or kwargs.get("target_sample")

    pca_prefix = out_path / "pca"
    plink_outputs = run_plink_pca(input_path, pca_prefix, plink_cmd=plink_cmd)
    eigenvec = _read_eigenvec(Path(plink_outputs["eigenvec"]))

    result_csv = out_path / "pca_results.csv"
    eigenvec.to_csv(result_csv, index=False, encoding="utf-8-sig")

    plot_df = _prepare_plot_df(eigenvec, group_file)
    interactive_html = _write_interactive_pca(plot_df, out_path / "pca_interactive.html")
    plot_png = _write_static_pca(plot_df, out_path / "pca_plot.png")

    group_csv = None
    nearest_csv = None
    if group_file:
        groups = _read_group_file(group_file)
        group_csv = out_path / "pca_groups.csv"
        groups.to_csv(group_csv, index=False, encoding="utf-8-sig")
    if focus_sample:
        nearest_csv = out_path / f"pca_{focus_sample}_nearest_farthest.csv"
        _nearest_farthest_summary(eigenvec, str(focus_sample)).to_csv(nearest_csv, index=False, encoding="utf-8-sig")

    return build_manifest(
        "pca_analysis",
        input_path,
        out_path,
        result_csv=result_csv,
        plot_html=interactive_html,
        plot_png=plot_png,
        extra={
            "plink_outputs": plink_outputs,
            "group_file": str(group_file) if group_file else None,
            "group_csv": str(group_csv) if group_csv else None,
            "nearest_farthest_csv": str(nearest_csv) if nearest_csv else None,
            "focus_sample": str(focus_sample) if focus_sample else None,
            "interactive_hover": ["SampleID", "PC1", "PC2", "Group"],
        },
    )


def run(input_path, out_dir="outputs", **kwargs):
    return pca_analysis(input_path,out_dir=out_dir,**kwargs)


if __name__ == "__main__":
    raise SystemExit("Import and call pca_analysis().")
