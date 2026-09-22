from __future__ import annotations
from core.runtime import contract

from pathlib import Path
import sys

_CODEX_ROOT = next(p for p in Path(__file__).resolve().parents if p.name == ".codex")
if str(_CODEX_ROOT / "tools") not in sys.path:
    sys.path.insert(0, str(_CODEX_ROOT / "tools"))
from typing import Any, Dict

import pandas as pd

from core.manifest import build_manifest
from tools.plink_wrapper import ensure_out_dir, run_plink_ibs, run_plink_ibd


def _read_genome(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+")
    if df.empty:
        raise ValueError(f"Empty PLINK genome file: {path}")
    return df


def _matrix_from_pairs(df: pd.DataFrame, value_col: str, id1_col: str = "IID1", id2_col: str = "IID2") -> pd.DataFrame:
    ids = sorted(set(df[id1_col].astype(str)) | set(df[id2_col].astype(str)))
    matrix = pd.DataFrame(index=ids, columns=ids, dtype=float)
    for _, row in df.iterrows():
        i = str(row[id1_col])
        j = str(row[id2_col])
        matrix.loc[i, j] = row[value_col]
        matrix.loc[j, i] = row[value_col]
    for x in ids:
        matrix.loc[x, x] = 1.0
    return matrix


def _write_interactive_heatmap(matrix: pd.DataFrame, out_html: Path, mode: str) -> Path:
    from dependency_resolver import ensure_python_package
    ensure_python_package("plotly", start=__file__)
    import plotly.graph_objects as go

    metric_label = "IBS" if mode == "ibs" else "IBD"
    z = matrix.to_numpy(dtype=float)
    x = [str(v) for v in matrix.columns]
    y = [str(v) for v in matrix.index]
    fig = go.Figure(
        data=go.Heatmap(
            z=z,
            x=x,
            y=y,
            colorbar={"title": metric_label},
            hovertemplate=(
                "Sample A: %{y}<br>"
                "Sample B: %{x}<br>"
                f"{metric_label}: %{{z:.5f}}<extra></extra>"
            ),
        )
    )
    fig.update_layout(
        title=f"Pairwise {metric_label} heatmap",
        xaxis_title="Sample B",
        yaxis_title="Sample A",
        template="plotly_white",
    )
    fig.update_yaxes(autorange="reversed")
    fig.write_html(str(out_html), include_plotlyjs=True, full_html=True)
    return out_html


def _write_static_heatmap(matrix: pd.DataFrame, out_png: Path, mode: str):
    try:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 7))
        image = ax.imshow(matrix.to_numpy(dtype=float), aspect="auto")
        ax.set_xticks(range(len(matrix.columns)))
        ax.set_xticklabels(matrix.columns, rotation=90, fontsize=7)
        ax.set_yticks(range(len(matrix.index)))
        ax.set_yticklabels(matrix.index, fontsize=7)
        label = "IBS" if mode == "ibs" else "IBD"
        ax.set_title(f"Pairwise {label} heatmap")
        fig.colorbar(image, ax=ax, label=label)
        fig.tight_layout()
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close(fig)
        return out_png
    except Exception:
        return None


@contract
def ibs_analysis(input_path, out_dir="outputs", **kwargs) -> Dict[str, Any]:
    out_path = ensure_out_dir(out_dir)
    plink_cmd = kwargs.get("plink_cmd")
    mode = str(kwargs.get("mode", "ibs")).lower()

    prefix = out_path / f"{mode}"
    if mode == "ibd":
        plink_outputs = run_plink_ibd(input_path, prefix, plink_cmd=plink_cmd)
    else:
        plink_outputs = run_plink_ibs(input_path, prefix, plink_cmd=plink_cmd)

    genome = _read_genome(Path(plink_outputs["genome"]))
    value_col = "PI_HAT" if mode == "ibd" else "DST"
    if value_col not in genome.columns:
        raise ValueError(f"PLINK genome output missing {value_col}")
    matrix = _matrix_from_pairs(genome, value_col)

    result_csv = out_path / f"{mode}_matrix.csv"
    matrix.to_csv(result_csv, encoding="utf-8-sig")

    heatmap_html = _write_interactive_heatmap(matrix, out_path / f"{mode}_heatmap_interactive.html", mode)
    heatmap_png = _write_static_heatmap(matrix, out_path / f"{mode}_heatmap.png", mode)

    return build_manifest(
        "ibs_analysis",
        input_path,
        out_path,
        result_csv=result_csv,
        plot_html=heatmap_html,
        plot_png=heatmap_png,
        extra={
            "mode": mode,
            "plink_outputs": plink_outputs,
            "matrix_value": value_col,
            "interactive_hover": ["Sample A", "Sample B", "IBS" if mode == "ibs" else "IBD"],
        },
    )


def run(input_path, out_dir="outputs", **kwargs):
    return ibs_analysis(input_path,out_dir=out_dir,**kwargs)


if __name__ == "__main__":
    raise SystemExit("Import and call ibs_analysis().")
