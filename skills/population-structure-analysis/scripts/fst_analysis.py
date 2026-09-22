from __future__ import annotations

from typing import Any, Dict, Sequence

import numpy as np
import pandas as pd

from core.manifest import build_manifest
from core.runtime import contract
from core.vcf_utils import vcf_to_variant_table
from tools.plink_wrapper import ensure_out_dir

RICE_CHROMOSOMES = [str(i) for i in range(1, 13)]


def _normalize_chrom(chrom) -> str:
    text = str(chrom).strip()
    if text.lower().startswith("chr"):
        text = text[3:]
    try:
        n = int(float(text))
        return str(n)
    except ValueError:
        return text


def _group_indices(samples: Sequence[str], group_a: Sequence[str], group_b: Sequence[str]):
    sample_set = set(samples)
    set_a = [s for s in group_a if s in sample_set]
    set_b = [s for s in group_b if s in sample_set]
    overlap = set(set_a) & set(set_b)
    if overlap:
        raise ValueError(f"Groups must be disjoint; overlapping samples: {sorted(overlap)}")
    if not set_a or not set_b:
        raise ValueError("Both sample groups must contain at least one VCF sample")
    return set_a, set_b


def _encode_gt(gt):
    text = str(gt).split(":", 1)[0].replace("|", "/")
    if text == "0/0":
        return 0.0
    if text in {"0/1", "1/0"}:
        return 1.0
    if text == "1/1":
        return 2.0
    return np.nan


def _weir_cockerham_components(*populations: np.ndarray) -> tuple[float, float, float]:
    """Return Weir-Cockerham a, b, c components for one diploid biallelic site."""
    called = [np.asarray(p, dtype=float)[~np.isnan(np.asarray(p, dtype=float))] for p in populations]
    if len(called) < 2 or any(len(p) == 0 for p in called):
        return np.nan, np.nan, np.nan

    n = np.asarray([len(p) for p in called], dtype=float)
    r = float(len(called))
    n_bar = n.mean()
    if n_bar <= 1.0:
        return np.nan, np.nan, np.nan

    p = np.asarray([p_.mean() / 2.0 for p_ in called], dtype=float)
    h = np.asarray([(p_ == 1.0).mean() for p_ in called], dtype=float)
    p_bar = np.sum(n * p) / (r * n_bar)
    h_bar = np.sum(n * h) / (r * n_bar)

    n_c = (r * n_bar - np.sum(n * n) / (r * n_bar)) / (r - 1.0)
    if n_c <= 0:
        return np.nan, np.nan, np.nan

    s2 = np.sum(n * (p - p_bar) ** 2) / ((r - 1.0) * n_bar)
    a = (n_bar / n_c) * (
        s2
        - (
            p_bar * (1.0 - p_bar)
            - ((r - 1.0) / r) * s2
            - h_bar / 4.0
        ) / (n_bar - 1.0)
    )
    b = (n_bar / (n_bar - 1.0)) * (
        p_bar * (1.0 - p_bar)
        - ((r - 1.0) / r) * s2
        - ((2.0 * n_bar - 1.0) / (4.0 * n_bar)) * h_bar
    )
    c = h_bar / 2.0
    return float(a), float(b), float(c)


def _weir_cockerham_fst(*populations: np.ndarray) -> float:
    """Per-variant Weir-Cockerham theta; negative finite estimates are retained."""
    a, b, c = _weir_cockerham_components(*populations)
    denominator = a + b + c
    if not np.isfinite(denominator) or denominator == 0:
        return np.nan
    theta = a / denominator
    return float(theta) if np.isfinite(theta) else np.nan


def _window_fst(site_df: pd.DataFrame, window: int, step: int, min_sites: int) -> pd.DataFrame:
    if window <= 0 or step <= 0:
        raise ValueError("window and step must both be > 0")
    if min_sites < 1:
        raise ValueError("min_sites must be >= 1")

    rows = []
    work = site_df.copy()
    work["CHROM_NORM"] = work["CHROM"].map(_normalize_chrom)
    for chrom in RICE_CHROMOSOMES:
        sub = work[work["CHROM_NORM"] == chrom].sort_values("POS")
        if sub.empty:
            continue
        max_pos = int(sub["POS"].max())
        for start in range(1, max_pos + 1, step):
            end = start + window - 1
            win = sub[(sub["POS"] >= start) & (sub["POS"] <= end)]
            valid = win[np.isfinite(win["A"]) & np.isfinite(win["B"]) & np.isfinite(win["C"])]
            n_sites = int(valid.shape[0])
            if n_sites >= min_sites:
                denom = float((valid["A"] + valid["B"] + valid["C"]).sum())
                weighted = float(valid["A"].sum() / denom) if np.isfinite(denom) and denom != 0 else np.nan
                mean_fst = float(valid["WeirCockerhamFst"].mean())
            else:
                weighted = np.nan
                mean_fst = np.nan
            rows.append({
                "CHROM": chrom,
                "START": start,
                "END": end,
                "MID": start + (window - 1) / 2.0,
                "N_SITES": n_sites,
                "WEIGHTED_FST": weighted,
                "MEAN_FST": mean_fst,
            })
    return pd.DataFrame(rows)


def _cumulative_positions(df: pd.DataFrame, pos_col: str) -> tuple[pd.DataFrame, list[tuple[str, float]]]:
    work = df.copy()
    work["CHROM_NORM"] = work["CHROM"].map(_normalize_chrom)
    offsets = {}
    ticks = []
    cumulative = 0.0
    for chrom in RICE_CHROMOSOMES:
        c = work[work["CHROM_NORM"] == chrom]
        if c.empty:
            continue
        max_pos = float(c[pos_col].max())
        offsets[chrom] = cumulative
        ticks.append((chrom, cumulative + max_pos / 2.0))
        cumulative += max_pos
    work["CUM_POS"] = [float(pos) + offsets.get(chrom, 0.0) for pos, chrom in zip(work[pos_col], work["CHROM_NORM"])]
    return work, ticks


def _plot_manhattan(df: pd.DataFrame, value_col: str, pos_col: str, out_prefix, title: str) -> tuple[str | None, str | None]:
    valid = df[np.isfinite(pd.to_numeric(df[value_col], errors="coerce"))].copy()
    if valid.empty:
        return None, None
    valid, ticks = _cumulative_positions(valid, pos_col)
    try:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(14, 5))
        for chrom in RICE_CHROMOSOMES:
            sub = valid[valid["CHROM_NORM"] == chrom]
            if sub.empty:
                continue
            ax.scatter(sub["CUM_POS"], sub[value_col], s=10, alpha=0.8, linewidths=0)
        ax.axhline(0, linewidth=0.8, linestyle="--")
        ax.set_xticks([x for _, x in ticks])
        ax.set_xticklabels([c for c, _ in ticks])
        ax.set_xlabel("Chromosome")
        ax.set_ylabel("Weir-Cockerham FST")
        ax.set_title(title)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        png = str(out_prefix) + ".png"
        pdf = str(out_prefix) + ".pdf"
        fig.savefig(png, dpi=300, bbox_inches="tight")
        fig.savefig(pdf, bbox_inches="tight")
        plt.close(fig)
        return png, pdf
    except Exception:
        return None, None


@contract
def fst_analysis(input_path, out_dir="outputs", **kwargs) -> Dict[str, Any]:
    out_path = ensure_out_dir(out_dir)
    group_a = kwargs.get("group_a") or []
    group_b = kwargs.get("group_b") or []
    if not group_a or not group_b:
        raise ValueError("group_a and group_b are required")
    sliding_window = bool(kwargs.get("sliding_window", False))
    raw_window = kwargs.get("window")
    raw_step = kwargs.get("step")
    # Supplying window/step explicitly also counts as an explicit request for sliding-window FST.
    if raw_window is not None or raw_step is not None:
        sliding_window = True
    window = int(raw_window if raw_window is not None else 500_000)
    step = int(raw_step if raw_step is not None else 100_000)
    min_sites = int(kwargs.get("min_sites", 1))

    table = vcf_to_variant_table(input_path)
    samples = [c for c in table.columns if c not in {"CHROM", "POS", "LocusID", "REF", "ALT"}]
    ga, gb = _group_indices(samples, group_a, group_b)

    records = []
    for _, row in table.iterrows():
        a_vals = np.asarray([_encode_gt(row[s]) for s in ga], dtype=float)
        b_vals = np.asarray([_encode_gt(row[s]) for s in gb], dtype=float)
        a, b, c = _weir_cockerham_components(a_vals, b_vals)
        denominator = a + b + c
        fst = a / denominator if np.isfinite(denominator) and denominator != 0 else np.nan
        records.append({
            "LocusID": row["LocusID"], "CHROM": row["CHROM"], "POS": int(row["POS"]),
            "WeirCockerhamFst": float(fst) if np.isfinite(fst) else np.nan,
            "A": a, "B": b, "C": c,
        })

    site_result = pd.DataFrame(records)
    site_csv = out_path / "fst_sites.csv"
    site_result.to_csv(site_csv, index=False, encoding="utf-8-sig")

    site_png, site_pdf = _plot_manhattan(
        site_result, "WeirCockerhamFst", "POS", out_path / "fst_manhattan",
        "Genome-wide per-site Weir-Cockerham FST",
    )

    window_csv = None
    window_png = None
    window_pdf = None
    if sliding_window:
        window_result = _window_fst(site_result, window, step, min_sites)
        window_csv = out_path / "fst_windows.csv"
        window_result.to_csv(window_csv, index=False, encoding="utf-8-sig")
        window_png, window_pdf = _plot_manhattan(
            window_result, "WEIGHTED_FST", "MID", out_path / "fst_window_manhattan",
            f"Genome-wide sliding-window FST ({window:,}-bp window, {step:,}-bp step)",
        )

    return build_manifest(
        "fst_analysis", input_path, out_path,
        result_csv=site_csv,
        window_csv=window_csv,
        plot_png=site_png,
        plot_pdf=site_pdf,
        window_plot_png=window_png,
        window_plot_pdf=window_pdf,
        extra={
            "estimator": "Weir-Cockerham theta (1984)",
            "sliding_window": sliding_window,
            "window_estimator": "sum(a) / sum(a+b+c)" if sliding_window else None,
            "group_a": list(ga), "group_b": list(gb),
            "window_bp": window if sliding_window else None,
            "step_bp": step if sliding_window else None,
            "min_sites": min_sites if sliding_window else None,
        },
    )


def run(input_path, out_dir="outputs", **kwargs):
    return fst_analysis(input_path, out_dir=out_dir, **kwargs)


if __name__ == "__main__":
    raise SystemExit("Import and call fst_analysis().")
