from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

_CODEX_ROOT = next(p for p in Path(__file__).resolve().parents if p.name == ".codex")
sys.path.insert(0, str(_CODEX_ROOT / "tools"))
from runtime_launcher import ensure_fixed_runtime
ensure_fixed_runtime(__file__)
from dependency_resolver import bootstrap_from_requirements
bootstrap_from_requirements(Path(__file__).resolve().parents[1] / "requirements.txt")

from _runtime import contract, build_manifest
import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle


mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["axes.unicode_minus"] = False


RICE_CHR_LEN = {
    "01": 43270923,
    "02": 35937250,
    "03": 36513819,
    "04": 35502694,
    "05": 29958434,
    "06": 31248788,
    "07": 29697621,
    "08": 28443022,
    "09": 23404648,
    "10": 23207287,
    "11": 29021106,
    "12": 27531856,
}


def skill_name(input_path, out_dir="outputs", **kwargs):
    return igci_evaluation(input_path, out_dir=out_dir, **kwargs)


def workspace_root():
    return Path.cwd()


def ensure_out_dir(out_dir: str | Path) -> Path:
    path = Path(out_dir)
    if not path.is_absolute():
        path = workspace_root() / path
    path.mkdir(parents=True, exist_ok=True)
    return path


def normalize_chrom(value):
    text = str(value).strip()
    if text.lower().startswith("chr"):text=text[3:]
    if text.isdecimal() and 1<=int(text)<=12:return f"{int(text):02d}"
    return text


def normalize_allele(value: Any) -> str:
    return str(value).strip().upper()


def is_missing_gt(gt: Any) -> bool:
    text = str(gt).strip()
    return text in {"./.", ".|.", "", "."}


def gt_state(gt: Any) -> str:
    text = str(gt).split(":")[0].strip()
    if text in {"0/0", "0|0"}:
        return "REF"
    if text in {"1/1", "1|1"}:
        return "ALT"
    if text in {"0/1", "1/0", "0|1", "1|0"}:
        return "HET"
    return "MISSING"


def load_panel(panel_file: str | Path) -> pd.DataFrame:
    """Load an indica/japonica diagnostic panel.

    Production Geng-BEVA panels use ``chr,POS,Geng,Xian``, where Geng is
    the japonica diagnostic allele and Xian is the indica diagnostic allele.
    The legacy ``chr,POS,Ref,Alt`` schema is also accepted for backward
    compatibility, with Ref interpreted as japonica and Alt as indica.
    """
    panel = pd.read_csv(panel_file)
    base_required = {"chr", "POS"}
    missing_base = base_required.difference(panel.columns)
    if missing_base:
        raise ValueError(f"Panel file missing required columns: {sorted(missing_base)}")

    if {"Geng", "Xian"}.issubset(panel.columns):
        japonica_col, indica_col = "Geng", "Xian"
        schema = "Geng/Xian"
    elif {"Ref", "Alt"}.issubset(panel.columns):
        japonica_col, indica_col = "Ref", "Alt"
        schema = "Ref/Alt (legacy)"
    else:
        raise ValueError(
            "Panel must contain either chr,POS,Geng,Xian (production schema) "
            "or chr,POS,Ref,Alt (legacy schema)."
        )

    panel = panel.copy()
    panel["CHROM"] = panel["chr"].map(normalize_chrom)
    panel["POS"] = pd.to_numeric(panel["POS"], errors="raise").astype(int)
    panel["JAPONICA_ALLELE"] = panel[japonica_col].map(normalize_allele)
    panel["INDICA_ALLELE"] = panel[indica_col].map(normalize_allele)

    valid = {"A", "C", "G", "T"}
    bad = ~panel["JAPONICA_ALLELE"].isin(valid) | ~panel["INDICA_ALLELE"].isin(valid)
    if bad.any():
        raise ValueError("Panel contains non-ACGT diagnostic alleles.")
    if (panel["JAPONICA_ALLELE"] == panel["INDICA_ALLELE"]).any():
        raise ValueError("Panel contains sites where japonica and indica diagnostic alleles are identical.")

    panel = panel[["CHROM", "POS", "JAPONICA_ALLELE", "INDICA_ALLELE"]]
    if panel.duplicated(["CHROM", "POS"]).any():
        duplicated = panel[panel.duplicated(["CHROM", "POS"], keep=False)]
        conflicting = (
            duplicated.groupby(["CHROM", "POS"])[["JAPONICA_ALLELE", "INDICA_ALLELE"]]
            .nunique()
            .max(axis=1)
            .gt(1)
            .any()
        )
        if conflicting:
            raise ValueError("Conflicting allele definitions for one panel coordinate")
        panel = panel.drop_duplicates(["CHROM", "POS"])
    panel.attrs["schema"] = schema
    return panel


def load_vcf(vcf_file: str | Path) -> pd.DataFrame:
    path=Path(vcf_file)
    raw = gzip.open(path,"rb").read() if path.suffix.lower()==".gz" else path.read_bytes()
    text: Optional[str] = None
    for encoding in ("utf-8-sig", "gb18030", "gbk", "cp936", "latin1"):
        try:
            text = raw.decode(encoding)
            break
        except UnicodeDecodeError:
            continue
    if text is None:
        raise ValueError(f"Unable to decode VCF file: {vcf_file}")

    header: Optional[List[str]] = None
    records: List[List[str]] = []
    for line in text.splitlines():
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            header = line.split("\t")
            continue
        if header is not None and line.strip():
            records.append(line.split("\t"))
    if header is None:
        raise ValueError("VCF header not found")
    vcf = pd.DataFrame(records, columns=header)
    vcf = vcf.rename(columns={"#CHROM": "CHROM"})
    vcf["CHROM"] = vcf["CHROM"].map(normalize_chrom)
    vcf["POS"] = vcf["POS"].astype(int)
    vcf["REF"] = vcf["REF"].map(normalize_allele)
    vcf["ALT"] = vcf["ALT"].map(normalize_allele)
    return vcf


def sample_columns(vcf: pd.DataFrame) -> List[str]:
    excluded = {
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
        "FORMAT",
    }
    return [c for c in vcf.columns if c not in excluded]


def orient_ancestry_genotype(
    japonica_allele: str,
    indica_allele: str,
    vcf_ref: str,
    vcf_alt: str,
    gt: Any,
) -> Dict[str, Any]:
    """Orient one VCF genotype to the panel-defined Geng/Xian alleles.

    The panel defines biological ancestry. VCF REF/ALT only defines local encoding.
    """
    raw_gt = str(gt).split(":")[0].strip()
    state = gt_state(raw_gt)
    result = {
        "VCF_GT": raw_gt,
        "ALIGNMENT": "missing",
        "ANCESTRY_GENOTYPE": "MISSING",
        "XIAN_ALLELE_DOSAGE": np.nan,
        "SCORE": np.nan,
    }
    if state == "MISSING":
        return result

    exact = japonica_allele == vcf_ref and indica_allele == vcf_alt
    swapped = japonica_allele == vcf_alt and indica_allele == vcf_ref
    if not exact and not swapped:
        result["ALIGNMENT"] = "mismatched"
        result["ANCESTRY_GENOTYPE"] = "UNALIGNABLE"
        return result

    result["ALIGNMENT"] = "exact" if exact else "swapped"
    if state == "HET":
        dosage = 1
        ancestry = "GENG_XIAN_HET"
    elif (exact and state == "ALT") or (swapped and state == "REF"):
        dosage = 2
        ancestry = "XIAN_HOM"
    else:
        dosage = 0
        ancestry = "GENG_HOM"

    result["ANCESTRY_GENOTYPE"] = ancestry
    result["XIAN_ALLELE_DOSAGE"] = dosage
    result["SCORE"] = dosage / 2.0
    return result


def support_score(japonica_allele: str, indica_allele: str, vcf_ref: str, vcf_alt: str, gt: Any) -> Tuple[Optional[float], str]:
    """Backward-compatible helper returning site IGCI contribution and alignment."""
    oriented = orient_ancestry_genotype(japonica_allele, indica_allele, vcf_ref, vcf_alt, gt)
    score = oriented["SCORE"]
    return (None if pd.isna(score) else float(score), str(oriented["ALIGNMENT"]))


def read_input(vcf_file: str | Path, panel_file: str | Path) -> Tuple[pd.DataFrame, pd.DataFrame, List[str]]:
    panel = load_panel(panel_file)
    vcf = load_vcf(vcf_file)
    samples = sample_columns(vcf)
    return vcf, panel, samples


def evaluate_sample(vcf: pd.DataFrame, panel: pd.DataFrame, sample: str, window=500_000, step=100_000) -> Dict[str, Any]:
    panel_dict = {
        (row.CHROM, int(row.POS)): (row.JAPONICA_ALLELE, row.INDICA_ALLELE)
        for row in panel.itertuples(index=False)
    }

    rows: List[Dict[str, Any]] = []
    for _, row in vcf.iterrows():
        key = (row.CHROM, int(row.POS))
        prior = panel_dict.get(key)
        if prior is None:
            continue

        japonica_allele, indica_allele = prior
        gt = dict(zip(str(row.FORMAT).split(":"), str(row[sample]).split(":"))).get("GT", "./.")
        oriented = orient_ancestry_genotype(japonica_allele, indica_allele, row.REF, row.ALT, gt)
        if pd.isna(oriented["SCORE"]):
            continue

        rows.append(
            {
                "CHROM": row.CHROM,
                "POS": int(row.POS),
                "GENG_ALLELE": japonica_allele,
                "XIAN_ALLELE": indica_allele,
                "VCF_REF": row.REF,
                "VCF_ALT": row.ALT,
                "VCF_GT": oriented["VCF_GT"],
                "SAMPLE": sample,
                "ALIGNMENT": oriented["ALIGNMENT"],
                "ANCESTRY_GENOTYPE": oriented["ANCESTRY_GENOTYPE"],
                "XIAN_ALLELE_DOSAGE": int(oriented["XIAN_ALLELE_DOSAGE"]),
                "SCORE": float(oriented["SCORE"]),
            }
        )

    score_df = pd.DataFrame(
        rows,
        columns=[
            "CHROM", "POS", "GENG_ALLELE", "XIAN_ALLELE",
            "VCF_REF", "VCF_ALT", "VCF_GT", "SAMPLE", "ALIGNMENT",
            "ANCESTRY_GENOTYPE", "XIAN_ALLELE_DOSAGE", "SCORE",
        ],
    )
    if score_df.empty:
        chrom_summary = pd.DataFrame(columns=["SAMPLE", "CHROM", "INFORMATIVE_SITES", "SUPPORT_SUM", "IGCI"])
        genome_summary = pd.DataFrame(
            [{"SAMPLE": sample, "INFORMATIVE_SITES": 0, "SUPPORT_SUM": 0.0, "IGCI": np.nan}]
        )
        window_df = compute_windows(score_df, window, step)
        return {
            "sample": sample,
            "score_table": score_df,
            "chrom_summary": chrom_summary,
            "genome_summary": genome_summary,
            "window_table": window_df,
        }

    chrom_summary = (
        score_df.groupby("CHROM", as_index=False)
        .agg(
            INFORMATIVE_SITES=("SCORE", "size"),
            SUPPORT_SUM=("SCORE", "sum"),
        )
        .assign(IGCI=lambda d: d["SUPPORT_SUM"] / d["INFORMATIVE_SITES"])
    )
    chrom_summary.insert(0, "SAMPLE", sample)

    genome_support = float(score_df["SCORE"].sum())
    genome_sites = int(score_df.shape[0])
    genome_igci = genome_support / genome_sites if genome_sites else np.nan
    genome_summary = pd.DataFrame(
        [
            {
                "SAMPLE": sample,
                "INFORMATIVE_SITES": genome_sites,
                "SUPPORT_SUM": genome_support,
                "IGCI": genome_igci,
            }
        ]
    )

    window_df = compute_windows(score_df, window, step)
    return {
        "sample": sample,
        "score_table": score_df,
        "chrom_summary": chrom_summary,
        "genome_summary": genome_summary,
        "window_table": window_df,
    }


def compute_windows(score_df: pd.DataFrame, window: int = 500_000, step: int = 100_000) -> pd.DataFrame:
    if window <= 0 or step <= 0 or step > window:raise ValueError("Require 0 < step <= window")
    rows: List[Dict[str, Any]] = []
    for chrom in RICE_CHR_LEN:
        chr_df = score_df[score_df["CHROM"] == chrom]
        max_pos = RICE_CHR_LEN.get(chrom)
        if max_pos is None:
            continue
        for start in np.arange(1, max_pos, step):
            start_i = int(start)
            end_i = min(start_i + window, max_pos + 1)
            sub = chr_df[(chr_df["POS"] >= start_i) & (chr_df["POS"] < end_i)]
            if sub.empty:
                score = np.nan
                count = 0
            else:
                score = float(sub["SCORE"].mean())
                count = int(sub.shape[0])
            rows.append(
                {
                    "CHROM": chrom,
                    "START": start_i,
                    "END": end_i,
                    "SCORE": score,
                    "COUNT": count,
                    "DISPLAY_END": min(start_i + step, max_pos + 1),
                }
            )
    return pd.DataFrame(rows)


def plot_heatmap(window_df: pd.DataFrame, sample: str, out_dir: Path, value_col: str = "SCORE") -> Optional[Path]:
    if window_df.empty:
        return None

    cmap = LinearSegmentedColormap.from_list(
        "igci",
        ["#8EB4D0", "#BDD7E7", "#FCE6C8", "#F4AF64"],
    )
    fig, ax = plt.subplots(figsize=(14, 5))
    norm = mpl.colors.Normalize(vmin=0, vmax=1)

    chr_list = [c for c in RICE_CHR_LEN if c in set(window_df["CHROM"].astype(str))]
    chr_map = {c: len(chr_list) - 1 - i for i, c in enumerate(chr_list)}

    for _, row in window_df.iterrows():
        chrom = row["CHROM"]
        y = chr_map[chrom]
        value = float(row[value_col])
        x0 = float(row["START"]) / 1e6
        width = (float(row["DISPLAY_END"]) - float(row["START"])) / 1e6
        rect = Rectangle(
            (x0, y - 0.4),
            width,
            0.8,
            facecolor=cmap(norm(value)) if np.isfinite(value) else "#D9D9D9",
            edgecolor="none",
        )
        ax.add_patch(rect)

    genome_max = max(RICE_CHR_LEN.values()) / 1e6
    ax.set_xlim(0, genome_max)
    ax.set_ylim(-0.5, len(chr_list) - 0.5)
    ax.set_yticks(range(len(chr_list)))
    ax.set_yticklabels(chr_list[::-1], fontsize=11)
    ax.set_xlabel("Physical position (Mb)", fontsize=12)
    ax.set_ylabel("Chromosome", fontsize=12)
    ax.set_title(f"{sample} single-sample indica introgression profile", fontsize=14)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cb = plt.colorbar(sm, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("IGCI (indica allele dosage proportion)")

    plt.tight_layout()
    out_png = out_dir / f"{sample}_IGCI_introgression_profile.png"
    out_pdf = out_dir / f"{sample}_IGCI_introgression_profile.pdf"
    plt.savefig(out_png, dpi=600, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    return out_png


@contract
def igci_evaluation(
    input_path,
    out_dir: str | Path = "outputs",
    **kwargs,
) -> Dict[str, Any]:
    input_path = Path(input_path)
    out_path = ensure_out_dir(out_dir)
    panel_file = kwargs.get("panel_file") or kwargs.get("panel")
    if panel_file is None:
        bundled_panel = Path(__file__).resolve().parents[1] / "assets" / "IJ_panel.csv"
        beside_vcf = input_path.parent / "IJ_panel.csv"
        if bundled_panel.exists():
            panel_file = bundled_panel
        elif beside_vcf.exists():
            panel_file = beside_vcf
        else:
            raise FileNotFoundError(
                "No IJ panel was supplied and no bundled assets/IJ_panel.csv was found. "
                "Provide --panel /path/to/IJ_panel.csv."
            )
    sample = kwargs.get("sample") or kwargs.get("sample_id") or kwargs.get("target_sample")

    vcf, panel, samples = read_input(input_path, panel_file)
    targets = [str(sample)] if sample else samples

    manifests: List[Dict[str, Any]] = []
    score_tables: List[pd.DataFrame] = []
    chrom_tables: List[pd.DataFrame] = []
    genome_tables: List[pd.DataFrame] = []

    for target in targets:
        if target not in samples:
            raise ValueError(f"Sample not found in VCF header: {target}")
        from urllib.parse import quote
        safe_target = quote(target, safe="") or "sample"
        result = evaluate_sample(vcf, panel, target, int(kwargs.get("window",500_000)), int(kwargs.get("step",100_000)))

        score_csv = out_path / f"{safe_target}_igci_scores.csv"
        chrom_csv = out_path / f"{safe_target}_igci_chromosome_summary.csv"
        genome_csv = out_path / f"{safe_target}_igci_genome_summary.csv"
        window_csv = out_path / f"{safe_target}_igci_windows.csv"

        result["score_table"].to_csv(score_csv, index=False, encoding="utf-8-sig")
        result["chrom_summary"].to_csv(chrom_csv, index=False, encoding="utf-8-sig")
        result["genome_summary"].to_csv(genome_csv, index=False, encoding="utf-8-sig")
        result["window_table"].to_csv(window_csv, index=False, encoding="utf-8-sig")

        plot_png = None
        if kwargs.get("plot",True):
            sample_plot_dir=out_path / safe_target
            sample_plot_dir.mkdir(exist_ok=True)
            plot_png = plot_heatmap(result["window_table"], safe_target, sample_plot_dir)

        manifests.append(
            {
                "sample": target,
                "score_csv": str(score_csv),
                "chromosome_csv": str(chrom_csv),
                "genome_csv": str(genome_csv),
                "window_csv": str(window_csv),
                "plot_png": str(plot_png) if plot_png else None,
                "plot_pdf": str(plot_png.with_suffix(".pdf")) if plot_png else None,
            }
        )
        score_tables.append(result["score_table"])
        chrom_tables.append(result["chrom_summary"])
        genome_tables.append(result["genome_summary"])

    combined_scores = pd.concat(score_tables, ignore_index=True) if score_tables else pd.DataFrame()
    combined_chrom = pd.concat(chrom_tables, ignore_index=True) if chrom_tables else pd.DataFrame()
    combined_genome = pd.concat(genome_tables, ignore_index=True) if genome_tables else pd.DataFrame()

    combined_scores_csv = out_path / "igci_all_samples_scores.csv"
    combined_chrom_csv = out_path / "igci_all_samples_chromosome_summary.csv"
    combined_genome_csv = out_path / "igci_all_samples_genome_summary.csv"

    if not combined_scores.empty:
        combined_scores.to_csv(combined_scores_csv, index=False, encoding="utf-8-sig")
    if not combined_chrom.empty:
        combined_chrom.to_csv(combined_chrom_csv, index=False, encoding="utf-8-sig")
    if not combined_genome.empty:
        combined_genome.to_csv(combined_genome_csv, index=False, encoding="utf-8-sig")

    outputs = {}
    for index, item in enumerate(manifests):
        for key in ("score_csv","chromosome_csv","genome_csv","window_csv","plot_png","plot_pdf"):
            if item.get(key): outputs[f"sample_{index}_{key}"] = item[key]
    for key,path,table in [("result_csv",combined_scores_csv,combined_scores),("chromosome_csv",combined_chrom_csv,combined_chrom),("summary_csv",combined_genome_csv,combined_genome)]:
        if not table.empty:outputs[key]=str(path)
    return {
        "module":"igci_evaluation", "status":"success", "outputs":outputs,
        "warnings": (["One or more samples have no informative sites; IGCI is missing, not zero."] if (combined_genome["INFORMATIVE_SITES"]==0).any() else []),
        "metrics":{"samples":len(targets),"window_bp":int(kwargs.get("window",500_000)),"step_bp":int(kwargs.get("step",100_000))},
        "input": str(input_path),
        "panel": str(panel_file),
        "out_dir": str(out_path),
        "samples": targets,
        "manifests": manifests,
        "combined_scores_csv": str(combined_scores_csv) if not combined_scores.empty else None,
        "combined_chromosome_csv": str(combined_chrom_csv) if not combined_chrom.empty else None,
        "combined_genome_csv": str(combined_genome_csv) if not combined_genome.empty else None,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="IGCI evaluation with IJ prior panel")
    parser.add_argument("--vcf", required=True, help="Input VCF file")
    parser.add_argument("--panel", default=None, help="IJ_panel.csv path")
    parser.add_argument("--out-dir", default="outputs", help="Output directory")
    parser.add_argument("--sample", default=None, help="Single sample to evaluate")
    parser.add_argument("--window",type=int,default=500000)
    parser.add_argument("--step",type=int,default=100000)
    parser.add_argument("--no-plot",action="store_true")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    result = igci_evaluation(
        args.vcf,
        out_dir=args.out_dir,
        panel_file=args.panel,
        sample=args.sample, window=args.window, step=args.step, plot=not args.no_plot,
    )
    print(
        pd.Series(
            {
                "input": result["input"],
                "panel": result["panel"],
                "out_dir": result["out_dir"],
                "samples": ", ".join(result["samples"]),
            }
        ).to_string()
    )
    return 0


def run(input_path,out_dir="outputs",**kwargs):
    return igci_evaluation(input_path,out_dir=out_dir,**kwargs)


if __name__ == "__main__":
    raise SystemExit(main())
