from __future__ import annotations
from core.runtime import contract

from pathlib import Path
from typing import Any, Dict

import pandas as pd

from core.manifest import build_manifest
from tools.plink_wrapper import ensure_out_dir


def _vcf_to_genotype_matrix(vcf_path: Path) -> pd.DataFrame:
    samples = []
    records = []
    with vcf_path.open("r", encoding="utf-8-sig") as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                parts = line.rstrip("\n").split("\t")
                samples = parts[9:]
                continue
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            vid = parts[2] if parts[2] != "." else f"{parts[0]}_{parts[1]}"
            gts = []
            for sample_field in parts[9:]:
                gt = sample_field.split(":", 1)[0]
                if gt in {"", ".", "./.", ".|."}:
                    gts.append(None)
                else:
                    gt = gt.replace("|", "/")
                    if gt == "0/0":
                        gts.append(0)
                    elif gt in {"0/1", "1/0"}:
                        gts.append(1)
                    elif gt == "1/1":
                        gts.append(2)
                    else:
                        gts.append(None)
            records.append((vid, gts))

    if not samples:
        raise ValueError("VCF header missing sample IDs")
    if not records:
        raise ValueError("No variants found in VCF")

    matrix = pd.DataFrame(index=samples)
    for vid, gts in records:
        matrix[vid] = gts
    return matrix


def _individual_diversity(matrix: pd.DataFrame) -> pd.DataFrame:
    diversity = pd.DataFrame(index=matrix.index)
    diversity["ObservedVariants"] = matrix.notna().sum(axis=1)
    diversity["HeterozygousCount"] = (matrix == 1).sum(axis=1)
    diversity["HomoAltCount"] = (matrix == 2).sum(axis=1)
    diversity["MissingCount"] = matrix.isna().sum(axis=1)
    diversity["AltAlleleDosageSum"] = diversity["HeterozygousCount"] + 2 * diversity["HomoAltCount"]
    return diversity.reset_index(names="SampleID")


@contract
def individual_diversity(input_path, out_dir="outputs", **kwargs) -> Dict[str, Any]:
    out_path = ensure_out_dir(out_dir)
    vcf_path = Path(input_path)
    matrix = _vcf_to_genotype_matrix(vcf_path)
    diversity = _individual_diversity(matrix)

    matrix_csv = out_path / "individual_genotype_matrix.csv"
    diversity_csv = out_path / "individual_genotype_summary.csv"
    matrix.to_csv(matrix_csv, encoding="utf-8-sig")
    diversity.to_csv(diversity_csv, index=False, encoding="utf-8-sig")

    heatmap_png = out_path / "individual_diversity_heatmap.png"
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns

        plot_data = diversity.set_index("SampleID")[["ObservedVariants", "HeterozygousCount", "HomoAltCount", "MissingCount", "AltAlleleDosageSum"]]
        plt.figure(figsize=(10, max(4, len(plot_data) * 0.25)))
        sns.heatmap(plot_data, cmap="viridis", linewidths=0.1)
        plt.title("Individual Genotype Summary")
        plt.tight_layout()
        plt.savefig(heatmap_png, dpi=300, bbox_inches="tight")
        plt.close()
    except Exception:
        heatmap_png = None

    return build_manifest(
        "individual_genotype_summary",
        input_path,
        out_path,
        result_csv=diversity_csv,
        heatmap_png=heatmap_png,
        extra={"genotype_matrix_csv": str(matrix_csv)},
    )


def run(input_path, out_dir="outputs", **kwargs):
    return individual_diversity(input_path,out_dir=out_dir,**kwargs)


if __name__ == "__main__":
    raise SystemExit("Import and call individual_diversity().")

