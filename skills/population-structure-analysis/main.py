from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_CODEX_ROOT = next(p for p in Path(__file__).resolve().parents if p.name == ".codex")
sys.path.insert(0, str(_CODEX_ROOT / "tools"))
from runtime_launcher import ensure_fixed_runtime
ensure_fixed_runtime(__file__)
from dependency_resolver import ensure_python_packages

# Core scientific packages are self-healed before analysis modules are imported.
# Plotly is resolved lazily only for PCA/IBS interactive output.
ensure_python_packages(["numpy", "pandas", "matplotlib"], start=__file__)

from scripts.fst_analysis import fst_analysis
from scripts.genotype_frequency import genotype_frequency
from scripts.ibs_analysis import ibs_analysis
from scripts.individual_diversity import individual_diversity
from scripts.nontransgene_vcf_export import nontransgene_vcf_export
from scripts.pca_analysis import pca_analysis
from scripts.pairwise_matrix import pairwise_matrix
from scripts.sample_name_explorer import sample_name_explorer
from scripts.sample_pairwise_compare import sample_pairwise_compare
from scripts.transgene_filter import transgene_filter
from scripts.transgene_event_summary import transgene_event_summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Population Structure Agent")
    sub = parser.add_subparsers(dest="command", required=True)

    def add_common(p):
        p.add_argument("input_path")
        p.add_argument("--out-dir", default="outputs")
        p.add_argument("--plink-cmd", default=None)

    for name in [
        "sample-names",
        "genotype-frequency",
        "transgene-filter",
        "transgene-events",
        "nontransgene-vcf",
        "ibs",
        "ibd",
        "diversity",
        "genotype-summary",
    ]:
        p = sub.add_parser(name)
        add_common(p)
        if name == "transgene-events":
            p.add_argument("--sample-id", default=None, help="Optional focus sample for sample-specific summary")

    p_pair = sub.add_parser("pairwise-compare")
    add_common(p_pair)
    p_pair.add_argument("--sample-a", required=True)
    p_pair.add_argument("--sample-b", required=True)

    p_matrix = sub.add_parser("pairwise-matrix")
    add_common(p_matrix)
    p_matrix.add_argument("--samples", nargs="+", required=True)

    p_pca = sub.add_parser("pca")
    add_common(p_pca)
    p_pca.add_argument("--group-file", default=None, help="CSV/TSV with SampleID and Group columns")
    p_pca.add_argument("--sample-id", default=None, help="Focus sample for nearest/farthest summary")

    p_fst = sub.add_parser("fst")
    add_common(p_fst)
    p_fst.add_argument("--group-a", nargs="+", required=True)
    p_fst.add_argument("--group-b", nargs="+", required=True)
    p_fst.add_argument("--note", default=None, help="Optional label for the comparison")
    p_fst.add_argument("--sliding-window", action="store_true", help="Also calculate sliding-window FST. If window/step are omitted, defaults are 500 kb/100 kb.")
    p_fst.add_argument("--window", type=int, default=None, help="Sliding-window size in bp; implies --sliding-window. Default when requested: 500000")
    p_fst.add_argument("--step", type=int, default=None, help="Sliding-window step in bp; implies --sliding-window. Default when requested: 100000")
    p_fst.add_argument("--min-sites", type=int, default=1, help="Minimum valid sites required for a sliding-window estimate")

    args = parser.parse_args()
    kwargs = vars(args).copy()
    cmd = kwargs.pop("command")

    manifest = run(command=cmd, **kwargs)

    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


def run(input_path, out_dir="outputs", **kwargs):
    command = kwargs.pop("command", "pca")
    functions = {
        "sample-names": sample_name_explorer, "genotype-frequency": genotype_frequency,
        "transgene-filter": transgene_filter, "transgene-events": transgene_event_summary,
        "nontransgene-vcf": nontransgene_vcf_export, "pca": pca_analysis,
        "ibs": ibs_analysis, "ibd": ibs_analysis, "diversity": individual_diversity, "genotype-summary": individual_diversity,
        "pairwise-compare": sample_pairwise_compare, "pairwise-matrix": pairwise_matrix,
        "fst": fst_analysis,
    }
    if command not in functions: raise ValueError(f"Unknown command: {command}")
    if command in {"ibs", "ibd"}: kwargs["mode"] = command
    return functions[command](input_path, out_dir=out_dir, **kwargs)


if __name__ == "__main__":
    raise SystemExit(main())

