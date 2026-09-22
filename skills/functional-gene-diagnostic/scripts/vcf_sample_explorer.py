from __future__ import annotations
from _runtime import contract

import argparse
import gzip
import json
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

import pandas as pd


def _ensure_out_dir(out_dir: str | Path) -> Path:
    path = Path(out_dir)
    if not path.is_absolute():
        path = Path.cwd() / path
    path.mkdir(parents=True, exist_ok=True)
    return path


def _open_text(path: Path):
    if path.suffix.lower() == ".gz" or str(path).lower().endswith(".vcf.gz"):
        return gzip.open(path, "rt", encoding="utf-8-sig")
    return path.open("r", encoding="utf-8-sig")


def _extract_sample_ids(vcf_path: str | Path) -> List[str]:
    vcf_path = Path(vcf_path)
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF not found: {vcf_path}")

    with _open_text(vcf_path) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= 9:
                    return []
                return parts[9:]
    raise ValueError(f"VCF header does not contain a #CHROM line: {vcf_path}")


def _apply_filters(sample_ids: Iterable[str], *, contains=None, prefix=None, regex=None) -> List[str]:
    items = list(sample_ids)
    if contains:
        tokens = [str(x) for x in contains] if isinstance(contains, (list, tuple, set)) else [str(contains)]
        items = [sid for sid in items if any(token in sid for token in tokens)]
    if prefix:
        tokens = [str(x) for x in prefix] if isinstance(prefix, (list, tuple, set)) else [str(prefix)]
        items = [sid for sid in items if any(sid.startswith(token) for token in tokens)]
    if regex:
        pattern = re.compile(str(regex))
        items = [sid for sid in items if pattern.search(sid)]
    return items


@contract
def vcf_sample_explorer(input_path, out_dir="outputs", **kwargs) -> Dict[str, Any]:
    out_path = _ensure_out_dir(out_dir)
    sample_ids = _extract_sample_ids(input_path)
    filtered_ids = _apply_filters(
        sample_ids,
        contains=kwargs.get("contains"),
        prefix=kwargs.get("prefix"),
        regex=kwargs.get("regex"),
    )

    result = pd.DataFrame(
        {
            "SampleIndex": range(1, len(sample_ids) + 1),
            "SampleID": sample_ids,
        }
    )
    if filtered_ids != sample_ids:
        result = result[result["SampleID"].isin(filtered_ids)].copy()

    result_csv = out_path / "vcf_sample_explorer.csv"
    result.to_csv(result_csv, index=False, encoding="utf-8-sig")

    summary_csv = out_path / "vcf_sample_explorer_summary.csv"
    summary = pd.DataFrame(
        [
            {
                "VCF": str(input_path),
                "TotalSamples": len(sample_ids),
                "MatchedSamples": len(result),
            }
        ]
    )
    summary.to_csv(summary_csv, index=False, encoding="utf-8-sig")

    return {
        "module": "vcf_sample_explorer",
        "input_path": str(input_path),
        "out_dir": str(out_path),
        "status": "success",
        "outputs": {
            "result_csv": str(result_csv),
            "summary_csv": str(summary_csv),
            "heatmap_png": None,
        },
        "filters": {
            "contains": kwargs.get("contains"),
            "prefix": kwargs.get("prefix"),
            "regex": kwargs.get("regex"),
        },
    }


def run(input_path, out_dir="outputs", **kwargs):
    return vcf_sample_explorer(input_path, out_dir=out_dir, **kwargs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract sample IDs from VCF header")
    parser.add_argument("input_path")
    parser.add_argument("--out-dir", default="outputs")
    parser.add_argument("--contains", action="append", default=None)
    parser.add_argument("--prefix", action="append", default=None)
    parser.add_argument("--regex", default=None)
    args = parser.parse_args()

    manifest = vcf_sample_explorer(
        args.input_path,
        out_dir=args.out_dir,
        contains=args.contains,
        prefix=args.prefix,
        regex=args.regex,
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))

