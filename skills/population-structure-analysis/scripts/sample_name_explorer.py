from __future__ import annotations
from core.runtime import contract

from typing import Any, Dict

from core.manifest import build_manifest
from core.vcf_utils import read_vcf_samples
from tools.plink_wrapper import ensure_out_dir


@contract
def sample_name_explorer(input_path, out_dir="outputs", **kwargs) -> Dict[str, Any]:
    out_path = ensure_out_dir(out_dir)
    samples = read_vcf_samples(input_path)
    result_csv = out_path / "sample_names.csv"
    import pandas as pd

    pd.DataFrame({"SampleID": samples}).to_csv(result_csv, index=False, encoding="utf-8-sig")
    return build_manifest("sample_name_explorer", input_path, out_path, result_csv=result_csv, extra={"sample_count": len(samples)})


def run(input_path, out_dir="outputs", **kwargs):
    return sample_name_explorer(input_path,out_dir=out_dir,**kwargs)


if __name__ == "__main__":
    raise SystemExit("Import and call sample_name_explorer().")

