from __future__ import annotations
import argparse
import json
from single_sample_func_table import single_sample_func_table


def run(input_path, out_dir="outputs", **kwargs):
    return single_sample_func_table(input_path, out_dir=out_dir, **kwargs)


def vcf_to_func_single_sample(vcf_path, out_dir="outputs", **kwargs):
    return run(vcf_path, out_dir=out_dir, **kwargs)


if __name__ == "__main__":
    p=argparse.ArgumentParser(description="Extract functional loci with PLINK 1.9 and annotate")
    p.add_argument("input_path");p.add_argument("--out-dir",default="outputs")
    p.add_argument("--functional-qtn-csv");p.add_argument("--plink-cmd");p.add_argument("--sample-id")
    p.add_argument("--breeding-preferences");p.add_argument("--environment")
    print(json.dumps(run(**vars(p.parse_args())),ensure_ascii=False,indent=2))
