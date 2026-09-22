import argparse
import json
import sys
from pathlib import Path

_CODEX_ROOT = next(p for p in Path(__file__).resolve().parents if p.name == ".codex")
sys.path.insert(0, str(_CODEX_ROOT / "tools"))
from runtime_launcher import ensure_fixed_runtime
ensure_fixed_runtime(__file__)
from dependency_resolver import bootstrap_from_requirements
bootstrap_from_requirements(Path(__file__).resolve().parent / "requirements.txt")

sys.path.insert(0,str(Path(__file__).resolve().parent/"scripts"))
from single_sample_func_table import run as single
from pop_func_table import run as population
from trait_category_func_extract import run as trait
from disadvantage_target_identify import run as targets
from vcf_sample_explorer import run as samples


def run(input_path,out_dir="outputs",**kwargs):
    mode=kwargs.pop("mode","single")
    return {"single":single,"population":population,"trait":trait,"targets":targets,"samples":samples}[mode](input_path,out_dir=out_dir,**kwargs)


if __name__=="__main__":
    p=argparse.ArgumentParser(description="Functional gene diagnostic")
    p.add_argument("mode",choices=["single","population","trait","targets","samples"])
    p.add_argument("input_path");p.add_argument("--out-dir",default="outputs")
    p.add_argument("--functional-qtn-csv");p.add_argument("--plink-cmd");p.add_argument("--sample-id")
    p.add_argument("--trait");p.add_argument("--category");p.add_argument("--breeding-preferences");p.add_argument("--environment")
    print(json.dumps(run(**vars(p.parse_args())),ensure_ascii=False,indent=2))
