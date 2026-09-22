from __future__ import annotations
import argparse
import json
from pathlib import Path
from _runtime import contract, ensure_out_dir, build_manifest
from _vcf import export_canonical_vcf
from transgene_engine import evaluate


@contract
def transgene_detection(input_path,out_dir="outputs",**kwargs):
    out=ensure_out_dir(out_dir)
    details,summary,count=evaluate(input_path,kwargs.get("sample_id"),low_dp=kwargs.get("low_dp",10),moderate_dp=kwargs.get("moderate_dp",50),high_dp=kwargs.get("high_dp",100))
    detail_csv=out/"transgene_locus_evidence.csv"
    summary_csv=out/"transgene_sample_summary.csv"
    details.to_csv(detail_csv,index=False,encoding="utf-8-sig")
    summary.to_csv(summary_csv,index=False,encoding="utf-8-sig")
    clean=None
    if kwargs.get("export_clean",False): clean=export_canonical_vcf(input_path,out/"canonical.nontransgene.vcf")
    m=build_manifest("transgene_detection",input_path,out,result_csv=summary_csv,summary_csv=summary_csv,evidence_csv=detail_csv,result_vcf=clean)
    m["metrics"]={"element_records":count,"samples":len(summary),"supported_samples":int(summary.HasTransgeneSupport.sum()),"low_dp":kwargs.get("low_dp",10),"moderate_dp":kwargs.get("moderate_dp",50),"high_dp":kwargs.get("high_dp",100)}
    if not count:m["warnings"].append("No transgene-element data records; transgene status cannot be assessed from this file.")
    return m


def run(input_path,out_dir="outputs",**kwargs):return transgene_detection(input_path,out_dir=out_dir,**kwargs)


if __name__=="__main__":
    p=argparse.ArgumentParser(description="Assay-specific transgene evidence; not calibrated probability")
    p.add_argument("input_path");p.add_argument("--out-dir",default="outputs");p.add_argument("--sample-id")
    p.add_argument("--export-clean",action="store_true")
    p.add_argument("--low-dp",type=float,default=10)
    p.add_argument("--moderate-dp",type=float,default=50)
    p.add_argument("--high-dp",type=float,default=100)
    print(json.dumps(run(**vars(p.parse_args())),ensure_ascii=False,indent=2))
