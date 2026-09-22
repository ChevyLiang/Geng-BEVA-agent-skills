from core.runtime import contract, build_manifest, ensure_out_dir
from core.transgene_bridge import engine


@contract
def transgene_event_summary(input_path,out_dir="outputs",**kwargs):
    out=ensure_out_dir(out_dir)
    details,summary,count=engine().evaluate(input_path,kwargs.get("sample_id"))
    result=out/"transgene_sample_summary.csv";evidence=out/"transgene_locus_evidence.csv"
    summary.to_csv(result,index=False,encoding="utf-8-sig");details.to_csv(evidence,index=False,encoding="utf-8-sig")
    m=build_manifest("transgene_event_summary",input_path,out,result_csv=result,evidence_csv=evidence)
    m["warnings"]=["Compatibility entry point; use the transgene-detection skill for transgene evidence."]
    m["metrics"]={"element_records":count}
    return m


def run(input_path,out_dir="outputs",**kwargs):return transgene_event_summary(input_path,out_dir=out_dir,**kwargs)
