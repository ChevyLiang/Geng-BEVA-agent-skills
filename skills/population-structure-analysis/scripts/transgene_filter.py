from scripts.transgene_event_summary import transgene_event_summary
from scripts.nontransgene_vcf_export import nontransgene_vcf_export
from core.runtime import contract


@contract
def transgene_filter(input_path,out_dir="outputs",**kwargs):
    summary=transgene_event_summary(input_path,out_dir=out_dir,**kwargs)
    cleaned=nontransgene_vcf_export(input_path,out_dir=out_dir)
    summary["outputs"].update(cleaned_vcf=cleaned["outputs"]["result_vcf"],removed_csv=cleaned["outputs"]["result_csv"])
    return summary


def run(input_path,out_dir="outputs",**kwargs):return transgene_filter(input_path,out_dir=out_dir,**kwargs)
