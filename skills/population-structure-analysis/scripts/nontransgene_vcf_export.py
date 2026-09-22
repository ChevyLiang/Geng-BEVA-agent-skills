from __future__ import annotations
from pathlib import Path
import pandas as pd
from core.runtime import contract, build_manifest, ensure_out_dir
from core.canonical_vcf import export_canonical_vcf, open_vcf, canonical_chrom


@contract
def nontransgene_vcf_export(input_path,out_dir="outputs",**kwargs):
    out=ensure_out_dir(out_dir)
    clean=export_canonical_vcf(input_path,out/"canonical.nontransgene.vcf")
    rows=[];kept=0
    with open_vcf(input_path) as src:
        for line in src:
            if line.startswith('#') or not line.strip():continue
            parts=line.rstrip().split('\t')
            if canonical_chrom(parts[0]):kept+=1
            else:rows.append({'CHROM':parts[0],'POS':parts[1],'LocusID':parts[2],'Reason':'standard_assay_transgene_record'})
    removed=out/"nontransgene_removed_variants.csv"
    pd.DataFrame(rows,columns=['CHROM','POS','LocusID','Reason']).to_csv(removed,index=False,encoding='utf-8-sig')
    m=build_manifest('nontransgene_vcf_export',input_path,out,result_vcf=clean,result_csv=removed)
    m['metrics']={'kept_records':kept,'removed_records':len(rows)}
    return m


def run(input_path,out_dir="outputs",**kwargs):return nontransgene_vcf_export(input_path,out_dir=out_dir,**kwargs)
