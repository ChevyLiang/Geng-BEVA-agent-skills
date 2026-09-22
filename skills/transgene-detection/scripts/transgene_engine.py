from __future__ import annotations
import gzip
import math
import re
from pathlib import Path
import pandas as pd

PATTERNS = [
    (r"(?:^|[^a-z0-9])(?:pcamv35s|pfmv35s|camv35s|fmv35s|35s)(?:$|[^a-z0-9])", "35S"),
    (r"(?:^|[^a-z0-9])(?:bt_resistant|bt)(?:$|[^a-z0-9])", "BT"),
    (r"(?:^|[^a-z0-9])cp4[_-]?epsps(?:$|[^a-z0-9])", "CP4_EPSPS"),
    (r"(?:^|[^a-z0-9])caoganlin(?:$|[^a-z0-9])", "caoganlin"),
    (r"(?:^|[^a-z0-9])tnos(?:$|[^a-z0-9])", "tNOS"),
    (r"(?:^|[^a-z0-9])nos(?:$|[^a-z0-9])", "Nos"),
    (r"(?:^|[^a-z0-9])ubiquitin(?:$|[^a-z0-9])", "ubiquitin"),
]
SUPPORT = {"low_support", "moderate_support", "high_support"}
LEVEL_ORDER = {"no_data":0, "insufficient_support":1, "low_support":2, "moderate_support":3, "high_support":4}


def open_vcf(path):
    path=Path(path)
    return gzip.open(path,"rt",encoding="utf-8-sig") if path.suffix.lower()==".gz" else path.open(encoding="utf-8-sig")


def canonical(chrom):
    text=re.sub(r"^chr","",str(chrom),flags=re.I)
    return text.isdecimal() and 1<=int(text)<=12


def canonical_transgene_element_name(value):
    for pattern,family in PATTERNS:
        if re.search(pattern,str(value),re.I): return family
    return None


def matches_transgene_record(chrom, vid):
    # Project-specific contract: non-1..12 data rows represent assay elements.
    # Metadata lines are never passed into this function.
    if canonical(chrom): return None
    return canonical_transgene_element_name(chrom) or canonical_transgene_element_name(vid) or str(chrom)


def is_called_genotype(gt):
    return str(gt).strip() in {"0/0","0|0","0/1","1/0","0|1","1|0","1/1","1|1"}


def classify_transgene_dp(dp, genotype=None, low_dp=10, moderate_dp=50, high_dp=100):
    if genotype is not None and not is_called_genotype(genotype): return "no_data"
    try: value=float(dp)
    except (TypeError,ValueError): return "no_data"
    if not math.isfinite(value) or value<0: return "no_data"
    if not (0 <= float(low_dp) <= float(moderate_dp) <= float(high_dp)):
        raise ValueError("Require 0 <= low_dp <= moderate_dp <= high_dp")
    if value>float(high_dp):return "high_support"
    if value>=float(moderate_dp):return "moderate_support"
    if value>=float(low_dp):return "low_support"
    return "insufficient_support"


def evaluate(vcf_path, sample_id=None, low_dp=10, moderate_dp=50, high_dp=100):
    samples=None;records=[];element_count=0
    with open_vcf(vcf_path) as src:
        for line in src:
            if line.startswith("##"):continue
            if line.startswith("#CHROM"):
                samples=line.rstrip("\r\n").split("\t")[9:]
                if len(samples)!=len(set(samples)): raise ValueError("Duplicate VCF sample IDs")
                if sample_id is not None and str(sample_id) not in samples: raise ValueError(f"Unknown sample ID: {sample_id}")
                continue
            if line.startswith("#") or not line.strip():continue
            if samples is None:raise ValueError("Missing VCF #CHROM header")
            fields=line.rstrip("\r\n").split("\t")
            if len(fields)!=9+len(samples):raise ValueError("VCF row width differs from header")
            chrom,pos,vid=fields[:3]
            family=matches_transgene_record(chrom,vid)
            if family is None:continue
            element_count+=1
            locus=vid if vid!='.' else f"{chrom}_{pos}"
            for sample,cell in zip(samples,fields[9:]):
                if sample_id is not None and sample!=str(sample_id):continue
                fm=dict(zip(fields[8].split(':'),cell.split(':')))
                gt=fm.get('GT','./.');dp=fm.get('DP')
                level=classify_transgene_dp(dp,gt,low_dp=low_dp,moderate_dp=moderate_dp,high_dp=high_dp)
                try:
                    number=float(dp)
                    if not math.isfinite(number) or number<0:number=float('nan')
                except (TypeError,ValueError):number=float('nan')
                records.append(dict(SampleID=sample,CHROM=chrom,POS=pos,RecordID=locus,TransgeneFamily=family,
                    Genotype=gt,SampleDP=number,EvidenceLevel=level,HasSupport=level in SUPPORT))
    if samples is None:raise ValueError('Missing VCF #CHROM header')
    cols=['SampleID','CHROM','POS','RecordID','TransgeneFamily','Genotype','SampleDP','EvidenceLevel','HasSupport']
    detail=pd.DataFrame(records,columns=cols)
    summaries=[]
    for sample in samples:
        if sample_id is not None and sample!=str(sample_id):continue
        rows=detail.loc[detail.SampleID==sample]
        valid=rows.loc[rows.EvidenceLevel!='no_data']
        hits=rows.loc[rows.HasSupport.astype(bool)]
        level=max(rows.EvidenceLevel,key=LEVEL_ORDER.get) if not rows.empty else 'no_transgene_information'
        summaries.append(dict(SampleID=sample,HasTransgeneInformation=element_count>0,
            HasTransgeneSupport=not hits.empty,EvidenceLevel=level,
            MaxValidDP=valid.SampleDP.max() if not valid.empty else float('nan'),
            SupportingRecordCount=len(hits),HighSupportRecordCount=int((rows.EvidenceLevel=='high_support').sum()),
            ModerateSupportRecordCount=int((rows.EvidenceLevel=='moderate_support').sum()),
            LowSupportRecordCount=int((rows.EvidenceLevel=='low_support').sum()),
            InsufficientSupportRecordCount=int((rows.EvidenceLevel=='insufficient_support').sum()),
            NoDataRecordCount=int((rows.EvidenceLevel=='no_data').sum()),
            SupportingFamilies=','.join(sorted(set(hits.TransgeneFamily))),
            SupportingRecords=','.join(hits.RecordID.astype(str))))
    return detail,pd.DataFrame(summaries),element_count


def summarize_transgene_sample_events(vcf_path):
    return evaluate(vcf_path)[1]


def summarize_transgene_element_events(vcf_path):
    return evaluate(vcf_path)[0]
