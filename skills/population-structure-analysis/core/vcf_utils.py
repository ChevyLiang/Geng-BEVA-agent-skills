from __future__ import annotations

import re
from .canonical_vcf import open_vcf
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

import pandas as pd

TRANSGENE_ELEMENT_PATTERNS = [
    (re.compile(r"(?i)(?:^|[^a-z0-9])(?:pcamv35s|pfmv35s|camv35s|fmv35s|35s)(?:$|[^a-z0-9])"), "35S"),
    (re.compile(r"(?i)(?:^|[^a-z0-9])(?:bt_resistant|bt)(?:$|[^a-z0-9])"), "BT"),
    (re.compile(r"(?i)(?:^|[^a-z0-9])cp4[_-]?epsps(?:$|[^a-z0-9])"), "CP4_EPSPS"),
    (re.compile(r"(?i)(?:^|[^a-z0-9])caoganlin(?:$|[^a-z0-9])"), "caoganlin"),
    (re.compile(r"(?i)(?:^|[^a-z0-9])tnos(?:$|[^a-z0-9])"), "tNOS"),
    (re.compile(r"(?i)(?:^|[^a-z0-9])nos(?:$|[^a-z0-9])"), "Nos"),
    (re.compile(r"(?i)(?:^|[^a-z0-9])ubiquitin(?:$|[^a-z0-9])"), "ubiquitin"),
]


def normalize_sample_id(value: str) -> str:
    return str(value).strip()


def normalize_chrom(value: str) -> str:
    text = str(value).strip()
    text = text.replace("chr", "").replace("Chr", "").replace("CHR", "")
    return text.lstrip("0") or "0"


def is_standard_chromosome(value: str) -> bool:
    chrom = normalize_chrom(value)
    try:
        return 1 <= int(chrom) <= 12
    except ValueError:
        return False


def is_transgene_element_name(value: str) -> bool:
    return canonical_transgene_element_name(value) is not None


def canonical_transgene_element_name(*args, **kwargs):
    from .transgene_bridge import engine
    return engine().canonical_transgene_element_name(*args, **kwargs)


def matches_transgene_record(*args, **kwargs):
    from .transgene_bridge import engine
    return engine().matches_transgene_record(*args, **kwargs)


def parse_info_field(info: str) -> Dict[str, str]:
    result: Dict[str, str] = {}
    text = str(info).strip()
    if not text or text == ".":
        return result
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
        else:
            result[item] = "True"
    return result


def read_vcf_samples(vcf_path: str | Path) -> List[str]:
    vcf_path = Path(vcf_path)
    with open_vcf(vcf_path) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                parts = line.rstrip("\n").split("\t")
                return parts[9:]
    raise ValueError(f"No #CHROM header found in {vcf_path}")


def vcf_header_and_records(vcf_path: str | Path):
    vcf_path = Path(vcf_path)
    with open_vcf(vcf_path) as fh:
        header = []
        for line in fh:
            if line.startswith("#"):
                header.append(line)
            else:
                yield header, line, fh
                return
    raise ValueError(f"No variant records found in {vcf_path}")


def vcf_to_variant_table(vcf_path: str | Path) -> pd.DataFrame:
    vcf_path = Path(vcf_path)
    samples = []
    records = []
    with open_vcf(vcf_path) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            chrom, pos, vid, ref, alt = parts[:5]
            locus_id = vid if vid != "." else f"{normalize_chrom(chrom)}_{pos}"
            sample_fields = parts[9:]
            genotypes = []
            for field in sample_fields:
                gt = field.split(":", 1)[0]
                if gt in {"", ".", "./.", ".|."}:
                    genotypes.append("./.")
                else:
                    genotypes.append(gt.replace("|", "/"))
            row = {
                "CHROM": chrom,
                "POS": pos,
                "LocusID": locus_id,
                "REF": ref,
                "ALT": alt,
            }
            row.update({sample: gt for sample, gt in zip(samples, genotypes)})
            records.append(row)
    return pd.DataFrame.from_records(records)


def _parse_format_map(sample_field: str, format_keys: Sequence[str]) -> Dict[str, str]:
    values = sample_field.split(":")
    return {key: value for key, value in zip(format_keys, values)}


def vcf_to_depth_table(vcf_path: str | Path) -> pd.DataFrame:
    vcf_path = Path(vcf_path)
    samples = []
    records = []
    with open_vcf(vcf_path) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            chrom, pos, vid, ref, alt, qual, flt, info = parts[:8]
            format_keys = parts[8].split(":")
            locus_id = vid if vid != "." else f"{normalize_chrom(chrom)}_{pos}"
            sample_fields = parts[9:]
            info_map = parse_info_field(info)
            sample_metrics = []
            for sample_name, sample_field in zip(samples, sample_fields):
                fm = _parse_format_map(sample_field, format_keys)
                gt = fm.get("GT", "./.")
                if gt in {"", ".", "./.", ".|."}:
                    gt = "./."
                else:
                    gt = gt.replace("|", "/")
                dp = fm.get("DP")
                ad = fm.get("AD")
                sample_metrics.append(
                    {
                        "SampleID": sample_name,
                        "Genotype": gt,
                        "DP": pd.to_numeric(dp, errors="coerce"),
                        "AD": ad,
                    }
                )
            records.append(
                {
                    "CHROM": chrom,
                    "POS": pos,
                    "LocusID": locus_id,
                    "REF": ref,
                    "ALT": alt,
                    "QUAL": qual,
                    "FILTER": flt,
                    "InfoDP": pd.to_numeric(info_map.get("DP"), errors="coerce"),
                    "sample_metrics": sample_metrics,
                }
            )
    return pd.DataFrame.from_records(records)


def summarize_variant_depths(vcf_path: str | Path) -> pd.DataFrame:
    depth_table = vcf_to_depth_table(vcf_path)
    records = []
    for _, row in depth_table.iterrows():
        sample_metrics = row["sample_metrics"]
        dps = [m["DP"] for m in sample_metrics if pd.notna(m["DP"])]
        gt_counts = {"REF": 0, "HET": 0, "ALT": 0, "MISSING": 0}
        alt_depth_total = 0.0
        ref_depth_total = 0.0
        for m in sample_metrics:
            gt = m["Genotype"]
            if gt in {"0/0", "0|0"}:
                gt_counts["REF"] += 1
            elif gt in {"0/1", "1/0", "0|1", "1|0"}:
                gt_counts["HET"] += 1
            elif gt in {"1/1", "1|1"}:
                gt_counts["ALT"] += 1
            else:
                gt_counts["MISSING"] += 1

            ad = m["AD"]
            if isinstance(ad, str) and ad not in {"", "."} and "," in ad:
                try:
                    ref_d, alt_d = ad.split(",", 1)
                    ref_depth_total += float(ref_d)
                    alt_depth_total += float(alt_d)
                except Exception:
                    pass

        dp_mean = float(pd.Series(dps).mean()) if dps else float("nan")
        dp_median = float(pd.Series(dps).median()) if dps else float("nan")
        dp_max = float(pd.Series(dps).max()) if dps else float("nan")
        dp_std = float(pd.Series(dps).std(ddof=0)) if len(dps) > 1 else float("nan")
        total_nonmissing = gt_counts["REF"] + gt_counts["HET"] + gt_counts["ALT"]
        alt_homo_fraction = gt_counts["ALT"] / total_nonmissing if total_nonmissing else float("nan")
        depth_ratio = dp_mean / dp_median if dp_median and pd.notna(dp_median) and dp_median != 0 else float("nan")
        records.append(
            {
                "CHROM": row["CHROM"],
                "POS": row["POS"],
                "LocusID": row["LocusID"],
                "REF": row["REF"],
                "ALT": row["ALT"],
                "InfoDP": row.get("InfoDP"),
                "SampleCount": len(sample_metrics),
                "MissingCount": gt_counts["MISSING"],
                "RefCount": gt_counts["REF"],
                "HetCount": gt_counts["HET"],
                "AltCount": gt_counts["ALT"],
                "AltHomoFraction": alt_homo_fraction,
                "DPMean": dp_mean,
                "DPMedian": dp_median,
                "DPMax": dp_max,
                "DPStd": dp_std,
                "DepthRatio": depth_ratio,
                "RefDepthTotal": ref_depth_total,
                "AltDepthTotal": alt_depth_total,
                "AltDepthFraction": alt_depth_total / (ref_depth_total + alt_depth_total) if (ref_depth_total + alt_depth_total) else float("nan"),
            }
        )
    return pd.DataFrame.from_records(records)


def summarize_transgene_element_depths(vcf_path: str | Path) -> pd.DataFrame:
    depth_table = vcf_to_depth_table(vcf_path)
    records = []
    for _, row in depth_table.iterrows():
        chrom = str(row["CHROM"])
        locus_id = str(row["LocusID"])
        if not matches_transgene_record(chrom, locus_id):
            continue
        sample_metrics = row["sample_metrics"]
        dps = [float(m["DP"]) for m in sample_metrics if pd.notna(m["DP"])]
        positive_dps = [dp for dp in dps if dp > 0]
        sample_dp_max = float(pd.Series(positive_dps).max()) if positive_dps else float("nan")
        records.append(
            {
                "CHROM": row["CHROM"],
                "POS": row["POS"],
                "LocusID": row["LocusID"],
                "REF": row["REF"],
                "ALT": row["ALT"],
                "InfoDP": pd.to_numeric(row.get("InfoDP"), errors="coerce"),
                "SampleCount": len(sample_metrics),
                "DPNonZeroCount": len(positive_dps),
                "DPMean": float(pd.Series(positive_dps).mean()) if positive_dps else float("nan"),
                "DPMedian": float(pd.Series(positive_dps).median()) if positive_dps else float("nan"),
                "DPMax": sample_dp_max,
                "DPGe10Count": int(sum(dp >= 10 for dp in positive_dps)),
                "DPGe20Count": int(sum(dp >= 20 for dp in positive_dps)),
                "DPGe50Count": int(sum(dp >= 50 for dp in positive_dps)),
                "DPGe100Count": int(sum(dp >= 100 for dp in positive_dps)),
            }
        )
    return pd.DataFrame.from_records(records)


def is_called_genotype(genotype):
    return str(genotype).strip() in {"0/0","0|0","0/1","1/0","0|1","1|0","1/1","1|1"}


def is_missing_genotype(genotype: str) -> bool:
    return not is_called_genotype(genotype)


def classify_transgene_dp(*args, **kwargs):
    from .transgene_bridge import engine
    return engine().classify_transgene_dp(*args, **kwargs)


def summarize_transgene_element_events(*args, **kwargs):
    from .transgene_bridge import engine
    return engine().summarize_transgene_element_events(*args, **kwargs)


def summarize_transgene_sample_events(*args, **kwargs):
    from .transgene_bridge import engine
    return engine().summarize_transgene_sample_events(*args, **kwargs)


def subset_vcf_variants(vcf_path: str | Path, variant_ids: Sequence[str]) -> pd.DataFrame:
    table = vcf_to_variant_table(vcf_path)
    ids = {str(x) for x in variant_ids}
    return table.loc[table["LocusID"].astype(str).isin(ids)].copy()
