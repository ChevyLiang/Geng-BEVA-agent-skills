from __future__ import annotations
import gzip
import re
from pathlib import Path


def open_vcf(path):
    path = Path(path)
    return gzip.open(path, "rt", encoding="utf-8-sig") if path.suffix.lower() == ".gz" else path.open(encoding="utf-8-sig")


def canonical_chrom(chrom):
    text = re.sub(r"^chr", "", str(chrom).strip(), flags=re.I)
    return text.isdecimal() and 1 <= int(text) <= 12


def export_canonical_vcf(input_path, output_path):
    source, target = Path(input_path).resolve(), Path(output_path).resolve()
    if source == target: raise ValueError("Input and output VCF must be different files")
    target.parent.mkdir(parents=True, exist_ok=True)
    header_seen = False
    with open_vcf(source) as src, target.open("w", encoding="utf-8", newline="\n") as dst:
        for line in src:
            if line.startswith("#"):
                if line.startswith("#CHROM"): header_seen = True
                dst.write(line)
            elif line.strip():
                if not header_seen: raise ValueError("VCF records precede #CHROM header")
                if canonical_chrom(line.split("\t", 1)[0]): dst.write(line)
    if not header_seen: raise ValueError("VCF #CHROM header missing")
    return target
