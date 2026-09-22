from __future__ import annotations

import subprocess
import re
from _runtime import ensure_out_dir
from _plink import resolve_plink, run_plink
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd

ARROW_UP = "increase"
ARROW_DOWN = "decrease"
NEUTRAL = "neutral"
UNKNOWN = "unknown"


def project_root() -> Path:
    return Path(__file__).resolve().parents[1]


def workspace_root():
    return Path.cwd()


def default_qtn_path() -> Path:
    return project_root() / "assets" / "RiceNavi_subset.csv"


def default_plink_cmd():
    return resolve_plink(project_root())


def ensure_out_dir(out_dir):
    from _runtime import ensure_out_dir as resolve_output
    return resolve_output(out_dir)


def normalize_text(value: Any) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def is_missing_allele(value: Any) -> bool:
    return normalize_text(value) in {"", "."}


def canonical_locus_id(value: Any) -> str:
    text = normalize_text(value)
    if not text:
        return ""
    parts = text.split("_", 1)
    chrom = parts[0]
    rest = parts[1] if len(parts) > 1 else ""
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    chrom = chrom.lstrip("0") or "0"
    return f"{chrom}_{rest}" if rest else chrom


def normalize_list(value: Any) -> List[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        values = value
    else:
        values = [value]
    result: List[str] = []
    for item in values:
        text = normalize_text(item)
        if text:
            result.append(text)
    return result


def call_to_state(call: str) -> str:
    call = normalize_text(call)
    if call in {"0/0", "0|0"}:
        return "REF"
    if call in {"1/1", "1|1"}:
        return "ALT"
    if call in {"0/1", "1/0", "0|1", "1|0"}:
        return "HET"
    return "MISSING"


def direction_from_annotation(alt_function, alt_function_token="", effect=""):
    # Directions describe the annotation phenotype, never desirability.
    text = normalize_text(alt_function).lower()
    tokens = {"increasing":"increase", "increase":"increase", "higher":"increase", "larger":"increase", "more":"increase",
              "decreasing":"decrease", "deceasing":"decrease", "decrease":"decrease", "lower":"decrease", "fewer":"decrease", "less":"decrease", "shorter":"decrease",
              "gain":"gain", "loss":"loss", "delay":"delay", "delaying":"delay", "promote":"promote", "promoting":"promote", "preventing":"prevent", "prevent":"prevent",
              "neutral":"neutral", "unknown":"unknown", "uncertainty":"unknown", "ambiguous":"unknown"}
    # Text has priority. Mixed directions and uncertainty remain unresolved.
    if text:
        if any(x in text for x in ("uncertain", "unknown", "ambiguous", "not ")):
            return UNKNOWN, "alt_function_text"
        found = {tokens[w] for w in re.findall(r"[a-z]+", text) if w in tokens}
        if len(found) > 1: return UNKNOWN, "alt_function_text"
        if found: return found.pop(), "alt_function_text"
    token = normalize_text(alt_function_token).lower()
    if token in tokens: return tokens[token], "alt_function_token"
    # Positive/Negative are desirability labels; never turn these into trait directions.
    secondary = normalize_text(effect).lower()
    if not text and not token and secondary in tokens:
        return tokens[secondary], "effect"
    return UNKNOWN, "none"


def allele_alignment_status(
    vcf_ref: Any,
    vcf_alt: Any,
    qtn_ref: Any,
    qtn_alt: Any,
) -> str:
    vcf_ref_text = normalize_text(vcf_ref)
    vcf_alt_text = normalize_text(vcf_alt)
    qtn_ref_text = normalize_text(qtn_ref)
    qtn_alt_text = normalize_text(qtn_alt)

    if not qtn_ref_text and not qtn_alt_text:
        return "unavailable"
    if not vcf_ref_text and not vcf_alt_text:
        return "unavailable"
    if vcf_ref_text == qtn_ref_text and vcf_alt_text == qtn_alt_text:
        return "exact"
    if vcf_ref_text == qtn_alt_text and vcf_alt_text == qtn_ref_text:
        return "swapped"
    if vcf_ref_text == qtn_ref_text and is_missing_allele(vcf_alt_text):
        return "ref_only"
    if vcf_alt_text == qtn_alt_text and is_missing_allele(vcf_ref_text):
        return "alt_only"
    if vcf_ref_text == qtn_ref_text or vcf_alt_text == qtn_alt_text:
        return "partial_ref"
    if vcf_ref_text == qtn_alt_text or vcf_alt_text == qtn_ref_text:
        return "partial_alt"
    return "mismatched"


def genotype_state_from_call(call: Any, alignment_status: str) -> str:
    state = call_to_state(call)
    if state == "MISSING":
        return state
    if alignment_status == "swapped":
        if state == "REF":
            return "ALT"
        if state == "ALT":
            return "REF"
    if alignment_status in {"exact", "swapped"}:
        return state
    return "AMBIGUOUS"


def carries_favorable_allele(genotype_state: str, favorable_group: str) -> str:
    if genotype_state in {"AMBIGUOUS", "MISSING"} or favorable_group not in {"REF", "ALT"}:
        return ""
    if genotype_state == "HET":
        return "yes"
    return "yes" if genotype_state == favorable_group else "no"


def collect_functional_ids(qtn_path: str | Path) -> List[str]:
    qtn = pd.read_csv(qtn_path)
    if "ID" not in qtn.columns:
        raise ValueError(f"Functional QTN table missing ID column: {qtn_path}")
    ids = (
        qtn["ID"]
        .dropna()
        .astype(str)
        .map(str.strip)
        .tolist()
    )
    ids = [item for item in ids if item and item != "."]
    if not ids:
        raise ValueError(f"No functional locus IDs found in: {qtn_path}")
    return list(dict.fromkeys(ids))


def write_extract_list(ids: Iterable[str], path: Path) -> Path:
    path.write_text("\n".join(ids) + "\n", encoding="utf-8")
    return path


def is_vcf_like(path: str | Path) -> bool:
    suffixes = [item.lower() for item in Path(path).suffixes]
    if not suffixes:
        return False
    if suffixes[-1] == ".vcf":
        return True
    if suffixes[-1] == ".bcf":
        return False
    return ".vcf" in suffixes and suffixes[-1] == ".gz"


def run_plink_extract(
    vcf_path: str | Path,
    extract_list: str | Path,
    out_prefix: str | Path,
    *,
    plink_cmd: Optional[str] = None,
) -> Path:
    out_prefix = Path(out_prefix)
    cmd = [
        resolve_plink(project_root(), plink_cmd),
        "--vcf",
        str(vcf_path),
        "--chr-set", "12", "no-xy", "no-mt",
        "--double-id", "--vcf-half-call", "missing", "--keep-allele-order",
        "--extract",
        str(extract_list),
        "--recode",
        "vcf",
        "--out",
        str(out_prefix),
    ]
    subprocess.run(cmd, check=True)
    recoded_vcf = out_prefix.with_suffix(".vcf")
    if not recoded_vcf.exists():
        raise FileNotFoundError(f"plink did not produce expected VCF: {recoded_vcf}")
    return recoded_vcf


def parse_vcf_to_long_table(vcf_path: str | Path) -> pd.DataFrame:
    vcf_path = Path(vcf_path)
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF not found: {vcf_path}")

    samples: List[str] = []
    records: List[Dict[str, Any]] = []

    from _vcf import open_vcf
    with open_vcf(vcf_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                parts = line.split("\t")
                samples = parts[9:]
                continue

            parts = line.split("\t")
            if len(parts) < 10:
                continue

            chrom, pos, vid, ref, alt, *_rest = parts[:8]
            locus_id = vid if vid and vid != "." else f"{chrom}_{pos}"
            sample_fields = parts[9:]
            for sample_name, sample_field in zip(samples, sample_fields):
                fields = dict(zip(parts[8].split(":"), sample_field.split(":")))
                genotype = normalize_text(fields.get("GT", "./."))
                if genotype in {"", ".", "./.", ".|."}:
                    genotype = "./."
                else:
                    genotype = genotype.replace("|", "/")
                records.append(
                    {
                        "SampleID": sample_name,
                        "Gene": locus_id,
                        "LocusID": locus_id,
                        "Genotype": genotype,
                        "CHROM": chrom,
                        "POS": pos,
                        "REF": ref,
                        "ALT": alt,
                    }
                )

    if not records:
        raise ValueError("No genotype records found in VCF")
    return pd.DataFrame.from_records(records)


def _chrom_key(value: Any) -> str:
    text = normalize_text(value)
    if text.lower().startswith("chr"):
        text = text[3:]
    text = text.lstrip("0") or "0"
    return text


def _functional_site_keys(qtn_path: str | Path) -> tuple[set[str], set[tuple[str, int]]]:
    qtn = pd.read_csv(qtn_path)
    if "ID" not in qtn.columns:
        raise ValueError(f"Functional QTN table missing ID column: {qtn_path}")
    canonical_ids: set[str] = set()
    coord_keys: set[tuple[str, int]] = set()
    for raw in qtn["ID"].dropna().astype(str):
        cid = canonical_locus_id(raw)
        if not cid:
            continue
        canonical_ids.add(cid)
        parts = cid.split("_", 1)
        if len(parts) == 2:
            try:
                coord_keys.add((_chrom_key(parts[0]), int(parts[1])))
            except ValueError:
                pass
    return canonical_ids, coord_keys


def extract_functional_vcf_to_long_table(vcf_path: str | Path, qtn_path: str | Path) -> pd.DataFrame:
    """Stream a VCF and retain only RiceNavi-derived functional loci.

    This is the default backend because it is cross-platform and launches no
    external process, which makes it robust in Windows/Codex sandboxes and in
    workspaces whose paths contain non-ASCII characters.
    """
    from _vcf import open_vcf

    canonical_ids, coord_keys = _functional_site_keys(qtn_path)
    samples: List[str] = []
    records: List[Dict[str, Any]] = []

    with open_vcf(vcf_path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line or line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                parts = line.split("\t")
                samples = parts[9:]
                continue
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 10:
                continue
            chrom, pos_text, vid, ref, alt = parts[:5]
            chrom_key = _chrom_key(chrom)
            if chrom_key not in {str(i) for i in range(1, 13)}:
                continue
            try:
                pos = int(pos_text)
            except ValueError:
                continue
            locus_id = vid if vid and vid != "." else f"{chrom}_{pos}"
            cid = canonical_locus_id(locus_id)
            if cid not in canonical_ids and (chrom_key, pos) not in coord_keys:
                continue
            fmt = parts[8].split(":")
            for sample_name, sample_field in zip(samples, parts[9:]):
                fields = dict(zip(fmt, sample_field.split(":")))
                genotype = normalize_text(fields.get("GT", "./."))
                if genotype in {"", ".", "./.", ".|."}:
                    genotype = "./."
                else:
                    genotype = genotype.replace("|", "/")
                records.append({
                    "SampleID": sample_name,
                    "Gene": locus_id,
                    "LocusID": locus_id,
                    "Genotype": genotype,
                    "CHROM": chrom,
                    "POS": pos,
                    "REF": ref,
                    "ALT": alt,
                })

    if not records:
        raise ValueError("No bundled functional loci were found in the VCF")
    return pd.DataFrame.from_records(records)


def maybe_prepare_functional_table(
    input_path: str | Path,
    out_dir: str | Path,
    *,
    functional_qtn_csv: str | Path | None = None,
    plink_cmd: Optional[str] = None,
) -> Tuple[pd.DataFrame, Dict[str, str]]:
    out_path = ensure_out_dir(out_dir)
    qtn_path = Path(functional_qtn_csv or default_qtn_path())
    details: Dict[str, str] = {"functional_qtn_csv": str(qtn_path)}
    in_path = Path(input_path)

    if is_vcf_like(in_path):
        long_df = extract_functional_vcf_to_long_table(in_path, qtn_path)
        long_csv = out_path / "functional_long_table.csv"
        long_df.to_csv(long_csv, index=False, encoding="utf-8-sig")
        details.update({
            "extraction_backend": "python_stream",
            "functional_long_csv": str(long_csv),
        })
        if plink_cmd:
            details["warning"] = "--plink-cmd is retained for compatibility but is not required by the default Python extraction backend."
        return long_df, details

    if in_path.suffix.lower() != ".csv":
        raise ValueError(
            "input_path must be a VCF-like file or a CSV functional long table"
        )

    long_df = pd.read_csv(in_path)
    required = {"SampleID", "Gene", "LocusID", "Genotype"}
    missing = required - set(long_df.columns)
    if missing:
        raise ValueError(f"Input CSV missing required columns: {sorted(missing)}")
    return long_df, details


def load_qtn_metadata(qtn_path: str | Path) -> pd.DataFrame:
    """Load the RiceNavi-derived table and normalize legacy column names.

    The bundled production table uses biological allele semantics rather than
    VCF REF/ALT semantics:
      - Nipponbare_Allele: the Nipponbare/reference allele in RiceNavi
      - Alternative_Allele: the non-Nipponbare alternative allele
      - Alternative_Allele_Function/Direction/Effect: properties of that
        Alternative_Allele

    Legacy Unfunctional/Functional/ALT_FUNCTION/alt_function/Effect files are
    accepted for backward compatibility, but are normalized immediately.
    """
    qtn = pd.read_csv(qtn_path, dtype=str, keep_default_na=False)
    legacy = {
        "Effect": "Alternative_Allele_Effect",
        "Unfunctional": "Nipponbare_Allele",
        "Functional": "Alternative_Allele",
        "ALT_FUNCTION": "Alternative_Allele_Function",
        "alt_function": "Alternative_Allele_Direction",
    }
    for old, new_name in legacy.items():
        if new_name not in qtn.columns and old in qtn.columns:
            qtn = qtn.rename(columns={old: new_name})
    required = {
        "ID",
        "Alternative_Allele_Effect",
        "Alternative_Allele_Function",
        "Alternative_Allele_Direction",
        "Nipponbare_Allele",
        "Alternative_Allele",
        "GENEname",
        "Category",
        "Trait",
    }
    missing = required - set(qtn.columns)
    if missing:
        raise ValueError(f"Functional QTN table missing required columns: {sorted(missing)}")
    qtn = qtn.copy()
    qtn["_canonical_id"] = qtn["ID"].astype(str).map(canonical_locus_id)
    qtn["_gene_key"] = qtn["GENEname"].astype(str).str.casefold()
    return qtn


def _is_heading_trait(meta: Dict[str, Any]) -> bool:
    text = " ".join(
        normalize_text(meta.get(key, ""))
        for key in ("Trait", "Category", "Alternative_Allele_Function")
    ).lower()
    return "heading" in text or "flowering" in text


def _default_favorable_rule(meta: Dict[str, Any], direction: str) -> Dict[str, str]:
    """Return the project-default favorable allele in RiceNavi coordinates.

    Project default rules are binary whenever RiceNavi provides a resolvable
    direction:
    1. Heading date is optimized for earlier heading. Alternative promote ->
       Alternative favorable / Nipponbare unfavorable; Alternative delay ->
       Nipponbare favorable / Alternative unfavorable.
    2. Otherwise, Alternative_Allele_Effect=Positive -> Alternative favorable
       / Nipponbare unfavorable.
    3. Alternative_Allele_Effect=Negative -> Nipponbare favorable /
       Alternative unfavorable.
    4. Non-heading Uncertainty remains unresolved rather than being invented.
    """
    result = {
        "BreedingDirection": "unspecified",
        "PreferenceSource": "none",
        "FavorableAllele": "",
        "FavorableAlleleGroup": "UNKNOWN",
        "FavorableRule": "unresolved",
    }
    nip = normalize_text(meta.get("NipponbareAllele", ""))
    alt = normalize_text(meta.get("AlternativeAllele", ""))
    function_text = normalize_text(meta.get("Alternative_Allele_Function", "")).lower()

    if _is_heading_trait(meta):
        if direction == "promote" or any(token in function_text for token in ("early heading", "earlier heading", "advance heading", "advancing heading", "promoting heading")):
            result.update(BreedingDirection="early_heading", PreferenceSource="default:early_heading",
                          FavorableAllele=alt, FavorableAlleleGroup="ALTERNATIVE", FavorableRule="early_heading_alternative")
            return result
        if direction == "delay" or any(token in function_text for token in ("delay heading", "delaying heading", "late heading", "later heading")):
            result.update(BreedingDirection="early_heading", PreferenceSource="default:early_heading",
                          FavorableAllele=nip, FavorableAlleleGroup="NIPPONBARE", FavorableRule="early_heading_nipponbare")
            return result

    effect = normalize_text(meta.get("Alternative_Allele_Effect", "")).casefold()
    if effect == "positive":
        result.update(BreedingDirection="database_positive", PreferenceSource="default:RiceNavi_positive",
                      FavorableAllele=alt, FavorableAlleleGroup="ALTERNATIVE", FavorableRule="RiceNavi_positive")
    elif effect == "negative":
        result.update(BreedingDirection="database_negative", PreferenceSource="default:RiceNavi_negative",
                      FavorableAllele=nip, FavorableAlleleGroup="NIPPONBARE", FavorableRule="RiceNavi_negative")
    return result

def annotate_long_table(input_df: pd.DataFrame, qtn_path: str | Path, breeding_preferences=None, environment=None) -> pd.DataFrame:
    qtn = load_qtn_metadata(qtn_path)
    preferences = load_preferences(breeding_preferences)
    records: List[Dict[str, Any]] = []

    # Preserve multiple gene/trait annotations at the same locus.
    expanded = []
    for _, original in input_df.iterrows():
        matches = qtn.loc[qtn["_canonical_id"] == canonical_locus_id(original["LocusID"])]
        if matches.empty:
            expanded.append((original, None))
        else:
            expanded.extend((original, annotation) for _, annotation in matches.iterrows())

    for row, annotation in expanded:
        locus_id = normalize_text(row["LocusID"])
        gene = normalize_text(row["Gene"])

        if annotation is None:
            meta = {
                "Trait": normalize_text(row.get("Trait", "")),
                "Category": normalize_text(row.get("Category", "")),
                "Alternative_Allele_Effect": "",
                "Alternative_Allele_Function": "",
                "Alternative_Allele_Direction": "",
                "NipponbareAllele": "",
                "AlternativeAllele": "",
                "ReferenceGene": gene,
            }
        else:
            ref = annotation
            meta = {
                "Trait": normalize_text(ref["Trait"]),
                "Category": normalize_text(ref["Category"]),
                "Alternative_Allele_Effect": normalize_text(ref["Alternative_Allele_Effect"]),
                "Alternative_Allele_Function": normalize_text(ref["Alternative_Allele_Function"]),
                "Alternative_Allele_Direction": normalize_text(ref["Alternative_Allele_Direction"]),
                "NipponbareAllele": normalize_text(ref["Nipponbare_Allele"]),
                "AlternativeAllele": normalize_text(ref["Alternative_Allele"]),
                "ReferenceGene": normalize_text(ref["GENEname"]),
            }

        genotype = normalize_text(row["Genotype"])
        vcf_ref = normalize_text(row.get("REF", ""))
        vcf_alt = normalize_text(row.get("ALT", ""))

        # Crucial correction: VCF REF/ALT may be frequency-oriented and must
        # never be interpreted as Nipponbare/non-Nipponbare without alignment.
        alignment_status = allele_alignment_status(
            vcf_ref, vcf_alt, meta.get("NipponbareAllele", ""), meta.get("AlternativeAllele", "")
        )
        annotation_state = genotype_state_from_call(genotype, alignment_status)
        rice_navi_state = {"REF": "NIPPONBARE_HOM", "ALT": "ALTERNATIVE_HOM", "HET": "HET",
                           "MISSING": "MISSING", "AMBIGUOUS": "AMBIGUOUS"}.get(annotation_state, "AMBIGUOUS")
        alt_dosage = {"NIPPONBARE_HOM": 0, "HET": 1, "ALTERNATIVE_HOM": 2}.get(rice_navi_state, "")

        direction, direction_source = direction_from_annotation(
            meta.get("Alternative_Allele_Function", ""),
            meta.get("Alternative_Allele_Direction", ""),
            meta.get("Alternative_Allele_Effect", ""),
        )
        breeding = evaluate_preference(locus_id, meta, direction, alignment_status, preferences, environment)
        favorable_group = breeding["FavorableAlleleGroup"]

        # Convert biological favorable group into the aligned annotation state
        # used by carrier logic, while keeping the explicit RiceNavi label.
        carrier_group = {"NIPPONBARE": "REF", "ALTERNATIVE": "ALT"}.get(favorable_group, "UNKNOWN")
        carrier = carries_favorable_allele(annotation_state, carrier_group)
        if carrier_group == "UNKNOWN" or annotation_state in {"MISSING", "AMBIGUOUS"}:
            favorable_genotype_status = "unresolved"
        elif annotation_state == "HET":
            favorable_genotype_status = "heterozygous"
        elif annotation_state == carrier_group:
            favorable_genotype_status = "homozygous_favorable"
        else:
            favorable_genotype_status = "homozygous_unfavorable"

        unfavorable_group = {"NIPPONBARE": "ALTERNATIVE", "ALTERNATIVE": "NIPPONBARE"}.get(favorable_group, "UNKNOWN")
        unfavorable_allele = meta.get("AlternativeAllele", "") if unfavorable_group == "ALTERNATIVE" else (meta.get("NipponbareAllele", "") if unfavorable_group == "NIPPONBARE" else "")
        if rice_navi_state == "HET":
            advantage_class = "HETEROZYGOUS"
        elif favorable_group == "UNKNOWN" or rice_navi_state in {"MISSING", "AMBIGUOUS"}:
            advantage_class = "UNRESOLVED"
        elif rice_navi_state == "NIPPONBARE_HOM":
            advantage_class = "NIPPONBARE_FAVORABLE" if favorable_group == "NIPPONBARE" else "NIPPONBARE_UNFAVORABLE"
        elif rice_navi_state == "ALTERNATIVE_HOM":
            advantage_class = "ALTERNATIVE_FAVORABLE" if favorable_group == "ALTERNATIVE" else "ALTERNATIVE_UNFAVORABLE"
        else:
            advantage_class = "UNRESOLVED"

        records.append({
            "SampleID": normalize_text(row["SampleID"]),
            "Gene": normalize_text(meta["ReferenceGene"] or row.get("Gene", "")),
            "Trait": meta["Trait"],
            "LocusID": locus_id,
            "Category": meta["Category"],
            "VCFGenotype": genotype,
            "VCFRefAllele": vcf_ref,
            "VCFAltAllele": vcf_alt,
            "AlignmentStatus": alignment_status,
            "NipponbareAllele": meta["NipponbareAllele"],
            "AlternativeAllele": meta["AlternativeAllele"],
            "RiceNaviGenotypeState": rice_navi_state,
            "AlternativeAlleleDosage": alt_dosage,
            "Alternative_Allele_Effect": meta["Alternative_Allele_Effect"],
            "Alternative_Allele_Function": meta["Alternative_Allele_Function"],
            "Alternative_Allele_Direction": meta["Alternative_Allele_Direction"],
            "EffectDirection": direction,
            "DirectionSource": direction_source,
            **breeding,
            "UnfavorableAllele": unfavorable_allele,
            "UnfavorableAlleleGroup": unfavorable_group,
            "CarriesFavorableAllele": carrier,
            "FavorableGenotypeStatus": favorable_genotype_status,
            "AdvantageClass": advantage_class,
            "GenotypeNote": (
                "VCF alleles were aligned to RiceNavi before biological interpretation; heterozygotes are reported separately and are not forced into favorable/unfavorable homozygous classes."
                if favorable_genotype_status == "heterozygous"
                else "VCF alleles were aligned to RiceNavi before biological interpretation."
            ),
            "CHROM": normalize_text(row.get("CHROM", "")),
            "POS": normalize_text(row.get("POS", "")),
        })

    result = pd.DataFrame.from_records(records)
    return result.reindex(columns=[
        "SampleID", "Gene", "Trait", "LocusID", "Category",
        "VCFGenotype", "VCFRefAllele", "VCFAltAllele", "AlignmentStatus",
        "NipponbareAllele", "AlternativeAllele", "RiceNaviGenotypeState", "AlternativeAlleleDosage",
        "Alternative_Allele_Effect", "Alternative_Allele_Function", "Alternative_Allele_Direction",
        "EffectDirection", "DirectionSource",
        "BreedingDirection", "BreedingContext", "PreferenceSource", "FavorableRule",
        "FavorableAllele", "FavorableVCFAlleleGroup", "FavorableAlleleGroup",
        "UnfavorableAllele", "UnfavorableAlleleGroup",
        "CarriesFavorableAllele", "FavorableGenotypeStatus", "AdvantageClass", "GenotypeNote",
        "CHROM", "POS",
    ])

def filter_annotated_table(
    df: pd.DataFrame,
    *,
    sample_ids: Any = None,
    traits: Any = None,
    categories: Any = None,
    query_logic: str = "and",
) -> pd.DataFrame:
    masks: List[pd.Series] = []
    if sample_ids is not None:
        masks.append(df["SampleID"].isin(normalize_list(sample_ids)))
    if traits is not None:
        masks.append(df["Trait"].isin(normalize_list(traits)))
    if categories is not None:
        masks.append(df["Category"].isin(normalize_list(categories)))

    if not masks:
        return df

    query_logic = normalize_text(query_logic).lower()
    if query_logic == "or":
        mask = masks[0].copy()
        for item in masks[1:]:
            mask |= item
    else:
        mask = masks[0].copy()
        for item in masks[1:]:
            mask &= item
    return df.loc[mask].copy()


def load_preferences(path):
    if path is None: return pd.DataFrame(columns=["Trait", "BreedingDirection", "LocusID", "FavorableAllele", "Environment"])
    pref = pd.read_csv(path, dtype=str, keep_default_na=False)
    if "Trait" not in pref: raise ValueError("Preference CSV requires Trait")
    for c in ["BreedingDirection", "LocusID", "FavorableAllele", "Environment"]:
        if c not in pref: pref[c] = ""
    allowed = {"", "increase", "decrease", "unspecified"}
    if not set(pref.BreedingDirection).issubset(allowed):
        raise ValueError("BreedingDirection must be increase, decrease, or unspecified")
    if ((pref.FavorableAllele != "") & (pref.LocusID == "")).any():
        raise ValueError("An explicit FavorableAllele requires LocusID")
    pref["_source"] = str(Path(path).resolve())
    return pref


def _trait_key(value):
    return re.sub(r"[^a-z0-9]+", " ", normalize_text(value).lower()).strip()


def evaluate_preference(locus_id, meta, direction, alignment, preferences, environment):
    default = _default_favorable_rule(meta, direction)
    result = {
        "BreedingDirection": default["BreedingDirection"],
        "BreedingContext": normalize_text(environment),
        "PreferenceSource": default["PreferenceSource"],
        "FavorableAllele": default["FavorableAllele"],
        "FavorableAlleleGroup": default["FavorableAlleleGroup"],
        "FavorableVCFAlleleGroup": "UNKNOWN",
        "FavorableRule": default["FavorableRule"],
    }

    # Explicit user preferences override all project defaults.
    candidates = pd.DataFrame()
    if not preferences.empty:
        candidates = preferences.loc[preferences.Trait.map(_trait_key) == _trait_key(meta["Trait"])]
        candidates = candidates.loc[(candidates.Environment == "") | (candidates.Environment == normalize_text(environment))]
        specific = candidates.loc[candidates.LocusID.map(canonical_locus_id) == canonical_locus_id(locus_id)]
        candidates = specific if not specific.empty else candidates.loc[candidates.LocusID == ""]
        env_specific = candidates.loc[(candidates.Environment != "") & (candidates.Environment == normalize_text(environment))]
        if not env_specific.empty:
            candidates = env_specific
        if len(candidates) > 1:
            raise ValueError(f"Multiple matching breeding preferences: {locus_id}, {meta['Trait']}")

    if not candidates.empty:
        pref = candidates.iloc[0]
        desired = pref.BreedingDirection or "unspecified"
        result.update(BreedingDirection=desired, PreferenceSource=pref._source, FavorableRule="explicit_preference")
        allele = normalize_text(pref.FavorableAllele)
        nip, alt = meta["NipponbareAllele"], meta["AlternativeAllele"]
        if allele:
            if allele not in {nip, alt}:
                raise ValueError(f"Preference allele does not match RiceNavi annotation at {locus_id}")
            group = "ALTERNATIVE" if allele == alt else "NIPPONBARE"
            result.update(FavorableAllele=allele, FavorableAlleleGroup=group)
        else:
            # Numeric trait override remains supported. If it cannot be resolved,
            # the explicit preference intentionally suppresses project defaults.
            text = normalize_text(meta["Alternative_Allele_Function"]).lower()
            match = re.fullmatch(r"(?:increasing|decreasing|deceasing|higher|lower|larger|fewer|shorter|more|less)\s+(.+)", text)
            compatible = match is not None and _trait_key(match.group(1)) == _trait_key(meta["Trait"])
            if compatible and direction in {"increase", "decrease"} and desired in {"increase", "decrease"}:
                group = "ALTERNATIVE" if direction == desired else "NIPPONBARE"
                allele = alt if group == "ALTERNATIVE" else nip
                result.update(FavorableAllele=allele, FavorableAlleleGroup=group)
            else:
                result.update(FavorableAllele="", FavorableAlleleGroup="UNKNOWN")

    if alignment in {"exact", "swapped"} and result["FavorableAlleleGroup"] in {"NIPPONBARE", "ALTERNATIVE"}:
        annotation_group = "ALT" if result["FavorableAlleleGroup"] == "ALTERNATIVE" else "REF"
        vcf_group = annotation_group if alignment == "exact" else ("REF" if annotation_group == "ALT" else "ALT")
        result["FavorableVCFAlleleGroup"] = vcf_group
    return result
