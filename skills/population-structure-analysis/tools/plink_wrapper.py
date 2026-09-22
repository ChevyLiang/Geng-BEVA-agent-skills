from __future__ import annotations

import subprocess
import tempfile
from .plink_resolver import resolve_plink, run_plink
from core.runtime import ensure_out_dir
from pathlib import Path
from typing import Dict, Iterable, Optional, Sequence


def project_root() -> Path:
    return Path(__file__).resolve().parents[1]


def workspace_root():
    return Path.cwd()


def default_plink_cmd():
    return resolve_plink(project_root())


def ensure_out_dir(out_dir):
    from core.runtime import ensure_out_dir as resolve_output
    return resolve_output(out_dir)


def _run_plink(args, *, plink_cmd=None):
    run_plink(args, project_root(), plink_cmd)


def _write_keep_file(out_prefix: str | Path, keep_samples: Sequence[str]) -> str:
    out_prefix = Path(out_prefix)
    keep_path = out_prefix.with_suffix(".keep.txt")
    lines = [f"{sample}\t{sample}" for sample in keep_samples]
    keep_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(keep_path)


def run_plink_pca(
    vcf_path: str | Path,
    out_prefix: str | Path,
    *,
    plink_cmd: Optional[str] = None,
    keep_samples: Optional[Sequence[str]] = None,
) -> Dict[str, str]:
    out_prefix = Path(out_prefix)
    args = [
        "--vcf",
        str(vcf_path),
    ]
    keep_path = None
    if keep_samples:
        keep_path = _write_keep_file(out_prefix, keep_samples)
        args.extend(["--keep", keep_path])
    args.extend(
        [
            "--chr-set", "12", "no-xy", "no-mt",
            "--double-id", "--vcf-half-call", "missing",
            "--pca",
            "20",
            "--out",
            str(out_prefix),
        ]
    )
    _run_plink(
        args,
        plink_cmd=plink_cmd,
    )
    outputs = {
        "eigenvec": str(out_prefix.with_suffix(".eigenvec")),
        "eigenval": str(out_prefix.with_suffix(".eigenval")),
        "log": str(out_prefix.with_suffix(".log")),
    }
    if keep_path:
        outputs["keep"] = keep_path
    return outputs


def run_plink_ibs(
    vcf_path: str | Path,
    out_prefix: str | Path,
    *,
    plink_cmd: Optional[str] = None,
    keep_samples: Optional[Sequence[str]] = None,
) -> Dict[str, str]:
    out_prefix = Path(out_prefix)
    args = [
        "--vcf",
        str(vcf_path),
    ]
    keep_path = None
    if keep_samples:
        keep_path = _write_keep_file(out_prefix, keep_samples)
        args.extend(["--keep", keep_path])
    args.extend(
        [
            "--chr-set", "12", "no-xy", "no-mt",
            "--double-id", "--vcf-half-call", "missing",
            "--genome",
            "--out",
            str(out_prefix),
        ]
    )
    _run_plink(
        args,
        plink_cmd=plink_cmd,
    )
    outputs = {
        "genome": str(out_prefix.with_suffix(".genome")),
        "log": str(out_prefix.with_suffix(".log")),
    }
    if keep_path:
        outputs["keep"] = keep_path
    return outputs


def run_plink_ibd(
    vcf_path: str | Path,
    out_prefix: str | Path,
    *,
    plink_cmd: Optional[str] = None,
    keep_samples: Optional[Sequence[str]] = None,
) -> Dict[str, str]:
    out_prefix = Path(out_prefix)
    args = [
        "--vcf",
        str(vcf_path),
    ]
    keep_path = None
    if keep_samples:
        keep_path = _write_keep_file(out_prefix, keep_samples)
        args.extend(["--keep", keep_path])
    args.extend(
        [
            "--chr-set", "12", "no-xy", "no-mt",
            "--double-id", "--vcf-half-call", "missing",
            "--genome",
            "--out",
            str(out_prefix),
        ]
    )
    _run_plink(
        args,
        plink_cmd=plink_cmd,
    )
    outputs = {
        "genome": str(out_prefix.with_suffix(".genome")),
        "log": str(out_prefix.with_suffix(".log")),
    }
    if keep_path:
        outputs["keep"] = keep_path
    return outputs
