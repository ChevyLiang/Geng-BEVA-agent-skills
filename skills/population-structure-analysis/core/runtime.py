from __future__ import annotations
import functools
import json
import os
from pathlib import Path


def ensure_out_dir(out_dir="outputs"):
    """Relative paths are relative to the caller's working directory."""
    path = Path(out_dir).expanduser().resolve()
    path.mkdir(parents=True, exist_ok=True)
    return path


def build_manifest(module, input_path, out_dir, *, notes=None, extra=None, **outputs):
    return {"schema_version": "1.0", "status": "success", "module": module,
            "input": str(Path(input_path).resolve()), "out_dir": str(Path(out_dir).resolve()),
            "outputs": {k: str(Path(v).resolve()) for k,v in outputs.items() if v is not None},
            "warnings": [notes] if notes else [], "metrics": {}, "extra": extra or {}}


def _json_default(value):
    if isinstance(value, Path): return str(value)
    if hasattr(value, "item"): return value.item()
    raise TypeError(f"Not JSON serializable: {type(value).__name__}")


def save_manifest(result, out_path, module):
    payload = json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False, default=_json_default)
    for filename in (module + ".manifest.json", "manifest.json"):
        target = out_path / filename
        tmp = target.with_suffix(target.suffix + ".tmp")
        tmp.write_text(payload + "\n", encoding="utf-8")
        os.replace(tmp, target)


def contract(fn):
    @functools.wraps(fn)
    def wrapped(input_path, out_dir="outputs", **kwargs):
        out_path = ensure_out_dir(out_dir)
        module = fn.__name__
        result = build_manifest(module, input_path, out_path)
        try:
            source = Path(input_path).expanduser().resolve()
            if not source.is_file(): raise FileNotFoundError(f"Input file not found: {source}")
            analysis_source = source
            if module in {"pca_analysis","ibs_analysis","genotype_frequency","individual_diversity","sample_pairwise_compare","pairwise_matrix","fst_analysis"}:
                from .canonical_vcf import export_canonical_vcf
                analysis_source = export_canonical_vcf(source,out_path/"_canonical"/"input.vcf")
            raw = fn(analysis_source, out_dir=out_path, **kwargs)
            if not isinstance(raw, dict): raise TypeError("Analysis must return a manifest dictionary")
            result.update(raw)
            result.update(schema_version="1.0", module=module, input=str(source), out_dir=str(out_path))
            if result.get("status") in {None, "ok"}: result["status"] = "success"
            result.setdefault("warnings", []); result.setdefault("metrics", {})
            result["outputs"] = {k: str(Path(v).resolve()) for k,v in result.get("outputs", result.get("output_files", {})).items() if v is not None}
            for key, value in result["outputs"].items():
                if not Path(value).is_file(): raise FileNotFoundError(f"Missing declared output {key}: {value}")
            save_manifest(result, out_path, module)
            return result
        except Exception as exc:
            result.update(status="error", outputs={}, error={"type": type(exc).__name__, "message": str(exc)})
            save_manifest(result, out_path, module)
            raise
    return wrapped
