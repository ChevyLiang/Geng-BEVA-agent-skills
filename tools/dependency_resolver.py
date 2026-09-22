from __future__ import annotations

"""Dependency verifier/repair bridge for Geng-BEVA.

Windows uses the fixed workspace virtual environment prepared by
``bootstrap_runtime.ps1``. Missing packages are repaired by synchronizing the
entire pinned runtime, never by silently omitting an analysis/visualization.
PLINK is not managed here.

A legacy package-level fallback remains only for non-Windows developer/test
execution, where the fixed Windows bootstrap is not applicable.
"""

import importlib
import os
from pathlib import Path
import subprocess
import sys
from typing import Iterable, Mapping

PACKAGE_IMPORT_NAMES: Mapping[str, str] = {
    "scikit-learn": "sklearn",
    "sklearn": "sklearn",
    "pyyaml": "yaml",
    "pillow": "PIL",
}


def codex_root(start: str | Path | None = None) -> Path:
    candidates = []
    if start is not None:
        p = Path(start).resolve()
        candidates.extend([p] + list(p.parents))
    candidates.extend([Path.cwd().resolve()] + list(Path.cwd().resolve().parents))
    here = Path(__file__).resolve()
    candidates.extend([here] + list(here.parents))
    for p in candidates:
        if p.name == ".codex" and (p / "skills").exists():
            return p
        if (p / ".codex" / "skills").exists():
            return p / ".codex"
    raise RuntimeError("Could not locate the .codex workspace root.")


def _import_name(package: str, import_name: str | None = None) -> str:
    if import_name:
        return import_name
    key = package.split("==", 1)[0].split(">=", 1)[0].split("<=", 1)[0].strip().lower()
    return PACKAGE_IMPORT_NAMES.get(key, key.replace("-", "_"))


def _can_import(name: str) -> bool:
    importlib.invalidate_caches()
    try:
        importlib.import_module(name)
        return True
    except Exception:
        return False


def _repair_windows_runtime(start: str | Path | None = None) -> None:
    root = codex_root(start)
    bootstrap = root / "tools" / "bootstrap_runtime.ps1"
    cmd = [
        "powershell.exe", "-NoProfile", "-ExecutionPolicy", "Bypass",
        "-File", str(bootstrap), "-ForceRepair", "-Quiet",
    ]
    completed = subprocess.run(cmd, check=False)
    if completed.returncode != 0:
        raise RuntimeError(
            "Geng-BEVA fixed runtime repair failed. Network, sandbox, or filesystem policy may be blocking it."
        )


def ensure_python_package(package: str, *, import_name: str | None = None, start: str | Path | None = None) -> str:
    name = _import_name(package, import_name)
    if _can_import(name):
        return "runtime" if os.environ.get("GENGBEVA_FIXED_RUNTIME") == "1" else "environment"

    if os.name == "nt":
        _repair_windows_runtime(start)
        if _can_import(name):
            return "repaired-runtime"
        raise RuntimeError(
            f"Required dependency {package!r} is still unavailable after fixed-runtime repair. "
            "The analysis is incomplete and must not be reported as successfully finished."
        )

    # Non-Windows developer/test fallback only.
    if os.environ.get("GENGBEVA_NO_AUTO_INSTALL", "").strip().lower() in {"1", "true", "yes"}:
        raise RuntimeError(f"Missing dependency {package!r}; auto-install is disabled.")
    cmd = [sys.executable, "-m", "pip", "install", package]
    subprocess.run(cmd, check=True)
    importlib.invalidate_caches()
    if not _can_import(name):
        raise RuntimeError(f"Dependency {package!r} remained unavailable after installation.")
    return "installed"


def ensure_python_packages(packages: Iterable[str], *, start: str | Path | None = None) -> dict[str, str]:
    return {pkg: ensure_python_package(pkg, start=start) for pkg in packages}


def bootstrap_from_requirements(requirements_file: str | Path) -> dict[str, str]:
    path = Path(requirements_file)
    if not path.exists():
        return {}
    packages = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line and not line.startswith("#"):
            packages.append(line)
    return ensure_python_packages(packages, start=path)
