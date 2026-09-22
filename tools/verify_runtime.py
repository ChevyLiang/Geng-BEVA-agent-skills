from __future__ import annotations

import importlib
import importlib.metadata
import json
import platform
import sys
from pathlib import Path

REQUIRED = {
    "numpy": "numpy",
    "pandas": "pandas",
    "matplotlib": "matplotlib",
    "plotly": "plotly",
    "scipy": "scipy",
    "scikit-learn": "sklearn",
    "seaborn": "seaborn",
}


def expected_versions(requirements: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for raw in requirements.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if "==" not in line:
            raise RuntimeError(f"Runtime requirement must be exactly pinned: {line}")
        name, version = line.split("==", 1)
        result[name.strip()] = version.strip()
    return result


def verify(requirements: Path) -> dict:
    expected = expected_versions(requirements)
    packages = {}
    errors = []
    for dist, import_name in REQUIRED.items():
        wanted = expected.get(dist)
        try:
            module = importlib.import_module(import_name)
            actual = importlib.metadata.version(dist)
            packages[dist] = {
                "import_name": import_name,
                "expected": wanted,
                "actual": actual,
                "ok": wanted == actual,
            }
            if wanted != actual:
                errors.append(f"{dist}: expected {wanted}, found {actual}")
        except Exception as exc:
            packages[dist] = {
                "import_name": import_name,
                "expected": wanted,
                "actual": None,
                "ok": False,
                "error": repr(exc),
            }
            errors.append(f"{dist}: import failed: {exc}")
    return {
        "ok": not errors,
        "python": sys.version.split()[0],
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "packages": packages,
        "errors": errors,
    }


def main() -> int:
    codex = Path(__file__).resolve().parents[1]
    req = codex / "runtime-requirements.txt"
    result = verify(req)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["ok"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
