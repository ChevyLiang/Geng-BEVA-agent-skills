from __future__ import annotations

"""Windows fixed-runtime launcher for direct Python entry-point calls.

Canonical Codex execution uses ``.codex/run-geng-beva.ps1``. This guard exists so
an accidental direct ``python main.py ...`` on Windows repairs the fixed runtime
and re-runs the same command with the managed Geng-BEVA interpreter instead of
continuing in an arbitrary system/Conda Python.
"""

import os
from pathlib import Path
import subprocess
import sys


def _codex_root(start: str | Path) -> Path:
    p = Path(start).resolve()
    for parent in [p] + list(p.parents):
        if parent.name == ".codex":
            return parent
        if (parent / ".codex" / "skills").exists():
            return parent / ".codex"
    raise RuntimeError("Could not locate .codex root.")


def ensure_fixed_runtime(start: str | Path) -> None:
    if os.name != "nt":
        return
    if os.environ.get("GENGBEVA_FIXED_RUNTIME") == "1":
        return

    codex = _codex_root(start)
    runtime_python = codex / "runtime" / "geng-beva" / ".venv" / "Scripts" / "python.exe"
    bootstrap = codex / "tools" / "bootstrap_runtime.ps1"

    try:
        same = runtime_python.exists() and Path(sys.executable).resolve() == runtime_python.resolve()
    except Exception:
        same = False
    if same:
        os.environ["GENGBEVA_FIXED_RUNTIME"] = "1"
        return

    cmd = [
        "powershell.exe",
        "-NoProfile",
        "-ExecutionPolicy",
        "Bypass",
        "-File",
        str(bootstrap),
        "-Quiet",
    ]
    completed = subprocess.run(cmd, check=False)
    if completed.returncode != 0 or not runtime_python.exists():
        raise RuntimeError(
            "Geng-BEVA could not prepare its fixed Windows runtime. "
            "Runtime bootstrap/repair must succeed before analysis; do not silently skip required outputs."
        )

    env = os.environ.copy()
    env["GENGBEVA_FIXED_RUNTIME"] = "1"
    completed = subprocess.run([str(runtime_python), *sys.argv], env=env, check=False)
    raise SystemExit(completed.returncode)
