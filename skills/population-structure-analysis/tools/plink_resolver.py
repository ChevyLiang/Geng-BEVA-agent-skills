from __future__ import annotations

import os
import platform
import re
import shutil
import subprocess
from pathlib import Path


def _validate_plink(path: Path) -> str:
    """Validate that *path* is a genomic PLINK 1.9 executable."""
    check = subprocess.run(
        [str(path), "--version"],
        capture_output=True,
        text=True,
        timeout=20,
        check=True,
    )
    version = check.stdout + check.stderr
    if not re.search(r"PLINK\s+v?1\.9", version, re.I):
        raise ValueError("Expected genomic PLINK 1.9 (not PuTTY plink or PLINK 2)")
    return version.strip().splitlines()[0] if version.strip() else "PLINK 1.9"


def _bundled_candidates(skill_root: str | Path | None) -> list[Path]:
    """Return repository-local PLINK locations without writing anywhere.

    The preferred Windows layout is::

        <skill_root>/tools/plink/windows-x64/plink.exe

    ``skill_root`` normally points to ``population-structure-analysis``.
    A few nearby legacy layouts are also accepted for compatibility.
    """
    if skill_root is None:
        root = Path(__file__).resolve().parents[1]
    else:
        root = Path(skill_root).expanduser().resolve()

    system = platform.system().lower()
    machine = platform.machine().lower()
    candidates: list[Path] = []

    if system == "windows":
        if machine in {"amd64", "x86_64"}:
            candidates.append(root / "tools" / "plink" / "windows-x64" / "plink.exe")
        else:
            candidates.append(root / "tools" / "plink" / "windows-x86" / "plink.exe")
        candidates.extend(
            [
                root / "tools" / "plink" / "plink.exe",
                root / "tools" / "plink.exe",
            ]
        )
    elif system == "linux":
        candidates.extend(
            [
                root / "tools" / "plink" / "linux-x64" / "plink",
                root / "tools" / "plink" / "plink",
            ]
        )
    elif system == "darwin":
        candidates.extend(
            [
                root / "tools" / "plink" / "macos" / "plink",
                root / "tools" / "plink" / "plink",
            ]
        )

    return candidates


def resolve_plink(skill_root=None, explicit=None, auto_install=False):
    """Resolve a usable PLINK 1.9 executable without modifying the system.

    Resolution order:
      1. explicit ``--plink-cmd``;
      2. repository-local bundled executable (Windows-first deployment);
      3. ``PLINK_BIN`` environment variable;
      4. ``plink`` / ``plink1.9`` on ``PATH``.

    ``auto_install`` is retained only for API compatibility with earlier
    versions. This resolver deliberately performs no automatic installation or
    writes to the user's home directory.
    """
    candidates: list[str | Path | None] = []
    if explicit:
        candidates.append(explicit)
    candidates.extend(_bundled_candidates(skill_root))
    configured = os.environ.get("PLINK_BIN")
    if configured:
        candidates.append(configured)
    candidates.extend([shutil.which("plink"), shutil.which("plink1.9")])

    failures: list[str] = []
    seen: set[str] = set()
    for candidate in candidates:
        if not candidate:
            continue
        candidate_str = str(candidate)
        found = shutil.which(candidate_str)
        path = Path(found or candidate_str).expanduser()
        try:
            path = path.resolve()
        except OSError:
            path = path.absolute()
        key = os.path.normcase(str(path))
        if key in seen:
            continue
        seen.add(key)
        if not path.is_file():
            continue
        try:
            _validate_plink(path)
            return str(path)
        except (OSError, ValueError, subprocess.SubprocessError) as exc:
            failures.append(f"{path}: {exc}")

    preferred = Path(skill_root).expanduser().resolve() if skill_root else Path(__file__).resolve().parents[1]
    preferred = preferred / "tools" / "plink" / "windows-x64" / "plink.exe"
    detail = (" Validation failures: " + "; ".join(failures)) if failures else ""
    raise RuntimeError(
        "No usable genomic PLINK 1.9 was found. On Windows, place plink.exe at "
        f"'{preferred}', or provide an existing executable through PLINK_BIN, PATH, or --plink-cmd. "
        "The skill does not auto-install PLINK or write to the user profile."
        + detail
    )


def run_plink(args, skill_root=None, explicit=None, auto_install=False):
    binary = resolve_plink(skill_root, explicit, auto_install=auto_install)
    subprocess.run([binary, *map(str, args)], check=True)
    return binary
