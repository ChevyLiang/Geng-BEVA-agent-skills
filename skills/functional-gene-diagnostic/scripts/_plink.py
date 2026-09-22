from __future__ import annotations

import os
import platform
import re
import shutil
import stat
import subprocess
import tempfile
import urllib.request
import zipfile
from pathlib import Path

PLINK_HOME = Path(os.environ.get("GENG_BEVA_TOOL_HOME", Path.home() / ".geng-beva" / "tools"))
PLINK_PAGE = "https://www.cog-genomics.org/plink/"
PLINK_FALLBACK_DATE = "20260909"


def _platform_key() -> tuple[str, str]:
    system = platform.system().lower()
    machine = platform.machine().lower()
    if system == "windows":
        return "windows", "plink_win64" if machine in {"amd64", "x86_64"} else "plink_win32"
    if system == "darwin":
        return "macos", "plink_mac"
    if system == "linux":
        if machine in {"x86_64", "amd64"}:
            return "linux", "plink_linux_x86_64"
        if machine in {"i386", "i686", "x86"}:
            return "linux", "plink_linux_i686"
    raise RuntimeError(f"Automatic PLINK installation is not supported on {platform.system()} {platform.machine()}.")


def _official_download_url(asset_prefix: str) -> str:
    """Resolve the current stable PLINK 1.9 ZIP from the official download page.

    Falls back to the repository-tested stable build if the page cannot be queried.
    """
    try:
        req = urllib.request.Request(PLINK_PAGE, headers={"User-Agent": "geng-beva-skills/1.0"})
        with urllib.request.urlopen(req, timeout=20) as response:
            html = response.read().decode("utf-8", errors="ignore")
        pattern = rf"https://s3\.amazonaws\.com/plink1-assets/{re.escape(asset_prefix)}_(\d{{8}})\.zip"
        matches = re.findall(pattern, html)
        if matches:
            # The first/lowest date is not guaranteed to be stable; choose newest matching official build.
            date = sorted(set(matches), reverse=True)[0]
            return f"https://s3.amazonaws.com/plink1-assets/{asset_prefix}_{date}.zip"
    except Exception:
        pass
    return f"https://s3.amazonaws.com/plink1-assets/{asset_prefix}_{PLINK_FALLBACK_DATE}.zip"


def _validate_plink(path: Path) -> str:
    check = subprocess.run([str(path), "--version"], capture_output=True, text=True, timeout=20, check=True)
    version = check.stdout + check.stderr
    if not re.search(r"PLINK\s+v?1\.9", version, re.I):
        raise ValueError("Expected genomic PLINK 1.9 (not PuTTY plink or PLINK 2)")
    return version.strip().splitlines()[0] if version.strip() else "PLINK 1.9"


def install_plink(cache_root: str | Path | None = None) -> str:
    os_name, asset_prefix = _platform_key()
    root = Path(cache_root).expanduser() if cache_root else PLINK_HOME
    install_dir = root / "plink-1.9" / os_name
    binary = install_dir / ("plink.exe" if os_name == "windows" else "plink")
    if binary.is_file():
        try:
            _validate_plink(binary)
            return str(binary.resolve())
        except Exception:
            binary.unlink(missing_ok=True)

    install_dir.mkdir(parents=True, exist_ok=True)
    url = _official_download_url(asset_prefix)
    with tempfile.TemporaryDirectory(prefix="geng_beva_plink_") as tmp:
        archive = Path(tmp) / "plink.zip"
        try:
            urllib.request.urlretrieve(url, archive)
        except Exception as exc:
            raise RuntimeError(
                "PLINK 1.9 is required but was not found, and automatic download failed. "
                f"Attempted official URL: {url}. Set PLINK_BIN or pass --plink-cmd to use an existing binary."
            ) from exc
        try:
            with zipfile.ZipFile(archive) as zf:
                member = next((n for n in zf.namelist() if Path(n).name.lower() in {"plink", "plink.exe"}), None)
                if member is None:
                    raise RuntimeError("The downloaded PLINK archive did not contain a PLINK executable.")
                extracted = Path(zf.extract(member, tmp))
                shutil.copy2(extracted, binary)
        except zipfile.BadZipFile as exc:
            raise RuntimeError("The downloaded PLINK archive was not a valid ZIP file.") from exc

    if os_name != "windows":
        binary.chmod(binary.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    _validate_plink(binary)
    return str(binary.resolve())


def resolve_plink(skill_root=None, explicit=None, auto_install=True):
    configured = explicit or os.environ.get("PLINK_BIN")
    cached_candidates = [
        PLINK_HOME / "plink-1.9" / "windows" / "plink.exe",
        PLINK_HOME / "plink-1.9" / "linux" / "plink",
        PLINK_HOME / "plink-1.9" / "macos" / "plink",
    ]
    candidates = [configured] if configured else [shutil.which("plink"), shutil.which("plink1.9"), *cached_candidates]
    failures = []

    for candidate in candidates:
        if not candidate:
            continue
        found = shutil.which(str(candidate))
        path = Path(found or candidate).expanduser().resolve()
        if not path.is_file():
            continue
        try:
            _validate_plink(path)
            return str(path)
        except (OSError, ValueError, subprocess.SubprocessError) as exc:
            failures.append(f"{path}: {exc}")

    if auto_install and not explicit and os.environ.get("GENG_BEVA_NO_AUTO_INSTALL", "").lower() not in {"1", "true", "yes"}:
        return install_plink()

    detail = "; ".join(failures)
    if detail:
        detail = " " + detail
    raise RuntimeError(
        "No usable genomic PLINK 1.9 was found. Automatic installation is disabled or unavailable. "
        "Set PLINK_BIN or pass --plink-cmd." + detail
    )


def run_plink(args, skill_root=None, explicit=None, auto_install=True):
    binary = resolve_plink(skill_root, explicit, auto_install=auto_install)
    subprocess.run([binary, *map(str, args)], check=True)
    return binary
