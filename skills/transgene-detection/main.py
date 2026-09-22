import sys
from pathlib import Path
import runpy

_CODEX_ROOT = next(p for p in Path(__file__).resolve().parents if p.name == ".codex")
sys.path.insert(0, str(_CODEX_ROOT / "tools"))
from runtime_launcher import ensure_fixed_runtime
ensure_fixed_runtime(__file__)
from dependency_resolver import bootstrap_from_requirements
bootstrap_from_requirements(Path(__file__).resolve().parent / "requirements.txt")

sys.path.insert(0,str(Path(__file__).resolve().parent/"scripts"))
from transgene_detection import run
if __name__=="__main__":
    runpy.run_path(str(Path(__file__).resolve().parent/"scripts/transgene_detection.py"),run_name="__main__")
