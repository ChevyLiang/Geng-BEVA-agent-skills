import importlib.util
import sys
from pathlib import Path


def engine():
    name="_geng_beva_transgene_engine"
    if name in sys.modules:return sys.modules[name]
    path=Path(__file__).resolve().parents[2]/"transgene-detection"/"scripts"/"transgene_engine.py"
    if not path.is_file():raise FileNotFoundError("Legacy transgene commands require sibling transgene-detection skill")
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module)
    return module
