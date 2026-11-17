"""Load one editable study without importing any other study's code."""
from pathlib import Path
import hashlib
import importlib.util
import sys
import types

_CURRENT = None


def load_study(directory, output=None):
    """Load config.py and model.py from an arbitrary study directory."""
    global _CURRENT
    directory = Path(directory).resolve()
    package = '_openmc_study_' + hashlib.sha256(str(directory).encode()).hexdigest()[:12]
    module = types.ModuleType(package)
    module.__path__ = [str(directory)]
    sys.modules[package] = module
    loaded = {}
    for name in ('config', 'model'):
        file = directory / (name + '.py')
        spec = importlib.util.spec_from_file_location(package + '.' + name, file)
        child = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = child
        spec.loader.exec_module(child)
        loaded[name] = child
    _CURRENT = types.SimpleNamespace(directory=directory, config=loaded['config'],
        model=loaded['model'], output=Path(output).resolve() if output else directory/getattr(loaded['config'], 'RESULTS_DIRECTORY', 'results'))
    return _CURRENT


def settings():
    if _CURRENT is None:
        raise RuntimeError('Load a study with load_study(directory) first')
    return _CURRENT.config


def output_directory():
    settings()
    return _CURRENT.output


def build_model(particles=1000, batches=10):
    """Delegate model construction entirely to the selected study."""
    settings()
    return _CURRENT.model.build_model(particles, batches)


def responses(model):
    """Use a study-defined response or the verified production response."""
    if hasattr(_CURRENT.model, 'responses'):
        return _CURRENT.model.responses(model)
    from .production import responses as production_responses
    return production_responses(model)


def register(model, response_definitions):
    from .config import ModelSpec
    names = tuple(sorted({n for m in model.geometry.get_all_materials().values()
                          for n in m.get_nuclides()}))
    c = settings()
    targets = tuple(c.TARGET_NUCLIDES) if c.TARGET_NUCLIDES else names
    if not set(targets).issubset(names):
        raise ValueError('TARGET_NUCLIDES contains nuclides absent from the model')
    return ModelSpec(c.STUDY_NAME, model, response_definitions, targets,
                     threads=c.TRANSPORT_THREADS)
