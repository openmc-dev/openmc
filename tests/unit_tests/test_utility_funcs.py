from pathlib import Path

import openmc
from openmc.utility_funcs import input_path


def test_input_path_empty_not_resolved():
    """An empty filename represents "no file" (e.g., an in-memory mesh or
    DAGMC universe read back from a summary/statepoint that was never
    written with a source file, see UnstructuredMesh.from_hdf5 and
    DAGMCUniverse.from_hdf5) and should stay empty rather than being
    resolved to the current working directory."""
    with openmc.config.patch('resolve_paths', True):
        assert input_path('') == Path()


def test_input_path_nonempty_resolved(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with openmc.config.patch('resolve_paths', True):
        assert input_path('foo.h5') == (tmp_path / 'foo.h5').resolve()
