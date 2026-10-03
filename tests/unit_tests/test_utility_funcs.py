from pathlib import Path

import pytest

import openmc
from openmc.utility_funcs import input_path, set_xml_input_path


class _PathLike:
    def __init__(self, filename):
        self.filename = filename

    def __fspath__(self):
        return self.filename


@pytest.mark.parametrize('filename', ['', _PathLike('')])
@pytest.mark.parametrize('resolve_paths', [True, False])
@pytest.mark.parametrize('xml_context', [True, False])
def test_input_path_empty_not_resolved(
    filename, resolve_paths, xml_context, tmp_path
):
    """An empty filename represents "no file" (e.g., an in-memory mesh or
    DAGMC universe read back from a summary/statepoint that was never
    written with a source file, see UnstructuredMesh.from_hdf5 and
    DAGMCUniverse.from_hdf5) and should stay empty rather than being
    resolved to the current working directory."""
    with openmc.config.patch('resolve_paths', resolve_paths):
        if xml_context:
            with set_xml_input_path(tmp_path / 'model' / 'model.xml'):
                assert input_path(filename) == Path()
        else:
            assert input_path(filename) == Path()


@pytest.mark.parametrize('filename', [
    'foo.h5', Path('foo.h5'), _PathLike('foo.h5')
])
@pytest.mark.parametrize('resolve_paths', [True, False])
@pytest.mark.parametrize('xml_context', [True, False])
def test_input_path_nonempty_resolved(
    filename, resolve_paths, xml_context, tmp_path, monkeypatch
):
    monkeypatch.chdir(tmp_path)
    with openmc.config.patch('resolve_paths', resolve_paths):
        if xml_context:
            with set_xml_input_path(tmp_path / 'model' / 'model.xml'):
                expected = tmp_path / 'model' / 'foo.h5'
                if not resolve_paths:
                    expected = Path('foo.h5')
                assert input_path(filename) == expected
        else:
            expected = tmp_path / 'foo.h5' if resolve_paths else Path('foo.h5')
            assert input_path(filename) == expected
