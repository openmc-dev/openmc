from pathlib import Path

import pytest

import openmc
from openmc.utility_funcs import input_path, set_xml_input_path


class _EmptyPath:
    def __fspath__(self):
        return ''


@pytest.mark.parametrize('filename', ['', _EmptyPath()])
def test_input_path_empty_not_resolved(filename):
    """Missing filenames stay unresolved, including path-like inputs."""
    with openmc.config.patch('resolve_paths', True):
        assert input_path(filename) == Path()


@pytest.mark.parametrize('filename', ['foo.h5', Path('foo.h5')])
def test_input_path_nonempty_resolved(filename, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with openmc.config.patch('resolve_paths', True):
        assert input_path(filename) == tmp_path / 'foo.h5'


def test_input_path_xml_context(tmp_path):
    with openmc.config.patch('resolve_paths', True):
        with set_xml_input_path(tmp_path / 'model.xml'):
            assert input_path('') == Path()
            assert input_path('foo.h5') == tmp_path / 'foo.h5'
