"""Tests that the Python API reads empty strings written by OpenMC's C++ layer.

An empty string has no zero-size representation in HDF5, so OpenMC stores it as
a single null byte. Older versions wrote nothing at all, leaving the dataset
absent from the file; both forms have to keep working.
"""

import h5py
import numpy as np
import pytest

import openmc


def write_string(group, name, value):
    """Write a fixed-length string the way OpenMC's C++ layer does."""
    data = np.array(value.encode(), dtype=f'S{max(len(value), 1)}')
    return group.create_dataset(name, data=data)


def write_string_attr(obj, name, value):
    """Write a fixed-length string attribute the way OpenMC's C++ layer does."""
    data = np.array(value.encode(), dtype=f'S{max(len(value), 1)}')
    obj.attrs.create(name, data)


@pytest.mark.parametrize('value', ['', 'a', 'some name', 'abc'])
def test_string_dataset_round_trip(run_in_tmpdir, value):
    """The on-disk contract: what OpenMC writes is what h5py reads back."""
    with h5py.File('strings.h5', 'w') as f:
        write_string(f, 'name', value)

    with h5py.File('strings.h5', 'r') as f:
        assert 'name' in f
        assert f['name'][()].decode() == value


@pytest.mark.parametrize('value', ['', 'a', '/path/to/inputs/'])
def test_string_attribute_round_trip(run_in_tmpdir, value):
    with h5py.File('strings.h5', 'w') as f:
        write_string_attr(f, 'path', value)

    with h5py.File('strings.h5', 'r') as f:
        assert 'path' in f.attrs
        assert f.attrs['path'].decode() == value


def test_string_vector_round_trip(run_in_tmpdir):
    """Vectors of strings are null-padded to the longest entry."""
    names = ['U235', '', 'H1']
    with h5py.File('strings.h5', 'w') as f:
        f.create_dataset('nuclides', data=np.array(
            [n.encode() for n in names], dtype='S5'))

    with h5py.File('strings.h5', 'r') as f:
        assert [n.decode() for n in f['nuclides'][()]] == names


# Each reader below is exercised twice: once against a file with the dataset
# present but empty (what OpenMC writes now) and once with it absent (what
# older versions wrote, and what the reader guards were added for).
@pytest.fixture(params=['present', 'absent'])
def empty_name(request):
    """Return a function that writes an empty name, or doesn't."""
    def write(group):
        if request.param == 'present':
            write_string(group, 'name', '')
    return write


def test_surface_from_hdf5_empty_name(run_in_tmpdir, empty_name):
    with h5py.File('surfaces.h5', 'w') as f:
        group = f.create_group('surface 1')
        empty_name(group)
        write_string(group, 'geom_type', 'csg')
        write_string(group, 'type', 'x-plane')
        write_string(group, 'boundary_type', 'vacuum')
        group.create_dataset('coefficients', data=np.array([3.0]))

        surface = openmc.Surface.from_hdf5(group)

    assert surface.name == ''
    assert surface.id == 1
    assert surface.boundary_type == 'vacuum'
    assert surface.x0 == 3.0


def test_material_from_hdf5_empty_name(run_in_tmpdir, empty_name):
    with h5py.File('materials.h5', 'w') as f:
        group = f.create_group('material 7')
        empty_name(group)
        group.attrs['depletable'] = 0
        group.create_dataset('atom_density', data=0.06022)
        group.create_dataset('nuclides', data=np.array([b'H1'], dtype='S3'))
        group.create_dataset('nuclide_densities', data=np.array([0.06022]))

        material = openmc.Material.from_hdf5(group)

    assert material.name == ''
    assert material.id == 7


def test_mesh_from_hdf5_empty_name(run_in_tmpdir, empty_name):
    with h5py.File('meshes.h5', 'w') as f:
        group = f.create_group('mesh 3')
        empty_name(group)
        write_string(group, 'type', 'regular')
        group.create_dataset('dimension', data=np.array([2, 2, 2]))
        group.create_dataset('lower_left', data=np.array([0.0, 0.0, 0.0]))
        group.create_dataset('width', data=np.array([1.0, 1.0, 1.0]))

        mesh = openmc.MeshBase.from_hdf5(group)

    assert mesh.name == ''
    assert mesh.id == 3
    assert mesh.dimension == (2, 2, 2)
