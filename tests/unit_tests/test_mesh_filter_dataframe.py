"""Mesh filter dataframes must match the mesh's own element indices.

`MeshFilter` and `MeshSurfaceFilter` build their index columns from
`mesh.indices` rather than reconstructing them per axis. These tests pin the
column names, the ordering, and the index base for every mesh type, so that a
change to `indices` cannot silently alter every user's dataframe.
"""

import numpy as np
import openmc
import pytest


@pytest.fixture
def regular_3d():
    mesh = openmc.RegularMesh()
    mesh.lower_left = (-1.0, -1.0, -1.0)
    mesh.upper_right = (1.0, 1.0, 1.0)
    mesh.dimension = (2, 3, 4)
    return mesh


@pytest.fixture
def regular_2d():
    mesh = openmc.RegularMesh()
    mesh.lower_left = (-1.0, -1.0)
    mesh.upper_right = (1.0, 1.0)
    mesh.dimension = (3, 2)
    return mesh


@pytest.fixture
def rectilinear():
    mesh = openmc.RectilinearMesh()
    mesh.x_grid = [-1.0, 0.0, 1.0]
    mesh.y_grid = [-1.0, 0.0, 2.0]
    mesh.z_grid = [-1.0, 1.0]
    return mesh


@pytest.fixture
def cylindrical():
    return openmc.CylindricalMesh(
        r_grid=[0.0, 1.0, 2.0],
        phi_grid=[0.0, np.pi, 2 * np.pi],
        z_grid=[-1.0, 0.0, 1.0],
    )


@pytest.fixture
def spherical():
    return openmc.SphericalMesh(
        r_grid=[0.0, 1.0, 2.0],
        theta_grid=[0.0, np.pi / 2, np.pi],
        phi_grid=[0.0, np.pi, 2 * np.pi],
    )


ALL_MESHES = ('regular_3d', 'regular_2d', 'rectilinear', 'cylindrical',
              'spherical')

EXPECTED_LABELS = {
    'regular_3d': ('x', 'y', 'z'),
    'regular_2d': ('x', 'y'),
    'rectilinear': ('x', 'y', 'z'),
    'cylindrical': ('r', 'phi', 'z'),
    'spherical': ('r', 'theta', 'phi'),
}


@pytest.mark.parametrize('mesh_name', ALL_MESHES)
def test_mesh_filter_columns(mesh_name, request):
    """Column names come from the mesh's own axis labels."""
    mesh = request.getfixturevalue(mesh_name)
    filt = openmc.MeshFilter(mesh)
    df = filt.get_pandas_dataframe(filt.num_bins, 1)

    key = f'mesh {mesh.id}'
    assert list(df.columns) == [(key, ax) for ax in EXPECTED_LABELS[mesh_name]]
    assert len(df) == filt.num_bins


@pytest.mark.parametrize('mesh_name', ALL_MESHES)
def test_mesh_filter_matches_mesh_indices(mesh_name, request):
    """Rows are exactly the mesh's element indices, in bin order.

    This is the property the implementation relies on. Structured meshes report
    1-based indices, so the first row is all ones and the last row is the
    dimension tuple.
    """
    mesh = request.getfixturevalue(mesh_name)
    filt = openmc.MeshFilter(mesh)
    df = filt.get_pandas_dataframe(filt.num_bins, 1)

    expected = list(mesh.indices)
    assert [tuple(row) for row in df.to_numpy()] == [tuple(i) for i in expected]

    # Structured meshes are 1-based
    assert tuple(df.to_numpy()[0]) == tuple([1] * len(EXPECTED_LABELS[mesh_name]))
    assert tuple(df.to_numpy()[-1]) == tuple(mesh.dimension)


@pytest.mark.parametrize('mesh_name', ALL_MESHES)
def test_mesh_surface_filter_columns(mesh_name, request):
    """Surface filter adds a 'surf' column after the axis columns."""
    mesh = request.getfixturevalue(mesh_name)
    filt = openmc.MeshSurfaceFilter(mesh)
    df = filt.get_pandas_dataframe(filt.num_bins, 1)

    key = f'mesh {mesh.id}'
    labels = EXPECTED_LABELS[mesh_name]
    assert list(df.columns) == [(key, ax) for ax in labels] + [(key, 'surf')]
    assert len(df) == filt.num_bins


@pytest.mark.parametrize('mesh_name', ALL_MESHES)
def test_mesh_surface_filter_bin_order(mesh_name, request):
    """Surface names cycle fastest, element indices slowest."""
    mesh = request.getfixturevalue(mesh_name)
    filt = openmc.MeshSurfaceFilter(mesh)
    df = filt.get_pandas_dataframe(filt.num_bins, 1)

    key = f'mesh {mesh.id}'
    labels = EXPECTED_LABELS[mesh_name]
    n_surfs = 4 * len(labels)

    surfs = list(df[(key, 'surf')])
    assert surfs[:4] == [f'{labels[0]}-min out', f'{labels[0]}-min in',
                         f'{labels[0]}-max out', f'{labels[0]}-max in']
    # The surface column repeats with period n_surfs
    assert surfs[:n_surfs] == surfs[n_surfs:2 * n_surfs]

    # The element index is constant across one element's surface bins
    first = df.iloc[:n_surfs][[(key, ax) for ax in labels]].to_numpy()
    assert (first == first[0]).all()


def test_unstructured_mesh_filter_is_zero_based():
    """Unstructured meshes stay 0-based with a single element_index column.

    This is the case the removed `idx_start` special case in MeshFilter existed
    for. UnstructuredMesh.indices is already 0-based, so delegating to it gives
    the same result without the isinstance check.
    """
    mesh = openmc.UnstructuredMesh('dummy.exo', 'moab')
    # Stand in for data that would normally come from a statepoint
    mesh._has_statepoint_data = True
    mesh.n_elements = 4
    mesh._volumes = np.ones(4)

    filt = openmc.MeshFilter(mesh)
    df = filt.get_pandas_dataframe(filt.num_bins, 1)

    key = f'mesh {mesh.id}'
    assert list(df.columns) == [(key, 'element_index')]
    assert [tuple(row) for row in df.to_numpy()] == [(0,), (1,), (2,), (3,)]
