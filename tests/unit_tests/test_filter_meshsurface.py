from math import pi
from pathlib import Path

import numpy as np
import openmc
import pytest


@pytest.mark.parametrize('mesh_type', ('regular', 'rectilinear', 'cylindrical', 'spherical'))
def test_axis_labels(mesh_type):
    if mesh_type == 'regular':
        mesh = openmc.RegularMesh()
        mesh.lower_left = np.asarray([0.0]*3)
        mesh.width = np.asarray([0.0]*3)
        mesh.dimension = (1, 1, 1)
    elif mesh_type == 'rectilinear':
        mesh = openmc.RectilinearMesh()
        mesh.x_grid = np.asarray([0.0, 1.0])
        mesh.y_grid = np.asarray([0.0, 1.0])
        mesh.z_grid = np.asarray([0.0, 1.0])
    elif mesh_type == 'cylindrical':
        r_grid = np.asarray([0.0, 1.0])
        z_grid = np.asarray([0.0, 1.0])
        p_grid = np.asarray([0.0, 1.0])
        mesh = openmc.CylindricalMesh(r_grid=r_grid, z_grid=z_grid, phi_grid=p_grid)
    elif mesh_type == 'spherical':
        r_grid = np.asarray([0.0, 1.0])
        t_grid = np.asarray([0.0, 1.0])
        p_grid = np.asarray([0.0, 1.0])
        mesh = openmc.SphericalMesh(r_grid=r_grid, theta_grid=t_grid, phi_grid=p_grid)
    else:
        raise ValueError(mesh_type)

    filt = openmc.MeshSurfaceFilter(mesh)
    assert len(filt.bins) == 12
    bin_names = [b[3] for b in filt.bins]
    if mesh_type in {'regular', 'rectilinear'}:
        assert bin_names[:8] == ['x-min out', 'x-min in', 'x-max out', 'x-max in',
                                 'y-min out', 'y-min in', 'y-max out', 'y-max in']
    if mesh_type in {'regular', 'rectilinear', 'cylindrical'}:
        assert bin_names[8:] == ['z-min out', 'z-min in', 'z-max out', 'z-max in']
    if mesh_type in {'cylindrical', 'spherical'}:
        assert bin_names[:4] == ['r-min out', 'r-min in', 'r-max out', 'r-max in']
    if mesh_type == 'spherical':
        assert bin_names[4:] == ['theta-min out', 'theta-min in', 'theta-max out', 'theta-max in',
                                 'phi-min out', 'phi-min in', 'phi-max out', 'phi-max in']
    if mesh_type == 'cylindrical':
        assert bin_names[4:8] == ['phi-min out', 'phi-min in', 'phi-max out', 'phi-max in']



def test_mesh_surface_labels(run_in_tmpdir):
    """Check that mesh surface labels use the mesh coordinate system."""
    sphere = openmc.Sphere(r=1.0, boundary_type='vacuum')
    model = openmc.Model(geometry=openmc.Geometry([openmc.Cell(region=-sphere)]))
    model.settings.run_mode = 'fixed source'
    model.settings.batches = 1
    model.settings.particles = 1
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point((0.1, 0.1, 0.1))
    )

    mesh_1d = openmc.RegularMesh()
    mesh_1d.dimension = (2,)
    mesh_1d.lower_left = (-1.0,)
    mesh_1d.upper_right = (1.0,)

    mesh_2d = openmc.RegularMesh()
    mesh_2d.dimension = (2, 1)
    mesh_2d.lower_left = (-1.0, -1.0)
    mesh_2d.upper_right = (1.0, 1.0)

    mesh_3d = openmc.RegularMesh()
    mesh_3d.dimension = (2, 1, 1)
    mesh_3d.lower_left = (-1.0, -1.0, -1.0)
    mesh_3d.upper_right = (1.0, 1.0, 1.0)

    rectilinear_mesh = openmc.RectilinearMesh()
    rectilinear_mesh.x_grid = (-1.0, 0.0, 1.0)
    rectilinear_mesh.y_grid = (-1.0, 1.0)
    rectilinear_mesh.z_grid = (-1.0, 1.0)
    cylindrical_mesh = openmc.CylindricalMesh(
        r_grid=(0.0, 0.5, 1.0),
        phi_grid=(0.0, 2*pi),
        z_grid=(-1.0, 1.0),
    )
    spherical_mesh = openmc.SphericalMesh(
        r_grid=(0.0, 0.5, 1.0),
        theta_grid=(0.0, pi),
        phi_grid=(0.0, 2*pi),
    )

    tally_meshes = {}
    for mesh in (
        mesh_1d, mesh_2d, mesh_3d, rectilinear_mesh, cylindrical_mesh,
        spherical_mesh
    ):
        tally = openmc.Tally()
        tally.filters = [openmc.MeshSurfaceFilter(mesh)]
        tally.scores = ['current']
        model.tallies.append(tally)
        tally_meshes[tally.id] = mesh

    model.run()
    output = Path('tallies.out').read_text()

    for tally_id, mesh in tally_meshes.items():
        header = f'TALLY {tally_id}'
        start = output.index(header) + len(header)
        end = output.find('TALLY', start)
        section = output[start:] if end == -1 else output[start:end]

        n_dim = mesh.n_dimension
        labels = [
            line for line in section.splitlines() if 'Mesh Index' in line
        ][:4*n_dim + 1]

        expected = []
        for axis in mesh.axis_labels:
            expected.extend([
                f'Outgoing, {axis}-min',
                f'Incoming, {axis}-min',
                f'Outgoing, {axis}-max',
                f'Incoming, {axis}-max',
            ])
        expected.append(expected[0])

        assert len(labels) == len(expected)
        for line, label in zip(labels, expected):
            assert label in line
