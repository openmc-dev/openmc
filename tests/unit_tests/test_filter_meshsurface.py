from math import pi
from pathlib import Path

import openmc


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
