import numpy as np
from math import pi

import openmc
import pytest

from tests.testing_harness import HashedPyAPITestHarness


class MeshSurfaceFilterTest(HashedPyAPITestHarness):

    def _compare_results(self):
        tally_meshes = {}
        with openmc.StatePoint(self.statepoint_name) as sp:
            for tally in sp.tallies.values():
                for filt in tally.filters:
                    if isinstance(filt, openmc.MeshSurfaceFilter):
                        tally_meshes[tally.id] = filt.mesh
                        break

        with open('tallies.out') as fh:
            text = fh.read()

        for tally_id, mesh in tally_meshes.items():
            # Snip the relevant text.
            header = f"TALLY {tally_id}"
            assert header in text
            text = text[text.find(header) + len(header):]
            end = text.find("TALLY")
            section = text if end == -1 else text[:end]
            # Define the relevant labels and ensure they wrap.
            nd = mesh.n_dimension
            labels = [ln for ln in section.splitlines() if "Mesh Index" in ln]
            labels = labels[:4*nd + 1]  # just the first 4, 8, or 12, plus one
            expected = []
            if isinstance(mesh, openmc.SphericalMesh):
                axis_labels = ('r', 'theta', 'phi')
            elif isinstance(mesh, openmc.CylindricalMesh):
                axis_labels = ('r', 'phi', 'z')
            else:
                axis_labels = ('x', 'y', 'z')
            for axis in axis_labels[:nd]:
                expected += [f'Outgoing, {axis}-min', f'Incoming, {axis}-min',
                             f'Outgoing, {axis}-max', f'Incoming, {axis}-max']
            expected.append(expected[0])
        assert len(labels) == len(expected)
        for line, exp in zip(labels, expected):
            assert exp in line, f"{exp} not in {line}"

        return super()._compare_results()


@pytest.fixture
def model():
    model = openmc.model.Model()

    fuel = openmc.Material()
    fuel.set_density('g/cm3', 10.0)
    fuel.add_nuclide('U235', 1.0)
    zr = openmc.Material()
    zr.set_density('g/cm3', 1.0)
    zr.add_nuclide('Zr90', 1.0)
    model.materials.extend([fuel, zr])

    box1 = openmc.model.RectangularPrism(10.0, 10.0)
    box2 = openmc.model.RectangularPrism(20.0, 20.0, boundary_type='reflective')
    top = openmc.ZPlane(z0=10.0, boundary_type='vacuum')
    bottom = openmc.ZPlane(z0=-10.0, boundary_type='vacuum')
    cell1 = openmc.Cell(fill=fuel, region=-box1 & +bottom & -top)
    cell2 = openmc.Cell(fill=zr, region=+box1 & -box2 & +bottom & -top)
    model.geometry = openmc.Geometry([cell1, cell2])

    model.settings.batches = 5
    model.settings.inactive = 0
    model.settings.particles = 1000

    # Create meshes
    mesh_1d = openmc.RegularMesh()
    mesh_1d.dimension = [5]
    mesh_1d.lower_left = [-7.5]
    mesh_1d.upper_right = [7.5]

    mesh_2d = openmc.RegularMesh()
    mesh_2d.dimension = [5, 5]
    mesh_2d.lower_left = [-7.5, -7.5]
    mesh_2d.upper_right = [7.5, 7.5]

    mesh_3d = openmc.RegularMesh()
    mesh_3d.dimension = [5, 5, 5]
    mesh_3d.lower_left = [-7.5, -7.5, -7.5]
    mesh_3d.upper_right = [7.5, 7.5, 7.5]
    dx = dy = dz = 15 / 5
    reg_mesh_exp_vols = np.full(mesh_3d.dimension, dx*dy*dz)
    np.testing.assert_equal(mesh_3d.volumes, reg_mesh_exp_vols)

    recti_mesh = openmc.RectilinearMesh()
    recti_mesh.x_grid = np.linspace(-7.5, 7.5, 18)
    recti_mesh.y_grid = np.linspace(-7.5, 7.5, 18)
    recti_mesh.z_grid = np.logspace(0, np.log10(7.5), 11)
    dx = dy = 15 / 17
    dz = np.diff(np.logspace(0, np.log10(7.5), 11))
    dxdy = np.full(recti_mesh.dimension[:2], dx*dy)
    recti_mesh_exp_vols = np.multiply.outer(dxdy, dz)
    np.testing.assert_allclose(recti_mesh.volumes, recti_mesh_exp_vols)

    cyl_mesh = openmc.CylindricalMesh(
        origin=(0, 0, -7.5),
        r_grid=np.linspace(0, 7.5, 18),
        phi_grid=np.linspace(0, 2*pi, 19),
        z_grid=np.linspace(0, 15, 17),
    )
    dr = 0.5 * np.diff(np.linspace(0, 7.5, 18)**2)
    dp = np.full(cyl_mesh.dimension[1], 2*pi / 18)
    dz = np.full(cyl_mesh.dimension[2], 15 / 16)
    drdp = np.outer(dr, dp)
    cyl_mesh_exp_vols = np.multiply.outer(drdp, dz)
    np.testing.assert_allclose(cyl_mesh.volumes, cyl_mesh_exp_vols)

    sph_mesh = openmc.SphericalMesh(
        r_grid=np.linspace(0, 7.5, 18),
        theta_grid=np.linspace(0, pi, 9),
        phi_grid=np.linspace(0, 2*pi, 19)
    )
    dr = np.diff(np.linspace(0, 7.5, 18)**3) / 3
    dt = np.diff(-np.cos(np.linspace(0, pi, 9)))
    dp = np.full(sph_mesh.dimension[2], 2*pi / 18)
    drdt = np.outer(dr, dt)
    sph_mesh_exp_vols = np.multiply.outer(drdt, dp)
    np.testing.assert_allclose(sph_mesh.volumes, sph_mesh_exp_vols)

    # Create filters
    reg_filters = [
        openmc.MeshFilter(mesh_1d),
        openmc.MeshFilter(mesh_2d),
        openmc.MeshFilter(mesh_3d),
        openmc.MeshFilter(recti_mesh),
        openmc.MeshFilter(cyl_mesh),
        openmc.MeshFilter(sph_mesh)
    ]
    surf_filters = [
        openmc.MeshSurfaceFilter(mesh_1d),
        openmc.MeshSurfaceFilter(mesh_2d),
        openmc.MeshSurfaceFilter(mesh_3d),
        openmc.MeshSurfaceFilter(recti_mesh),
        openmc.MeshSurfaceFilter(cyl_mesh),
        openmc.MeshSurfaceFilter(sph_mesh)
    ]

    # Create tallies
    for f1, f2 in zip(reg_filters, surf_filters):
        tally = openmc.Tally()
        tally.filters = [f1]
        tally.scores = ['total']
        model.tallies.append(tally)
        tally = openmc.Tally()
        tally.filters = [f2]
        tally.scores = ['current']
        model.tallies.append(tally)

    return model


def test_filter_mesh(model):
    harness = MeshSurfaceFilterTest('statepoint.5.h5', model)
    harness.main()
