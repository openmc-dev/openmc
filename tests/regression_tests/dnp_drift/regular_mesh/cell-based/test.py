import openmc
import pytest
import numpy as np

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def dnp_drift_model():
    model = openmc.Model()

    # Material
    mat = openmc.Material()
    mat.add_nuclide("U235", 0.2)
    mat.add_nuclide("U238", 0.8)
    mat.add_nuclide("O16", 3.0)
    mat.add_nuclide("H1", 2.0)
    mat.set_density("g/cm3", 5.0)
    model.materials = openmc.Materials([mat])

    # Create mesh
    dim = 2
    lower_left = (0.0, 0.0, 0.0)
    upper_right = (10.0, 10.0, 10.0)
    mesh = openmc.RegularMesh()
    mesh.lower_left = lower_left
    mesh.upper_right = upper_right
    mesh.dimension = (dim, dim, dim)

    # Create the velocity field
    velocity_field = openmc.VelocityField(
        mesh=mesh,
        values=np.array([1.0]*3*8).reshape(-1,3),
        mapping="cell"
    )

    # Create geometry
    box = openmc.model.RectangularParallelepiped(
        xmin=lower_left[0], xmax=upper_right[0],
        ymin=lower_left[1], ymax=upper_right[1],
        zmin=lower_left[2], zmax=upper_right[2],
        boundary_type="reflective"
    )
    cell = openmc.Cell(fill=mat, region=-box)
    model.geometry = openmc.Geometry([cell])

    # DNP drift
    dnp_drift = openmc.DNPDrift(
        velocity_field=velocity_field,
        boundary_map={"inlet": [1], "outlet": [2], "wall": [3, 4]},
        physical_group_map={
            "face_ids": [
                0, 24, 12, 36, 7, 31, 19, 43, 15, 39, 21, 45, 2,
                26, 8,32, 4, 16, 10, 22, 29, 41, 35, 47
            ],
            "physical_groups": [
                1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
                3, 3, 3, 3, 3, 3, 3, 3
            ]},
        integrator="RK4",
        integrator_dt=0.1,
        recycling=True,
        external_travel_time=2.0
    )

    # Settings
    settings = openmc.Settings()
    settings.batches = 20
    settings.inactive = 2
    settings.particles = 4000
    settings.dnp_drift = dnp_drift
    spatial_dist = openmc.stats.Box(lower_left, upper_right)
    settings.source = openmc.IndependentSource(
        space=spatial_dist, constraints={"fissionable": True})
    model.settings = settings

    # Add tallies
    mesh_filter = openmc.MeshFilter(mesh)
    mesh_tally = openmc.Tally(name="total reaction rate")
    mesh_tally.filters = [mesh_filter]
    mesh_tally.scores = ["total"]
    tallies = openmc.Tallies([mesh_tally])
    model.tallies = tallies

    return model


def test_dnp_drift_regular_mesh_cell_based(dnp_drift_model):
    harness = PyAPITestHarness('statepoint.20.h5', dnp_drift_model)
    harness.main()
