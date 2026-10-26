"""Test the dnp_drift Python object."""

import openmc
import pytest
import numpy as np


@pytest.fixture()
def mesh():
    """2x2x2 regular mesh"""
    dim = 2
    lower_left = (0.0, 0.0, 0.0)
    upper_right = (5.0, 5.0, 5.0)
    mesh = openmc.RegularMesh()
    mesh.lower_left = lower_left
    mesh.upper_right = upper_right
    mesh.dimension = (dim, dim, dim)
    return mesh


def test_xml_serialization(mesh, run_in_tmpdir):
    """Test XML serialization."""
    velocity_field = openmc.VelocityField(
        mesh=mesh, values=np.tile([0, 0, 10], (8, 1)), mapping="cell"
    )

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
        external_travel_time=2.0,
    )

    # Export settings
    settings = openmc.Settings()
    settings.dnp_drift = dnp_drift
    settings.export_to_xml()

    # Read settings from xml file
    read_settings = openmc.Settings.from_xml()

    # Check consistency
    assert read_settings.dnp_drift.velocity_field == velocity_field
    assert read_settings.dnp_drift.boundary_map == {
        "inlet": [1],
        "outlet": [2],
        "wall": [3, 4],
    }
    assert read_settings.dnp_drift.physical_group_map == {
        "face_ids": [
            0, 24, 12, 36, 7, 31, 19, 43, 15, 39, 21, 45, 2,
            26, 8,32, 4, 16, 10, 22, 29, 41, 35, 47
        ],
        "physical_groups": [
            1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3
        ]}
    assert read_settings.dnp_drift.integrator == "RK4"
    assert read_settings.dnp_drift.integrator_dt == 0.1
    assert read_settings.dnp_drift.external_travel_time == 2.0


def test_invalid_mesh_type():
    """Check one non-implemented mesh type."""
    with pytest.raises(TypeError):
        velocity_field = openmc.VelocityField(
            openmc.SphericalMesh(r_grid=[0.0, 1.0]), []
        )


def test_inconsistent_values_and_cells_sizes(mesh):
    """Check inconsistency in number of values / mesh cells."""
    with pytest.raises(ValueError):
        velocity_field = openmc.VelocityField(
            mesh=mesh, values=np.tile([0, 0, 10], (9, 1))
        )


def test_missing_boundary_map_keys(mesh):
    """Check missing required boundary map keys."""
    velocity_field = openmc.VelocityField(mesh=mesh, values=np.tile([0, 0, 10], (8, 1)))

    with pytest.raises(ValueError):
        openmc.DNPDrift(
            velocity_field=velocity_field,
            boundary_map={"inlet": [1]},
            physical_group_map={"face_ids": [0], "physical_groups": [1]}
        )


def test_invalid_boundary_map_key(mesh):
    """Check invalid boundary map keys."""
    velocity_field = openmc.VelocityField(mesh=mesh, values=np.tile([0, 0, 10], (8, 1)))

    with pytest.raises(ValueError):
        openmc.DNPDrift(
            velocity_field=velocity_field,
            boundary_map={"inlet": [1], "outlet": [2], "dnp_drift": [3]},
            physical_group_map={"face_ids": [0], "physical_groups": [1]}
        )


def test_invalid_integrator_dt(mesh):
    """Check negative integrator time step value."""
    velocity_field = openmc.VelocityField(mesh=mesh, values=np.tile([0, 0, 10], (8, 1)))

    with pytest.raises(ValueError):
        openmc.DNPDrift(
            velocity_field=velocity_field,
            boundary_map={"inlet": [1], "outlet": [2]},
            physical_group_map={"face_ids": [0], "physical_groups": [1]},
            integrator_dt=-0.1
        )


def test_invalid_external_travel_time(mesh):
    """Check negative external travel time."""
    velocity_field = openmc.VelocityField(mesh=mesh, values=np.tile([0, 0, 10], (8, 1)))

    with pytest.raises(ValueError):
        openmc.DNPDrift(
            velocity_field=velocity_field,
            boundary_map={"inlet": [1], "outlet": [2]},
            physical_group_map={"face_ids": [0], "physical_groups": [1]},
            external_travel_time=-1.0
        )
