import openmc
import pytest

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def box_model():
    model = openmc.model.Model()

    # Define materials
    water = openmc.Material()
    water.add_nuclide("H1", 2.0)
    water.add_nuclide("O16", 1.0)
    water.add_s_alpha_beta("c_H_in_H2O")
    water.set_density("g/cc", 1.0)

    fuel = openmc.Material()
    fuel.add_nuclide("U235", 1.0)
    fuel.set_density("g/cc", 4.5)

    # Define geometry
    x_periodic_rot = openmc.XPlane(
        surface_id=1,
        x0=0.0,
        boundary_type="transformation",
        transformation={
            "direction": [
                0.0, 1.0, 0.0, 0.0,
                -1.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0,
            ],
            "position": [
                0.0, 1.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0,
            ]
        },
    )
    y_periodic_rot = openmc.YPlane(
        surface_id=3,
        y0=0.0,
        boundary_type="transformation",
        transformation={
            "direction": [
                0.0, -1.0, 0.0, 0.0,
                1.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0,
            ],
            "position": [
                0.0, 0.0, 0.0, 0.0,
                1.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0,
            ]
        },
    )

    x_reflective = openmc.XPlane(
        surface_id=2,
        x0=5.0,
        boundary_type="transformation",
        transformation={
            "direction": [
                -1.0, 0.0, 0.0, 0.0,
                0.0, 1.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0,
            ]
        },
    )
    y_reflective = openmc.YPlane(
        surface_id=4,
        y0=5.0,
        boundary_type="transformation",
        transformation={
            "direction": [
                1.0, 0.0, 0.0, 0.0,
                0.0, -1.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0,
            ]
        },
    )

    z_periodic_trans_1 = openmc.ZPlane(
        surface_id=5,
        z0=-5.0,
        boundary_type="transformation",
        transformation={
            "position": [
                1.0, 0.0, 0.0, 0.0,
                0.0, 1.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 10.0,
            ]
        },
    )
    z_periodic_trans_2 = openmc.Plane(
        surface_id=6,
        a=0,
        b=0,
        c=1,
        d=5.0,
        boundary_type="transformation",
        transformation={
            "position": [
                1.0, 0.0, 0.0, 0.0,
                0.0, 1.0, 0.0, 0.0,
                0.0, 0.0, 1.0, -10.0,
            ]
        },
    )
    z_cyl = openmc.ZCylinder(surface_id=7, x0=2.5, y0=0.0, r=2.0)

    outside_cyl = openmc.Cell(
        1,
        fill=water,
        region=(
            +x_periodic_rot
            & -x_reflective
            & +y_periodic_rot
            & -y_reflective
            & +z_periodic_trans_1
            & -z_periodic_trans_2
            & +z_cyl
        ),
    )
    inside_cyl = openmc.Cell(
        2,
        fill=fuel,
        region=(
            +y_periodic_rot
            & +z_periodic_trans_1
            & -z_periodic_trans_2
            & -z_cyl
        ),
    )
    root_universe = openmc.Universe(0, cells=(outside_cyl, inside_cyl))
    model.geometry = openmc.Geometry(root_universe)

    # Define settings
    model.settings.particles = 1000
    model.settings.batches = 4
    model.settings.inactive = 0
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Box((0, 0, 0), (5, 5, 0))
    )

    return model


def test_transformation_mixed(box_model):
    harness = PyAPITestHarness("statepoint.4.h5", box_model)
    harness.main()
