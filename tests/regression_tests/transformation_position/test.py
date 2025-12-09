import openmc
import pytest
import numpy as np

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def hexagon_model():
    model = openmc.model.Model()

    fuel = openmc.Material()
    fuel.add_nuclide('U235', 1.0)
    fuel.set_density('g/cc', 4.5)

    surface_1 = openmc.XPlane(
        x0=8.660254037844386,
        boundary_type="transformation",
        transformation={"position":
            [1.0, 0.0, 0.0, 2 * -8.660254037844386,
             0.0, 1.0, 0.0, 0.0,
             0.0, 0.0, 1.0, 0.0]}
    )
    surface_2 = openmc.XPlane(
        x0=-8.660254037844386,
        boundary_type="transformation",
        transformation={"position":
            [1.0, 0.0, 0.0, 2 * 8.660254037844386,
             0.0, 1.0, 0.0, 0.0,
             0.0, 0.0, 1.0, 0.0]}
    )
    surface_3 = openmc.Plane(
        a=0.5773502691896257,
        b=1.0,
        c=0.0,
        d=10.0,
        boundary_type="transformation",
        transformation={"position":
            [1.0, 0.0, 0.0, -5.0,
             0.0, 1.0, 0.0, -8.660254037844386,
             0.0, 0.0, 1.0, 0.0]}
    )
    surface_4 = openmc.Plane(
        a=-0.5773502691896257,
        b=1.0,
        c=0.0,
        d=10.0,
        boundary_type="transformation",
        transformation={"position":
            [1.0, 0.0, 0.0, 5.0,
             0.0, 1.0, 0.0, -8.660254037844386,
             0.0, 0.0, 1.0, 0.0]}
    )
    surface_5 = openmc.Plane(
        a=-0.5773502691896257,
        b=1.0,
        c=0.0,
        d=-10.0,
        boundary_type="transformation",
        transformation={"position":
            [1.0, 0.0, 0.0, -5.0,
             0.0, 1.0, 0.0, 8.660254037844386,
             0.0, 0.0, 1.0, 0.0]}
    )
    surface_6 = openmc.Plane(
        a=0.5773502691896257,
        b=1.,
        c=0.0,
        d=-10.0,
        boundary_type="transformation",
        transformation={"position":
            [1.0, 0.0, 0.0, 5.0,
             0.0, 1.0, 0.0, 8.660254037844386,
             0.0, 0.0, 1.0, 0.0]}
    )

    region = (
        -surface_1
        & +surface_2
        & -surface_3
        & -surface_4
        & +surface_5
        & +surface_6
    )
    cell = openmc.Cell(fill=fuel, region=region)

    model.geometry = openmc.Geometry([cell])

    model.settings.particles = 1000
    model.settings.batches = 5
    model.settings.inactive = 0

    return model


def test_transformation_position(hexagon_model):
    harness = PyAPITestHarness('statepoint.5.h5', hexagon_model)
    harness.main()
