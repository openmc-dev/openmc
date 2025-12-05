import openmc
import pytest
import numpy as np

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def box_model():
    model = openmc.model.Model()

    fuel = openmc.Material()
    fuel.add_nuclide('U235', 1.0)
    fuel.set_density('g/cc', 4.5)

    normal_x = np.array([1.0, 0.0, 0.0])
    dir_transform_xplane = np.identity(3) - 2 * np.outer(
        normal_x, normal_x
    )
    dir_transform_xplane = np.append(dir_transform_xplane, np.zeros((3,1)), axis=1).flatten()
    
    normal_y = np.array([0.0, 1.0, 0.0])
    dir_transform_yplane = np.identity(3) - 2 * np.outer(
        normal_y, normal_y
    )
    dir_transform_yplane = np.append(dir_transform_yplane, np.zeros((3,1)), axis=1).flatten()
    
    normal_z = np.array([0.0, 0.0, 1.0])
    dir_transform_zplane = np.identity(3) - 2 * np.outer(
        normal_z, normal_z
    )
    dir_transform_zplane = np.append(dir_transform_zplane, np.zeros((3,1)), axis=1).flatten()
    
    neg_xplane = openmc.XPlane(
        x0=-5.0,
        boundary_type="transformation",
        transformation={"direction": dir_transform_xplane}
    )
    pos_xplane = openmc.XPlane(
        x0=5.0,
        boundary_type="transformation",
        transformation={"direction": dir_transform_xplane}
    )

    neg_yplane = openmc.YPlane(
        y0=-5.0,
        boundary_type="transformation",
        transformation={"direction": dir_transform_yplane}
    )
    pos_yplane = openmc.YPlane(
        y0=5.0,
        boundary_type="transformation",
        transformation={"direction": dir_transform_yplane}
    )

    neg_zplane = openmc.ZPlane(
        z0=-5.0,
        boundary_type="transformation",
        transformation={"direction": dir_transform_zplane}
    )
    pos_zplane = openmc.ZPlane(
        z0=5.0,
        boundary_type="transformation",
        transformation={"direction": dir_transform_zplane}
    )

    region = (
        +neg_xplane
        & -pos_xplane
        & +neg_yplane
        & -pos_yplane
        & +neg_zplane
        & -pos_zplane
    )
    cell = openmc.Cell(fill=fuel, region=region)

    model.geometry = openmc.Geometry([cell])

    model.settings.particles = 1000
    model.settings.batches = 10
    model.settings.inactive = 5

    source = openmc.IndependentSource()
    source.space = openmc.stats.Box([-4.0, -4.0, -4.0], [4.0, 4.0, 4.0])
    model.settings.source = source

    return model


def test_transformation_direction(box_model):
    harness = PyAPITestHarness('statepoint.10.h5', box_model)
    harness.main()
