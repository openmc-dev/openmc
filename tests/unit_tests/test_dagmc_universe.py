import openmc
import pytest


def test_bounding_box():

    u = openmc.DAGMCUniverse("tests/regression_tests/dagmc/universes/dagmc.h5m")

    ll, ur = u.bounding_box
    assert ll == pytest.approx((-25., -25., -25))
    assert ur == pytest.approx((25., 25., 25))


def test_bounding_region():

    u = openmc.DAGMCUniverse("tests/regression_tests/dagmc/universes/dagmc.h5m")

    region = u.bounding_region()  # should default to bounded_type='box'
    assert isinstance(region, openmc.Region)
    assert len(region) == 6
    assert region[0].surface.type == 'x-plane'
    assert region[1].surface.type == 'x-plane'
    assert region[2].surface.type == 'y-plane'
    assert region[3].surface.type == 'y-plane'
    assert region[4].surface.type == 'z-plane'
    assert region[5].surface.type == 'z-plane'
    assert region[0].surface.boundary_type == 'vacuum'
    assert region[1].surface.boundary_type == 'vacuum'
    assert region[2].surface.boundary_type == 'vacuum'
    assert region[3].surface.boundary_type == 'vacuum'
    assert region[4].surface.boundary_type == 'vacuum'
    assert region[5].surface.boundary_type == 'vacuum'

    region = u.bounding_region(bounded_type='sphere', boundary_type='reflective')
    assert isinstance(region, openmc.Region)
    assert isinstance(region, openmc.Halfspace)
    assert region.surface.type == 'sphere'
    assert region.surface.boundary_type == 'reflective'
