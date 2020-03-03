import numpy as np
import pytest
import openmc

from tests.unit_tests import assert_unbounded


@pytest.fixture
def reset():
    openmc.reset_auto_ids()


def test_union(reset):
    s1 = openmc.XPlane(x0=5, surface_id=1)
    s2 = openmc.XPlane(x0=-5, surface_id=2)
    region = +s1 | -s2
    assert isinstance(region, openmc.Union)

    # Check bounding box
    assert_unbounded(region)

    # __contains__
    assert (6, 0, 0) in region
    assert (-6, 0, 0) in region
    assert (0, 0, 0) not in region

    # string representation
    assert str(region) == '(1 | -2)'

    # Combining region with intersection
    s3 = openmc.YPlane(surface_id=3)
    reg2 = region & +s3
    assert (6, 1, 0) in reg2
    assert (6, -1, 0) not in reg2
    assert str(reg2) == '((1 | -2) 3)'

    # translate method
    regt = region.translate((2.0, 0.0, 0.0))
    assert (-4, 0, 0) in regt
    assert (6, 0, 0) not in regt
    assert (8, 0, 0) in regt

    # rotate method
    regr = region.rotate((0., 90., 0.))
    assert (-4, 0, 0) not in regr
    assert (0, 0, 6) in regr
    assert (0, 0, -6) in regr
    assert (0, 0, 3) not in regr


def test_intersection(reset):
    s1 = openmc.XPlane(x0=5, surface_id=1)
    s2 = openmc.XPlane(x0=-5, surface_id=2)
    region = -s1 & +s2
    assert isinstance(region, openmc.Intersection)

    # Check bounding box
    ll, ur = region.bounding_box
    assert ll == pytest.approx((-5, -np.inf, -np.inf))
    assert ur == pytest.approx((5, np.inf, np.inf))

    # __contains__
    assert (6, 0, 0) not in region
    assert (-6, 0, 0) not in region
    assert (0, 0, 0) in region

    # string representation
    assert str(region) == '(-1 2)'

    # Combining region with union
    s3 = openmc.YPlane(surface_id=3)
    reg2 = region | +s3
    assert (-6, 2, 0) in reg2
    assert (-6, -2, 0) not in reg2
    assert str(reg2) == '((-1 2) | 3)'

    # translate method
    regt = region.translate((2.0, 0.0, 0.0))
    assert (-4, 0, 0) not in regt
    assert (6, 0, 0) in regt
    assert (8, 0, 0) not in regt

    # rotate method
    regr = region.rotate((0., 90., 0.))
    assert (-4, 0, 0) in regr
    assert (0, 0, 6) not in regr
    assert (0, 0, -6) not in regr
    assert (0, 0, 3) in regr


def test_complement(reset):
    zcyl = openmc.ZCylinder(r=1., surface_id=1)
    z0 = openmc.ZPlane(-5., surface_id=2)
    z1 = openmc.ZPlane(5., surface_id=3)
    outside = +zcyl | -z0 | +z1
    inside = ~outside
    outside_equiv = ~(-zcyl & +z0 & -z1)
    inside_equiv = ~outside_equiv

    # Check bounding box
    for region in (inside, inside_equiv):
        ll, ur = region.bounding_box
        assert ll == pytest.approx((-1., -1., -5.))
        assert ur == pytest.approx((1., 1., 5.))
    assert_unbounded(outside)
    assert_unbounded(outside_equiv)

    # string represention
    assert str(inside) == '~(1 | -2 | 3)'

    # evaluate method
    assert (0, 0, 0) in inside
    assert (0, 0, 0) not in outside
    assert (0, 0, 6) not in inside
    assert (0, 0, 6) in outside

    # translate method
    inside_t = inside.translate((1.0, 1.0, 1.0))
    ll, ur = inside_t.bounding_box
    assert ll == pytest.approx((0., 0., -4.))
    assert ur == pytest.approx((2., 2., 6.))

    # rotate method
    inside_r = inside.rotate((90., 0., 0.))
    ll, ur = inside_r.bounding_box
    assert (.5, 2, 0) in inside_r
    assert (0, 0, 6) not in inside_r
    assert ll == pytest.approx((-1., -5., -1.))
    assert ur == pytest.approx((1., 5., 1.))


def test_get_surfaces():
    s1 = openmc.XPlane()
    s2 = openmc.YPlane()
    s3 = openmc.ZPlane()
    region = (+s1 & -s2) | +s3

    # Make sure get_surfaces() returns all surfaces
    surfs = set(region.get_surfaces().values())
    assert not (surfs ^ {s1, s2, s3})

    inverse = ~region
    surfs = set(inverse.get_surfaces().values())
    assert not (surfs ^ {s1, s2, s3})


def test_extend_clone():
    s1 = openmc.XPlane()
    s2 = openmc.YPlane()
    s3 = openmc.ZPlane()
    s4 = openmc.ZCylinder()

    # extend intersection
    r1 = +s1 & -s2
    r1 &= +s3 & -s4
    assert r1[:] == [+s1, -s2, +s3, -s4]

    # extend union
    r2 = +s1 | -s2
    r2 |= +s3 | -s4
    assert r2[:] == [+s1, -s2, +s3, -s4]

    # clone methods
    r3 = r1.clone()
    assert len(r3) == len(r1)
    r4 = r2.clone()
    assert len(r4) == len(r2)

    r5 = ~r1
    r6 = r5.clone()


def test_from_expression(reset):
    # Create surface dictionary
    s1 = openmc.ZCylinder(surface_id=1)
    s2 = openmc.ZPlane(-10., surface_id=2)
    s3 = openmc.ZPlane(10., surface_id=3)
    surfs = {1: s1, 2: s2, 3: s3}

    r = openmc.Region.from_expression('-1 2 -3', surfs)
    assert isinstance(r, openmc.Intersection)
    assert r[:] == [-s1, +s2, -s3]

    r = openmc.Region.from_expression('+1 | -2 | +3', surfs)
    assert isinstance(r, openmc.Union)
    assert r[:] == [+s1, -s2, +s3]

    r = openmc.Region.from_expression('~(-1)', surfs)
    assert r == +s1

    # Since & has higher precendence than |, the resulting region should be an
    # instance of Union
    r = openmc.Region.from_expression('1 -2 | 3', surfs)
    assert isinstance(r, openmc.Union)
    assert isinstance(r[0], openmc.Intersection)
    assert r[0][:] == [+s1, -s2]

    # ...but not if we use parentheses
    r = openmc.Region.from_expression('1 (-2 | 3)', surfs)
    assert isinstance(r, openmc.Intersection)
    assert isinstance(r[1], openmc.Union)
    assert r[1][:] == [-s2, +s3]
