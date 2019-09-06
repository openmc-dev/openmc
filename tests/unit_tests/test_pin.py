"""
Tests for constructing Pin universes
"""

import numpy as np
import pytest

import openmc
from openmc.model import pin


def get_pin_radii(pin_univ):
    """Return a sorted list of all radii from pin"""
    rads = set()

    for cell in pin_univ.get_all_cells().values():
        surfs = cell.region.get_surfaces().values()
        rads.update(set(s.r for s in surfs))

    return list(sorted(rads))


@pytest.fixture
def pin_mats():
    fuel = openmc.Material(name="UO2")
    fuel.volume = 100
    clad = openmc.Material(name="zirc")
    clad.volume = 100
    water = openmc.Material(name="water")
    return fuel, clad, water


@pytest.fixture
def good_radii():
    return (0.4, 0.42)


def test_failure(pin_mats, good_radii):
    """Check for various failure modes"""
    good_surfaces = [openmc.ZCylinder(r=r) for r in good_radii]
    # Bad material type
    with pytest.raises(TypeError):
        pin(good_surfaces, [mat.name for mat in pin_mats])

    # Incorrect lengths
    with pytest.raises(ValueError, match="length"):
        pin(good_surfaces[:len(pin_mats) - 2], pin_mats)

    # Non-positive radii
    rad = [openmc.ZCylinder(r=-0.1)] + good_surfaces[1:]
    with pytest.raises(ValueError, match="index 0"):
        pin(rad, pin_mats)

    # Non-increasing radii
    surfs = tuple(reversed(good_surfaces))
    with pytest.raises(ValueError, match="index 1"):
        pin(surfs, pin_mats)

    # Bad orientation
    surfs = [openmc.XCylinder(r=good_surfaces[0].r)] + good_surfaces[1:]
    with pytest.raises(TypeError, match="surfaces"):
        pin(surfs, pin_mats)

    # Passing cells argument
    with pytest.raises(SyntaxError, match="Cells"):
        pin(surfs, pin_mats, cells=[])


def test_pins_of_universes(pin_mats, good_radii):
    """Build a pin with a Universe in one ring"""
    u1 = openmc.Universe(cells=[openmc.Cell(fill=pin_mats[1])])
    new_items = pin_mats[:1] + (u1, ) + pin_mats[2:]
    new_pin = pin(
        [openmc.ZCylinder(r=r) for r in good_radii], new_items,
        subdivisions={0: 2}, divide_vols=True)
    assert len(new_pin.cells) == len(pin_mats) + 1


@pytest.mark.parametrize(
    "surf_type", [openmc.ZCylinder, openmc.XCylinder, openmc.YCylinder])
def test_subdivide(pin_mats, good_radii, surf_type):
    """Test the subdivision with various orientations"""
    surfs = [surf_type(r=r) for r in good_radii]
    fresh = pin(surfs, pin_mats, name="fresh pin")
    assert len(fresh.cells) == len(pin_mats)
    assert fresh.name == "fresh pin"

    # subdivide inner region
    N = 5
    div0 = pin(surfs, pin_mats, {0: N})
    assert len(div0.cells) == len(pin_mats) + N - 1

    # Check volume of fuel material
    for mid, mat in div0.get_all_materials().items():
        if mat.name == "UO2":
            assert mat.volume == pytest.approx(100 / N)

    # check volumes of new rings
    radii = get_pin_radii(div0)
    bounds = [0] + radii[:N]
    sqrs = np.square(bounds)
    assert np.all(sqrs[1:] - sqrs[:-1] == pytest.approx(good_radii[0] ** 2 / N))

    # subdivide non-inner most region
    new_pin = pin(surfs, pin_mats, {1: N})
    assert len(new_pin.cells) == len(pin_mats) + N - 1

    # Check volume of clad material
    for mid, mat in div0.get_all_materials().items():
        if mat.name == "zirc":
            assert mat.volume == pytest.approx(100 / N)

    # check volumes of new rings
    radii = get_pin_radii(new_pin)
    sqrs = np.square(radii[:N + 1])
    assert np.all(sqrs[1:] - sqrs[:-1] == pytest.approx(
        (good_radii[1] ** 2 - good_radii[0] ** 2) / N))
