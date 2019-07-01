"""
Tests for constructing Pin universes
"""

import numpy
import pytest

import openmc
from openmc.model import Pin


@pytest.fixture
def pin_mats():
    fuel = openmc.Material(name="UO2")
    clad = openmc.Material(name="zirc")
    water = openmc.Material(name="water")
    return fuel, clad, water


@pytest.fixture
def good_radii():
    return (0.4, 0.42)


def test_failure(pin_mats, good_radii):
    """Check for various failure modes"""
    # Bad material type
    with pytest.raises(TypeError):
        Pin.from_radii(good_radii, [mat.name for mat in pin_mats])

    # Incorrect lengths
    with pytest.raises(ValueError, match="length") as exec_info:
        Pin.from_radii(good_radii[: len(pin_mats) - 2], pin_mats)

    # Non-positive radii
    rad = (-0.1,) + good_radii[1:]
    with pytest.raises(ValueError, match="index 0") as exec_info:
        Pin.from_radii(rad, pin_mats)

    # Non-increasing radii
    rad = tuple(reversed(good_radii))
    with pytest.raises(ValueError, match="index 1") as exec_info:
        Pin.from_radii(rad, pin_mats)

    # Bad orientation
    with pytest.raises(ValueError, match="Orientation") as exec_info:
        Pin.from_radii(good_radii, pin_mats, orientation="fail")


def test_from_radii(pin_mats, good_radii):
    name = "test pin"
    p = Pin.from_radii(good_radii, pin_mats, name=name)
    assert len(p.cells) == len(pin_mats)
    assert p.name == name
    assert p.radii == good_radii


def test_subdivide(pin_mats, good_radii):
    surfs = [openmc.ZCylinder(r=r) for r in good_radii]
    pin = Pin(surfs, pin_mats)
    assert pin.radii == good_radii
    assert len(pin.cells) == len(pin_mats)

    # subdivide inner region
    N = 5
    pin.subdivide_ring(0, N)
    assert len(pin.radii) == len(good_radii) + N - 1
    assert len(pin.cells) == len(pin_mats) + N - 1
    # check volumes of new rings
    bounds = (0,) + pin.radii[:N]
    sqrs = numpy.square(bounds)
    assert sqrs[1:] - sqrs[:-1] == pytest.approx(good_radii[0] ** 2 / N)

    # subdivide non-inner most region
    new_pin = Pin.from_radii(good_radii, pin_mats)
    new_pin.subdivide_ring(1, N)
    assert len(new_pin.radii) == len(good_radii) + N - 1
    assert len(new_pin.cells) == len(pin_mats) + N - 1
    # check volumes of new rings
    bounds = new_pin.radii[:N + 1]
    sqrs = numpy.square(bounds)
    assert sqrs[1:] - sqrs[:-1] == pytest.approx(
        (good_radii[1] ** 2 - good_radii[0] ** 2) / N)
