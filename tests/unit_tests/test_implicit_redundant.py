"""
Tests for Geometry.remove_redundant_surfaces with ImplicitSurface.

The fix ensures that is_equal() is used within each coefficient-key bucket,
so that ImplicitSurface objects with identical transforms but different
isovalue or function are NOT incorrectly declared redundant.
"""

import numpy as np
import pytest

import openmc
from openmc.surface import ImplicitSurface
from openmc.implicit import X, Y, Z, Sin, Cos, Cached


# ==============================================================================
# Helpers
# ==============================================================================

def make_geometry(*cells):
    """Wrap cells in a minimal geometry."""
    universe = openmc.Universe(cells=list(cells))
    return openmc.Geometry(universe)


def make_cell(region, fill=None):
    mat = openmc.Material()
    mat.add_nuclide('U235', 1.0)
    mat.set_density('g/cm3', 16.0)
    return openmc.Cell(region=region, fill=fill or mat)


# ==============================================================================
# Standard surfaces — existing behaviour preserved
# ==============================================================================

def test_identical_spheres_are_redundant():
    """Two spheres with the same radius are redundant."""
    s1 = openmc.Sphere(r=5.0, boundary_type='vacuum')
    s2 = openmc.Sphere(r=5.0, boundary_type='vacuum')
    c1 = make_cell(-s1)
    c2 = make_cell(-s2)
    geom = make_geometry(c1, c2)
    redundant = geom.remove_redundant_surfaces()
    assert len(redundant) == 1
    assert s2.id in redundant
    assert redundant[s2.id] is s1

def test_different_spheres_are_not_redundant():
    """Two spheres with different radii are not redundant."""
    s1 = openmc.Sphere(r=5.0, boundary_type='vacuum')
    s2 = openmc.Sphere(r=6.0, boundary_type='vacuum')
    c1 = make_cell(-s1)
    c2 = make_cell(-s2)
    geom = make_geometry(c1, c2)
    redundant = geom.remove_redundant_surfaces()
    assert len(redundant) == 0

def test_no_redundant_surfaces():
    """Geometry with all unique surfaces returns empty dict."""
    s = openmc.Sphere(r=5.0, boundary_type='vacuum')
    geom = make_geometry(make_cell(-s))
    assert geom.remove_redundant_surfaces() == {}


# ==============================================================================
# ImplicitSurface
# ==============================================================================


def test_same_function_same_isovalue_is_redundant():
    """Two implicit surfaces sharing the same function object and isovalue
    are genuinely redundant."""
    func = X()**2 + Y()**2 + Z()**2
    outer = openmc.Sphere(r=11.0, boundary_type='vacuum')
    s1 = ImplicitSurface(function=func, isovalue=25.)
    s2 = ImplicitSurface(function=func, isovalue=25.)
    c1 = make_cell(-s1 & -outer)
    c2 = make_cell(-s2 & -outer)
    geom = make_geometry(c1, c2)
    redundant = geom.remove_redundant_surfaces()
    assert len(redundant) == 1
    assert s2.id in redundant
    assert redundant[s2.id] is s1

def test_same_function_different_isovalue_is_not_redundant():
    """Same function, different isovalue — different surfaces.
    This was the bug: both surfaces had identical _coefficients keys
    but different isovalues, yet were incorrectly declared redundant."""
    func = X()**2 + Y()**2 + Z()**2
    outer = openmc.Sphere(r=11.0, boundary_type='vacuum')
    s1 = ImplicitSurface(function=func, isovalue=25.)  # r = 5
    s2 = ImplicitSurface(function=func, isovalue=100.) # r = 10
    c1 = make_cell(-s1 & -outer)
    c2 = make_cell(-s2 & -outer)
    geom = make_geometry(c1, c2)
    redundant = geom.remove_redundant_surfaces()
    assert len(redundant) == 0

def test_different_function_objects_not_redundant():
    """Equivalent but distinct function objects are not redundant —
    is_equal uses Python object identity for the function."""
    outer = openmc.Sphere(r=11.0, boundary_type='vacuum')
    s1 = ImplicitSurface(function=X()**2 + Y()**2 + Z()**2, isovalue=25.)
    s2 = ImplicitSurface(function=X()**2 + Y()**2 + Z()**2, isovalue=25.)
    c1 = make_cell(-s1 & -outer)
    c2 = make_cell(-s2 & -outer)
    geom = make_geometry(c1, c2)
    redundant = geom.remove_redundant_surfaces()
    assert len(redundant) == 0

def test_different_transforms_not_redundant():
    """Same function and isovalue but different x0 — not redundant."""
    func = X()**2 + Y()**2 + Z()**2
    outer = openmc.Sphere(r=11.0, boundary_type='vacuum')
    s1 = ImplicitSurface(function=func, isovalue=25., x0=0.)
    s2 = ImplicitSurface(function=func, isovalue=25., x0=1.)
    c1 = make_cell(-s1 & -outer)
    c2 = make_cell(-s2 & -outer)
    geom = make_geometry(c1, c2)
    redundant = geom.remove_redundant_surfaces()
    assert len(redundant) == 0

def test_tpms_same_function_object_is_redundant():
    """Two TPMS surfaces sharing the same underlying function object
    and isovalue are correctly identified as redundant."""
    func = Sin(X()) * Cos(Z()) + Sin(Y()) * Cos(X()) + Sin(Z()) * Cos(Y())
    outer = openmc.Sphere(r=11.0, boundary_type='vacuum')
    s1 = ImplicitSurface(function=func, isovalue=0.)
    s2 = ImplicitSurface(function=func, isovalue=0.)
    c1 = make_cell(-s1 & -outer)
    c2 = make_cell(-s2 & -outer)
    geom = make_geometry(c1, c2)
    redundant = geom.remove_redundant_surfaces()
    assert len(redundant) == 1
    assert redundant[s2.id] is s1


# ==============================================================================
# Mixed geometry
# ==============================================================================

def test_standard_and_implicit_both_handled():
    """In a geometry with both standard and implicit surfaces,
    redundancy detection works correctly for both types independently."""
    # Standard redundant pair
    sphere_a1 = openmc.Sphere(r=5.0)
    sphere_a2 = openmc.Sphere(r=5.0)

    # Standard non-redundant pair
    sphere_b = openmc.Sphere(r=8.0)

    # Implicit non-redundant pair (different isovalue)
    func = X()**2 + Y()**2 + Z()**2
    outer = openmc.Sphere(r=11.0, boundary_type='vacuum')
    impl1 = ImplicitSurface(function=func, isovalue=25.)
    impl2 = ImplicitSurface(function=func, isovalue=64.)

    c1 = make_cell(-sphere_a1 & -outer)
    c2 = make_cell(-sphere_a2 & -outer)
    c3 = make_cell(-sphere_b  & -outer)
    c4 = make_cell(-impl1     & -outer)
    c5 = make_cell(-impl2     & -outer)
    geom = make_geometry(c1, c2, c3, c4, c5)

    redundant = geom.remove_redundant_surfaces()

    # Only sphere_a2 is redundant
    assert len(redundant) == 1
    assert sphere_a2.id in redundant
    assert impl1.id not in redundant
    assert impl2.id not in redundant