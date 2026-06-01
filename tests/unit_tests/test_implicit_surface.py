"""
Tests for openmc.ImplicitSurface and openmc.TPMS.

Run with:  pytest test_implicit_surface.py -v
"""

import warnings

import numpy as np
import pytest

import openmc.implicit as impl
from openmc.implicit import (
    X, Y, Z,
    Add, Mul, Scale, Sin, Cos, Cached,
    ImplicitFunction,
)
from openmc.surface import ImplicitSurface, TPMS


# ==============================================================================
# Construction
# ==============================================================================


def test_valid_minimal():
    """Only the function is required; all other params default correctly."""
    func = X()**2 + Y()**2 + Z()**2
    surf = ImplicitSurface(function=func)
    assert surf.isovalue == 0.0
    assert surf.x0 == 0.0
    assert surf.y0 == 0.0
    assert surf.z0 == 0.0
    assert np.allclose(surf.get_rotation_matrix(), np.eye(3))

def test_valid_full_params():
    """All twelve transform parameters plus isovalue are stored."""
    surf = ImplicitSurface(
        function=X(), isovalue=3.14,
        x0=1., y0=2., z0=3.,
        a=1., b=0., c=0.,
        d=0., e=1., f=0.,
        g=0., h=0., i=1.,
    )
    assert surf.isovalue == pytest.approx(3.14)
    assert surf.x0 == pytest.approx(1.)
    assert surf.y0 == pytest.approx(2.)
    assert surf.z0 == pytest.approx(3.)

def test_isovalue_stored_as_float():
    """isovalue is always stored as float even when passed as int."""
    surf = ImplicitSurface(function=X(), isovalue=5)
    assert isinstance(surf.isovalue, float)

def test_default_rotation_is_identity():
    """Default a-i coefficients encode the 3x3 identity."""
    surf = ImplicitSurface(function=X())
    assert np.allclose(surf.get_rotation_matrix(), np.eye(3))

def test_invalid_boundary_vacuum():
    """Vacuum boundary type raises ValueError."""
    with pytest.raises(ValueError, match="transmission"):
        ImplicitSurface(function=X(), boundary_type='vacuum')

def test_invalid_boundary_reflective():
    """Reflective boundary type raises ValueError."""
    with pytest.raises(ValueError, match="transmission"):
        ImplicitSurface(function=X(), boundary_type='reflective')

def test_invalid_rotation_non_orthogonal():
    """Non-orthogonal a-i matrix raises ValueError."""
    with pytest.raises(ValueError, match="rotation"):
        ImplicitSurface(function=X(),
                        a=2., b=0., c=0.,
                        d=0., e=1., f=0.,
                        g=0., h=0., i=1.)

def test_invalid_rotation_det_minus_one():
    """Orthogonal matrix with det = -1 (reflection) raises ValueError."""
    with pytest.raises(ValueError, match="rotation"):
        ImplicitSurface(function=X(),
                        a=-1., b=0., c=0.,
                        d=0.,  e=1., f=0.,
                        g=0.,  h=0., i=1.)

def test_invalid_function_type_string():
    """Passing a string as function raises TypeError."""
    with pytest.raises(TypeError):
        ImplicitSurface(function="not a function")

def test_invalid_function_type_none():
    """Passing None as function raises TypeError."""
    with pytest.raises(TypeError):
        ImplicitSurface(function=None)

def test_repr_does_not_crash():
    """repr() returns a non-empty string that mentions Function and Isovalue."""
    surf = ImplicitSurface(function=X()**2 + Y()**2 + Z()**2, isovalue=25.)
    s = repr(surf)
    assert len(s) > 0
    assert "Function" in s
    assert "Isovalue" in s


# ==============================================================================
# evaluate
# ==============================================================================

# Sphere: f = x^2 + y^2 + z^2, isovalue = R^2
R = 5.0

@pytest.fixture
def sphere():
    func = X()**2 + Y()**2 + Z()**2
    return ImplicitSurface(function=func, isovalue=R**2)

def test_interior_negative(sphere):
    """Origin is inside the sphere → evaluate < 0."""
    assert sphere.evaluate((0., 0., 0.)) < 0.

def test_exterior_positive(sphere):
    """Point clearly outside sphere → evaluate > 0."""
    assert sphere.evaluate((10., 0., 0.)) > 0.

def test_on_surface_zero(sphere):
    """Point exactly on sphere surface → evaluate ≈ 0."""
    assert sphere.evaluate((R, 0., 0.)) == pytest.approx(0., abs=1e-10)
    assert sphere.evaluate((0., R, 0.)) == pytest.approx(0., abs=1e-10)
    assert sphere.evaluate((0., 0., R)) == pytest.approx(0., abs=1e-10)

def test_isovalue_shifts_zero_level():
    """Larger isovalue makes previously-exterior points interior."""
    func = X()**2 + Y()**2 + Z()**2
    s_small = ImplicitSurface(function=func, isovalue=25.)   # r = 5
    s_large = ImplicitSurface(function=func, isovalue=100.)  # r = 10
    pt = (7., 0., 0.)                                        # between the two radii
    assert s_small.evaluate(pt) > 0.
    assert s_large.evaluate(pt) < 0.

def test_translation_shifts_surface():
    """x0 = 5 moves the sphere centre to (5, 0, 0)."""
    func = X()**2 + Y()**2 + Z()**2
    surf = ImplicitSurface(function=func, isovalue=16., x0=5.)
    assert surf.evaluate((9., 0., 0.)) == pytest.approx(0., abs=1e-10)   # on surface
    assert surf.evaluate((5., 0., 0.)) < 0.                              # centre → inside
    assert surf.evaluate((5., 0., 0.)) == pytest.approx(-16., abs=1e-10) # centre value
    assert surf.evaluate((0., 0., 0.)) > 0.                              # outside

def test_rotation_applied_correctly():
    """90° z-rotation of the x-plane (f = x) makes f = y in world coords.

    Rotation matrix [[0,-1,0],[1,0,0],[0,0,1]] maps r_local.x = -y,
    so evaluate(r) = X().evaluate(r_local) = -y.
    """
    surf = ImplicitSurface(function=X(), isovalue=0.)
    rotated = surf.rotate([0., 0., 90.])
    assert rotated.evaluate((0.,  1., 0.)) == pytest.approx( 1., abs=1e-10)
    assert rotated.evaluate((0., -1., 0.)) == pytest.approx(-1., abs=1e-10)
    assert rotated.evaluate((5.,  0., 0.)) == pytest.approx( 0., abs=1e-10)  # on surface

def test_matrix_rotation_applied_correctly():
    """Same as precedent but with a matrix"""
    surf = ImplicitSurface(function=X(), isovalue=0.)
    rmat = [[0,-1,0],[1,0,0],[0,0,1]]
    rotated = surf.rotate(rmat)
    assert rotated.evaluate((0.,  1., 0.)) == pytest.approx( 1., abs=1e-10)
    assert rotated.evaluate((0., -1., 0.)) == pytest.approx(-1., abs=1e-10)
    assert rotated.evaluate((5.,  0., 0.)) == pytest.approx( 0., abs=1e-10)  # on surface

# ==============================================================================
# translate
# ==============================================================================

def test_updates_origin_components():
    """translate([1, 2, 3]) increments x0, y0, z0 correctly."""
    surf = ImplicitSurface(function=X())
    t = surf.translate([1., 2., 3.])
    assert t.x0 == pytest.approx(1.)
    assert t.y0 == pytest.approx(2.)
    assert t.z0 == pytest.approx(3.)

def test_returns_new_object_by_default():
    """Default inplace=False returns a different Python object."""
    surf = ImplicitSurface(function=X())
    t = surf.translate([1., 0., 0.])
    assert t is not surf

def test_original_unchanged():
    """Original surface is not modified when inplace=False."""
    surf = ImplicitSurface(function=X())
    surf.translate([5., 0., 0.])
    assert surf.x0 == pytest.approx(0.)

def test_inplace_modifies_same_object():
    """inplace=True modifies and returns the same object."""
    surf = ImplicitSurface(function=X())
    result = surf.translate([3., 0., 0.], inplace=True)
    assert result is surf
    assert surf.x0 == pytest.approx(3.)

def test_zero_vector_returns_clone():
    """translate([0,0,0]) returns a clone (new object), not self.

    This differs from the base-class behaviour where the same object is
    returned for zero translations.
    """
    surf = ImplicitSurface(function=X())
    result = surf.translate([0., 0., 0.])
    assert result is not surf

def test_preserves_function_reference():
    """translate does not copy or modify the ImplicitFunction DAG."""
    func = X()**2 + Y()**2 + Z()**2
    surf = ImplicitSurface(function=func)
    t = surf.translate([1., 0., 0.])
    assert t.function is func

def test_evaluate_consistency():
    """After translate([5,0,0]), sphere centre moves accordingly."""
    func = X()**2 + Y()**2 + Z()**2
    surf = ImplicitSurface(function=func, isovalue=16.)
    t = surf.translate([5., 0., 0.])
    assert t.evaluate((9., 0., 0.)) == pytest.approx(0., abs=1e-10)    # on surface
    assert t.evaluate((5.,  0., 0.)) == pytest.approx(-16., abs=1e-10) # centre
    assert t.evaluate((0.,  0., 0.)) > 0.                              # outside


# ==============================================================================
# rotate
# ==============================================================================

def test_euler_angles_update_rotation_matrix():
    """rotate([0, 0, 90]) applies Rz(90°) to the rotation matrix."""
    surf = ImplicitSurface(function=X())
    rotated = surf.rotate([0., 0., 90.])
    R90z = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    expected = R90z.T
    assert np.allclose(rotated.get_rotation_matrix(), expected, atol=1e-10)

def test_rotation_matrix_passed_directly():
    """A 3x3 ndarray can be passed directly as the rotation."""
    R90z = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    surf = ImplicitSurface(function=X())
    rotated = surf.rotate(R90z)
    assert np.allclose(rotated.get_rotation_matrix(), R90z.T, atol=1e-10)

@pytest.mark.parametrize("angles", [
    [30., 0., 0.], [0., 45., 0.], [0., 0., 90.],
    [30., 45., 60.], [15., -30., 75.],
])
def test_result_is_always_valid_rotation(angles):
    """Rotated surface always has an orthogonal matrix with det = +1."""
    surf = ImplicitSurface(function=X())
    rotated = surf.rotate(angles)
    R = rotated.get_rotation_matrix()
    assert np.allclose(R @ R.T, np.eye(3), atol=1e-10), \
        f"Not orthogonal for angles {angles}"
    assert np.isclose(np.linalg.det(R), 1.0, atol=1e-10), \
        f"det ≠ 1 for angles {angles}"

def test_preserves_function_reference():
    """rotate does not copy or modify the ImplicitFunction DAG."""
    func = X()
    surf = ImplicitSurface(function=func)
    assert surf.rotate([0., 0., 45.]).function is func

def test_pivot_translates_origin():
    """Rotating about a non-origin pivot moves the surface origin correctly.

    Sphere centred at (5,0,0); 90° z-rotation about origin → centre at (0,5,0).
    """
    func = X()**2 + Y()**2 + Z()**2
    surf = ImplicitSurface(function=func, isovalue=4., x0=5.)
    rotated = surf.rotate([0., 0., 90.], pivot=(0., 0., 0.))
    assert rotated.x0 == pytest.approx(0., abs=1e-10)
    assert rotated.y0 == pytest.approx(5., abs=1e-10)
    assert rotated.z0 == pytest.approx(0., abs=1e-10)


# ==============================================================================
# is_equal and clone
# ==============================================================================


def test_equal_to_itself():
    func = X()**2 + Y()**2 + Z()**2
    surf = ImplicitSurface(function=func, isovalue=25.)
    assert surf.is_equal(surf)

def test_equal_same_function_object_same_params():
    """Two surfaces sharing the same Python function object and params."""
    func = X()
    s1 = ImplicitSurface(function=func, isovalue=1.)
    s2 = ImplicitSurface(function=func, isovalue=1.)
    assert s1.is_equal(s2)

def test_not_equal_different_isovalue():
    func = X()
    s1 = ImplicitSurface(function=func, isovalue=1.)
    s2 = ImplicitSurface(function=func, isovalue=2.)
    assert not s1.is_equal(s2)

def test_not_equal_different_translation():
    func = X()
    s1 = ImplicitSurface(function=func, x0=0.)
    s2 = ImplicitSurface(function=func, x0=1.)
    assert not s1.is_equal(s2)

def test_not_equal_different_function_objects():
    """is_equal uses Python object identity for the function — two
    equivalent but distinct objects are NOT equal."""
    s1 = ImplicitSurface(function=X())
    s2 = ImplicitSurface(function=X())
    assert not s1.is_equal(s2)

def test_clone_has_new_id():
    surf = ImplicitSurface(function=X())
    clone = surf.clone()
    assert clone.id != surf.id

def test_clone_shares_function_reference():
    """clone preserves the Python function object, not a deep copy."""
    func = X()**2 + Y()**2 + Z()**2
    surf = ImplicitSurface(function=func, isovalue=25.)
    assert surf.clone().function is func

def test_clone_independent_coefficients():
    """Translating a clone does not affect the original."""
    surf = ImplicitSurface(function=X())
    clone = surf.clone()
    clone.translate([5., 0., 0.], inplace=True)
    assert surf.x0 == pytest.approx(0.)


# ==============================================================================
# bounding_box
# ==============================================================================

def test_always_infinite():
    """bounding_box returns an infinite box on both sides — no analytical
    bound is available without a user-supplied bbox."""
    surf = ImplicitSurface(function=X()**2 + Y()**2 + Z()**2, isovalue=25.)
    bb_pos = surf.bounding_box('+')
    bb_neg = surf.bounding_box('-')
    assert np.all(bb_pos.lower_left == -np.inf)
    assert np.all(bb_neg.lower_left == -np.inf)
    assert np.all(bb_pos.upper_right == +np.inf)
    assert np.all(bb_neg.upper_right == +np.inf)


# ==============================================================================
# XML round-trip
# ==============================================================================

def test_roundtrip_basic():
    """to_xml_element + from_xml_element produces a surface that evaluates
    identically."""
    func = X()**2 + Y()**2 + Z()**2
    surf = ImplicitSurface(function=func, isovalue=25.)
    elem = surf.to_xml_element()
    restored = ImplicitSurface.from_xml_element(elem)
    for pt in [(3., 4., 0.), (0., 0., 0.), (5., 0., 0.)]:
        assert restored.evaluate(pt) == pytest.approx(surf.evaluate(pt), abs=1e-10)

def test_roundtrip_isovalue_preserved():
    """isovalue survives the XML round-trip exactly."""
    surf = ImplicitSurface(function=X(), isovalue=3.14)
    restored = ImplicitSurface.from_xml_element(surf.to_xml_element())
    assert restored.isovalue == pytest.approx(3.14)

def test_roundtrip_transform_params_preserved():
    """All 12 transform coefficients survive the round-trip."""
    surf = ImplicitSurface(function=X(), isovalue=0.,
                            x0=1., y0=2., z0=3.,
                            a=1., b=0., c=0.,
                            d=0., e=1., f=0.,
                            g=0., h=0., i=1.)
    restored = ImplicitSurface.from_xml_element(surf.to_xml_element())
    assert np.allclose(restored._get_base_coeffs(),
                        surf._get_base_coeffs(), atol=1e-10)

def test_roundtrip_with_cached_nodes():
    """Surface with Cached nodes evaluates identically after round-trip."""
    cx = Cached(2. * X())
    surf = ImplicitSurface(function=Sin(cx) * Cos(cx), isovalue=0.)
    restored = ImplicitSurface.from_xml_element(surf.to_xml_element())
    for pt in [(0., 0., 0.), (0.5, 0., 0.), (1., 0., 0.)]:
        assert restored.evaluate(pt) == pytest.approx(surf.evaluate(pt), abs=1e-10)

def test_roundtrip_produces_implicit_surface_not_tpms():
    """from_xml_element always produces ImplicitSurface (not TPMS) because
    TPMS is not a direct subclass of Surface and is absent from
    _SURFACE_CLASSES."""
    surf = TPMS.from_pitch_isovalue("gyroid", 1.0, 0.0)
    restored = ImplicitSurface.from_xml_element(surf.to_xml_element())
    assert type(restored) is ImplicitSurface

def test_single_use_cached_warns():
    """Cached node used only once triggers a UserWarning during serialisation."""
    surf = ImplicitSurface(function=Sin(Cached(X())))
    with pytest.warns(UserWarning, match="only used once"):
        surf.to_xml_element()

def test_reused_cached_no_warning():
    """Cached node used twice (shared) triggers no UserWarning."""
    cx = Cached(X())
    surf = ImplicitSurface(function=Sin(cx) + Cos(cx))
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        surf.to_xml_element()   # must not raise


# ==============================================================================
# TPMS
# ==============================================================================

# ── construction ──────────────────────────────────────────────────────────

@pytest.mark.parametrize("name", [
    "primitive", "schwarz_p",
    "gyroid",    "schoen-g",
    "diamond",   "schwarz_d",
])
def test_valid_names_produce_implicit_surface(name):
    assert isinstance(TPMS.from_pitch_isovalue(name, 1.0, 0.0), ImplicitSurface)

def test_unknown_name_raises():
    """Unknown TPMS name raises NotImplementedError."""
    with pytest.raises(NotImplementedError):
        TPMS.from_pitch_isovalue("screugneugneu", 1.0, 0.0)

def test_isovalue_stored_correctly():
    surf = TPMS.from_pitch_isovalue("primitive", 1.0, 0.5)
    assert surf.isovalue == pytest.approx(0.5)

# ── primitive (Schwartz-P) ────────────────────────────────────────────────

def test_primitive_at_origin():
    """cos(0)+cos(0)+cos(0) = 3.0, isovalue=0 → evaluate = 3."""
    surf = TPMS.from_pitch_isovalue("primitive", 1.0, 0.0)
    assert surf.evaluate((0., 0., 0.)) == pytest.approx(3.0, abs=1e-10)

def test_primitive_periodicity():
    """Primitive is periodic with pitch L along each axis."""
    L = 1.0
    surf = TPMS.from_pitch_isovalue("primitive", L, 0.0)
    for pt in [(0.1, 0.2, 0.3), (0.3, 0.4, 0.1)]:
        v_orig    = surf.evaluate(pt)
        v_shift_x = surf.evaluate((pt[0] + L, pt[1],     pt[2]))
        v_shift_y = surf.evaluate((pt[0],     pt[1] + L, pt[2]))
        v_shift_z = surf.evaluate((pt[0],     pt[1],     pt[2] + L))
        assert v_orig == pytest.approx(v_shift_x, abs=1e-10)
        assert v_orig == pytest.approx(v_shift_y, abs=1e-10)
        assert v_orig == pytest.approx(v_shift_z, abs=1e-10)

def test_pitch_scales_spatial_frequency():
    """Doubling the pitch halves the spatial frequency:
    evaluate_2L(x, 0, 0) == evaluate_L(x/2, 0, 0)."""
    s1 = TPMS.from_pitch_isovalue("primitive", 1.0, 0.0)
    s2 = TPMS.from_pitch_isovalue("primitive", 2.0, 0.0)
    for x in [0.1, 0.2, 0.35]:
        assert s2.evaluate((x, 0., 0.)) == pytest.approx(
            s1.evaluate((x / 2., 0., 0.)), abs=1e-10)

def test_isovalue_shifts_evaluate():
    """Increasing isovalue by Δ decreases evaluate by Δ at every point."""
    s0 = TPMS.from_pitch_isovalue("primitive", 1.0, 0.0)
    s1 = TPMS.from_pitch_isovalue("primitive", 1.0, 1.0)
    for pt in [(0., 0., 0.), (0.1, 0.2, 0.3)]:
        assert s0.evaluate(pt) - s1.evaluate(pt) == pytest.approx(1.0, abs=1e-10)

# ── gyroid ────────────────────────────────────────────────────────────────

def test_gyroid_at_origin():
    """Gyroid f = sin·cos + sin·cos + sin·cos = 0 at the origin."""
    surf = TPMS.from_pitch_isovalue("gyroid", 1.0, 0.0)
    assert surf.evaluate((0., 0., 0.)) == pytest.approx(0., abs=1e-10)

def test_gyroid_at_quarter_pitch():
    """Gyroid at (L/4, 0, 0) should be 1.0.

    Analytical:
        x_scaled = 2π*(L/4)/L = π/2,  y_scaled = 0,  z_scaled = 0
        sin(π/2)*cos(0) + sin(0)*cos(π/2) + sin(0)*cos(0) = 1 + 0 + 0 = 1. 
    """
    L = 1.0
    surf = TPMS.from_pitch_isovalue("gyroid", L, 0.0)
    assert surf.evaluate((L / 4., 0., 0.)) == pytest.approx(1.0, abs=1e-10)

def test_gyroid_periodicity():
    """Gyroid is periodic with pitch L along x."""
    L = 1.0
    surf = TPMS.from_pitch_isovalue("gyroid", L, 0.0)
    for pt in [(0.1, 0.2, 0.3), (0.4, 0.1, 0.5)]:
        assert surf.evaluate(pt) == pytest.approx(
            surf.evaluate((pt[0] + L, pt[1], pt[2])), abs=1e-10)

# ── diamond —──────────────────────────────────────────────────────────────

def test_diamond_at_origin():
    """Diamond f = sin(x)*cos(y-z) + sin(y+z)*cos(x) = 0 at the origin."""
    surf = TPMS.from_pitch_isovalue("diamond", 1.0, 0.0)
    assert surf.evaluate((0., 0., 0.)) == pytest.approx(0., abs=1e-10)

def test_diamond_at_eighth_pitch():
    """Diamond at (L/8, 0, 0) should equal sin(π/4) = √2/2.

    Analytical:
        x_scaled = 2π*(L/8)/L = π/4,  y_scaled = 0,  z_scaled = 0
        sin(π/4)*cos(0-0) + sin(0+0)*cos(π/4) = sin(π/4) + 0 = √2/2.
    """
    L = 1.0
    surf = TPMS.from_pitch_isovalue("diamond", L, 0.0)
    expected = np.sin(np.pi / 4.)   # ≈ 0.7071
    assert surf.evaluate((L / 8., 0., 0.)) == pytest.approx(expected, abs=1e-10)

def test_diamond_periodicity():
    """Diamond is periodic with pitch L along x."""
    L = 1.0
    surf = TPMS.from_pitch_isovalue("diamond", L, 0.0)
    for pt in [(0.1, 0.2, 0.3), (0.4, 0.1, 0.5)]:
        assert surf.evaluate(pt) == pytest.approx(
            surf.evaluate((pt[0] + L, pt[1], pt[2])), abs=1e-10)