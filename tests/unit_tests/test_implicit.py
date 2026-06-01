"""
Temporary pytest suite for openmc/implicit_function.py
"""
import pytest
import numpy as np
import lxml.etree as ET

from openmc.implicit import (
    ImplicitFunction,
    X, Y, Z, Constant,
    Add, Sub, Mul, Div, Scale, Neg, Pow,
    Sin, Cos, Sqrt, Exp, Log, Abs,
    Cached,
)

# Canonical test point used throughout
PT = (1.0, 2.0, 3.0)


# ── helpers ────────────────────────────────────────────────────────────────

def xml_roundtrip(func: ImplicitFunction, point=PT) -> float:
    """Serialise → deserialise → evaluate.  Must match the original."""
    elem     = func.to_xml_element([])
    restored = ImplicitFunction.from_xml_element(elem)
    return restored.evaluate(point)


# ── terminals ──────────────────────────────────────────────────────────────

def test_x():            assert X().evaluate(PT)          == 1.0
def test_y():            assert Y().evaluate(PT)          == 2.0
def test_z():            assert Z().evaluate(PT)          == 3.0
def test_constant():     assert Constant(5.0).evaluate(PT) == 5.0
def test_constant_int(): assert Constant(3).evaluate(PT)  == 3.0


# ── arithmetic ─────────────────────────────────────────────────────────────

def test_add():   assert Add(X(), Y()).evaluate(PT)       == 3.0
def test_sub():   assert Sub(X(), Y()).evaluate(PT)       == -1.0
def test_mul():   assert Mul(X(), Y()).evaluate(PT)       == 2.0
def test_div():   assert Div(Y(), X()).evaluate(PT)       == 2.0
def test_scale(): assert Scale(X(), 3.0).evaluate(PT)    == 3.0
def test_neg():   assert Neg(X()).evaluate(PT)            == -1.0
def test_pow_int():   assert Pow(Y(), 2).evaluate(PT)    == 4.0


# ── transcendentals ────────────────────────────────────────────────────────

def test_sin():  assert Sin(Constant(np.pi / 2)).evaluate(PT) == pytest.approx(1.0)
def test_cos():  assert Cos(Constant(0.0)).evaluate(PT)       == pytest.approx(1.0)
def test_sqrt(): assert Sqrt(Constant(4.0)).evaluate(PT)      == pytest.approx(2.0)
def test_exp():  assert Exp(Constant(0.0)).evaluate(PT)       == pytest.approx(1.0)
def test_log():  assert Log(Constant(np.e)).evaluate(PT)      == pytest.approx(1.0)
def test_abs_pos(): assert Abs(X()).evaluate(PT)              == 1.0
def test_abs_neg(): assert Abs(Neg(X())).evaluate(PT)        == 1.0


# ── domain guards ──────────────────────────────────────────────────────────

def test_div_zero():
    with pytest.raises(ValueError, match="zero"):
        Div(X(), Constant(0.0)).evaluate(PT)

def test_sqrt_negative():
    with pytest.raises(ValueError, match="negative"):
        Sqrt(Constant(-1.0)).evaluate(PT)

def test_log_negative():
    with pytest.raises(ValueError, match="negative"):
        Log(Constant(-1.0)).evaluate(PT)

def test_pow_float():
    with pytest.raises(TypeError):
        X() ** 2.0

def test_pow_non_scalar_via_operator():
    with pytest.raises(TypeError):
        X() ** Y()


# ── operator overloading ───────────────────────────────────────────────────

def test_add():          assert (X() + Y()).evaluate(PT)    == 3.0
def test_radd():         assert (1.0 + X()).evaluate(PT)    == 2.0
def test_sub():          assert (X() - Y()).evaluate(PT)    == -1.0
def test_rsub():         assert (3.0 - X()).evaluate(PT)    == 2.0
def test_mul_func():     assert (X() * Y()).evaluate(PT)    == 2.0
def test_mul_scalar():   assert (X() * 3.0).evaluate(PT)   == 3.0
def test_rmul_scalar():  assert (3.0 * X()).evaluate(PT)   == 3.0
def test_div():          assert (Y() / X()).evaluate(PT)    == 2.0
def test_rdiv():         assert (2.0 / Y()).evaluate(PT)   == 1.0
def test_neg():          assert (-X()).evaluate(PT)         == -1.0
def test_pow():          assert (Y() ** 2).evaluate(PT)    == 4.0

def test_scalar_mul_produces_Scale(): assert isinstance(2 * X(), Scale)
def test_func_mul_produces_Mul():     assert isinstance(X() * Y(), Mul)

def test_nested():
    f = X()**2 + Y()**2 + Z()**2
    assert f.evaluate(PT) == pytest.approx(14.0)


# ── XML tag correctness ────────────────────────────────────────────────────

@pytest.mark.parametrize("node,expected_tag", [
    (X(),             "x"),
    (Y(),             "y"),
    (Z(),             "z"),
    (Constant(1.0),   "constant"),
    (Add(X(), Y()),   "add"),
    (Sub(X(), Y()),   "sub"),
    (Mul(X(), Y()),   "mul"),
    (Div(X(), Y()),   "div"),
    (Scale(X(), 2.0), "scale"),
    (Neg(X()),        "neg"),
    (Pow(X(), 2),     "pow"),
    (Sin(X()),        "sin"),
    (Cos(X()),        "cos"),
    (Sqrt(X()),       "sqrt"),
    (Exp(X()),        "exp"),
    (Log(X()),        "log"),
    (Abs(X()),        "abs"),
])
def test_tag(node, expected_tag):
    assert node.to_xml_element([]).tag == expected_tag

def test_constant_value_attr():
    assert float(Constant(3.14).to_xml_element([]).get("value")) == pytest.approx(3.14)

def test_scale_value_attr():
    assert float(Scale(X(), 2.5).to_xml_element([]).get("value")) == pytest.approx(2.5)

def test_pow_value_attr():
    assert float(Pow(X(), 3).to_xml_element([]).get("value")) == pytest.approx(3.0)

def test_binary_node_has_two_children():
    for binary in [Add, Sub, Mul, Div]:
        elem = binary(X(), Y()).to_xml_element([])
        assert len(elem) == 2

def test_unary_node_has_one_child():
    for unary in [Neg, Sin, Cos, Abs, Exp, Log, Sqrt]:
        elem = unary(X()).to_xml_element([])
        assert len(elem) == 1
    for unary in [Scale, Pow]:
        elem = unary(X(), 1).to_xml_element([])
        assert len(elem) == 1

def test_wrong_tag():
    bad_elem = ET.Element("tanh")  # not in the dispatch table
    with pytest.raises(ValueError, match="Unknown tag 'tanh'"):
        ImplicitFunction.from_xml_element(bad_elem)

# ── Cached serialisation ───────────────────────────────────────────────────

def test_single_use_emits_to_cache():
    elem = Cached(X()).to_xml_element([])
    assert elem.tag == "to_cache"
    assert elem.get("id") == "0"

def test_to_cache_wraps_child():
    elem = Cached(X()).to_xml_element([])
    assert len(elem) == 1
    assert elem[0].tag == "x"

def test_shared_node_first_to_cache_then_from_cache():
    node     = Cached(X())
    _cached  = []
    first    = node.to_xml_element(_cached)
    second   = node.to_xml_element(_cached)
    assert first.tag  == "to_cache"
    assert second.tag == "from_cache"
    assert first.get("id") == second.get("id") == "0"

def test_distinct_cached_nodes_get_distinct_ids():
        node1 = Cached(X())
        node2 = Cached(Y())
        _cached = []
        ex = node1.to_xml_element(_cached)
        ey = node2.to_xml_element(_cached)
        assert ex.get("id") == "0"
        assert ey.get("id") == "1"

def test_cached_in_expression():
    cx   = Cached(X())
    elem = (cx + cx).to_xml_element([])
    assert elem.tag    == "add"
    assert elem[0].tag == "to_cache"
    assert elem[1].tag == "from_cache"
    assert elem[0].get("id") == elem[1].get("id")


# ── XML round-trip ─────────────────────────────────────────────────────────

@pytest.mark.parametrize("func,expected", [
    (X(),                        1.0),
    (Y(),                        2.0),
    (Z(),                        3.0),
    (Constant(7.0),              7.0),
    (Add(X(), Y()),              3.0),
    (Sub(Z(), X()),              2.0),
    (Mul(X(), Y()),              2.0),
    (Div(Y(), X()),              2.0),
    (Scale(Z(), 2.0),            6.0),
    (Neg(X()),                  -1.0),
    (Pow(Y(), 2),                4.0),
    (Sin(Constant(0.0)),         0.0),
    (Cos(Constant(0.0)),         1.0),
    (Sqrt(Constant(9.0)),        3.0),
    (Exp(Constant(0.0)),         1.0),
    (Log(Constant(1.0)),         0.0),
    (Abs(Neg(Y())),              2.0),
    (X()**2 + Y()**2 + Z()**2,  14.0),
])
def test_roundtrip(func, expected):
    assert xml_roundtrip(func) == pytest.approx(expected)

def test_roundtrip_cached_single():
    cx = Cached(X())
    f  = Sin(cx) + Cos(cx)
    assert xml_roundtrip(f) == pytest.approx(np.sin(1.0) + np.cos(1.0))

def test_roundtrip_cached_shared():
    cx = Cached(2.0 * X())
    f  = cx * cx                   # (2*1)^2 = 4
    assert xml_roundtrip(f) == pytest.approx(4.0)

def test_roundtrip_deeply_nested():
    f = Sin(Cos(X() + Constant(np.pi)))
    assert xml_roundtrip(f) == pytest.approx(f.evaluate(PT))


# ── TPMS-like integration smoke test ──────────────────────────────────────

@pytest.fixture
def gyroid():
    L   = 1.0
    k   = 2 * np.pi / L
    cx  = Cached(k * X())
    cy  = Cached(k * Y())
    cz  = Cached(k * Z())
    return Sin(cx)*Cos(cz) + Sin(cy)*Cos(cx) + Sin(cz)*Cos(cy)

def test_gyroid_evaluates(gyroid):
    # value at origin should be 0 for standard gyroid
    assert gyroid.evaluate((0.0, 0.0, 0.0)) == pytest.approx(0.0)

def test_gyroid_roundtrip(gyroid):
    pt = (0.1, 0.2, 0.3)
    assert xml_roundtrip(gyroid, pt) == pytest.approx(gyroid.evaluate(pt))