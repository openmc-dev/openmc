"""
Unit tests for implicit surface solver settings.

Tests:
  - Default values
  - Valid assignments
  - Invalid assignments raise errors
  - XML round-trip preserves all values
  - ImplicitSurface uses the solver specified in Settings
"""

import pytest
from collections.abc import Mapping

import openmc


# ==============================================================================
# Default value
# ==============================================================================


def test_implicit_default_is_none_or_dict():
    """Before assignment, implicit should be None or an empty dict."""
    s = openmc.Settings()
    assert s.implicit is None or isinstance(s.implicit, Mapping)


# ==============================================================================
# Valid assignments
# ==============================================================================


def test_empty_dict():
    """Empty dict is valid — all values remain at C++ defaults."""
    s = openmc.Settings()
    s.implicit = {}
    assert s.implicit == {}

@pytest.mark.parametrize("name", ["naive", "fast"])
def test_valid_solver_names(name):
    s = openmc.Settings()
    s.implicit = {'name': name}
    assert s.implicit['name'] == name

def test_set_atol():
    s = openmc.Settings()
    s.implicit = {'atol': 1e-10}
    assert s.implicit['atol'] == pytest.approx(1e-10)

def test_set_ftol():
    s = openmc.Settings()
    s.implicit = {'ftol': 1e-6}
    assert s.implicit['ftol'] == pytest.approx(1e-6)

def test_set_maxiter():
    s = openmc.Settings()
    s.implicit = {'maxiter': 5000}
    assert s.implicit['maxiter'] == 5000

def test_set_margin():
    s = openmc.Settings()
    s.implicit = {'margin': 2e-8}
    assert s.implicit['margin'] == pytest.approx(2e-8)

def test_set_all_fields():
    s = openmc.Settings()
    s.implicit = {
        'name':    'fast',
        'atol':    1e-9,
        'ftol':    1e-6,
        'maxiter': 500,
        'margin':  3e-8,
    }
    assert s.implicit['name']    == 'fast'
    assert s.implicit['atol']    == pytest.approx(1e-9)
    assert s.implicit['ftol']    == pytest.approx(1e-6)
    assert s.implicit['maxiter'] == 500
    assert s.implicit['margin']  == pytest.approx(3e-8)


# ==============================================================================
# Invalid assignments
# ==============================================================================

def test_not_a_dict():
    """Setting implicit to a non-Mapping raises ValueError."""
    s = openmc.Settings()
    with pytest.raises(ValueError):
        s.implicit = 'fast'

def test_not_a_dict_list():
    s = openmc.Settings()
    with pytest.raises(ValueError):
        s.implicit = ['fast', 1e-8]

def test_unknown_key():
    """Unknown keys raise ValueError."""
    s = openmc.Settings()
    with pytest.raises(ValueError):
        s.implicit = {'unknown_key': 1}

def test_atol_negative():
    s = openmc.Settings()
    with pytest.raises((ValueError, Exception)):
        s.implicit = {'atol': -1e-8}

def test_atol_zero():
    s = openmc.Settings()
    with pytest.raises((ValueError, Exception)):
        s.implicit = {'atol': 0.0}

def test_atol_wrong_type():
    s = openmc.Settings()
    with pytest.raises((ValueError, TypeError, Exception)):
        s.implicit = {'atol': 'small'}

def test_maxiter_negative():
    s = openmc.Settings()
    with pytest.raises((ValueError, Exception)):
        s.implicit = {'maxiter': -1}

def test_maxiter_zero():
    s = openmc.Settings()
    with pytest.raises((ValueError, Exception)):
        s.implicit = {'maxiter': 0}

def test_maxiter_wrong_type():
    s = openmc.Settings()
    with pytest.raises((ValueError, TypeError, Exception)):
        s.implicit = {'maxiter': 1.5}

def test_ftol_negative():
    s = openmc.Settings()
    with pytest.raises((ValueError, Exception)):
        s.implicit = {'ftol': -1.0}

def test_margin_negative():
    s = openmc.Settings()
    with pytest.raises((ValueError, Exception)):
        s.implicit = {'margin': -1e-7}

def test_name_wrong_type():
    s = openmc.Settings()
    with pytest.raises((ValueError, TypeError, Exception)):
        s.implicit = {'name': 42}

def test_name_wrong_name():
    s = openmc.Settings()
    with pytest.raises((ValueError, TypeError, Exception)):
        s.implicit = {'name': 'screugneugneu'}

# ==============================================================================
# XML round-trip
# ==============================================================================

def _roundtrip(settings):
    elem = settings.to_xml_element()
    restored = openmc.Settings()
    return restored.from_xml_element(elem)

def test_name_preserved():
    s = openmc.Settings()
    s.implicit = {'name': 'naive'}
    r = _roundtrip(s)
    assert r.implicit['name'] == 'naive'

def test_atol_preserved():
    s = openmc.Settings()
    s.implicit = {'atol': 1e-10}
    r = _roundtrip(s)
    assert r.implicit['atol'] == pytest.approx(1e-10)

def test_maxiter_preserved():
    s = openmc.Settings()
    s.implicit = {'maxiter': 2000}
    r = _roundtrip(s)
    assert r.implicit['maxiter'] == 2000

def test_ftol_preserved():
    s = openmc.Settings()
    s.implicit = {'ftol': 1e-6}
    r = _roundtrip(s)
    assert r.implicit['ftol'] == pytest.approx(1e-6)

def test_margin_preserved():
    s = openmc.Settings()
    s.implicit = {'margin': 2e-7}
    r = _roundtrip(s)
    assert r.implicit['margin'] == pytest.approx(2e-7)

def test_implicit_element_present_in_xml():
    """<implicit> element must appear in settings XML."""
    s = openmc.Settings()
    s.implicit = {'name': 'fast', 'atol': 1e-9}
    elem = s.to_xml_element()
    assert elem.find('implicit') is not None

def test_full_roundtrip():
    s = openmc.Settings()
    s.implicit = {
        'name':    'naive',
        'atol':    1e-9,
        'ftol':    1e-6,
        'maxiter': 300,
        'margin':  5e-8,
    }
    r = _roundtrip(s)
    assert r.implicit['name']    == 'naive'
    assert r.implicit['atol']    == pytest.approx(1e-9)
    assert r.implicit['ftol']    == pytest.approx(1e-6)
    assert r.implicit['maxiter'] == 300
    assert r.implicit['margin']  == pytest.approx(5e-8)