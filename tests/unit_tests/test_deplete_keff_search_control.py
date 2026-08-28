"""Unit tests for openmc.deplete.keff_search_control."""

import numpy as np
import pytest

from openmc.deplete.keff_search_control import _KeffSearchControl


class MockOperator:
    """Minimal operator recording calls to _update_materials_and_nuclides."""

    def __init__(self, calls):
        self.calls = calls

    def _update_materials_and_nuclides(self, vec):
        self.calls.append(('update_materials', [v.copy() for v in vec]))


@pytest.fixture
def control_and_calls(monkeypatch):
    calls = []
    operator = MockOperator(calls)
    control = _KeffSearchControl(
        operator, lambda x: None, x0=0.0, x1=1.0, bracket=[0.0, 2.0])

    def fake_search():
        calls.append(('search', None))
        return 0.5

    def fake_update_vec(x):
        calls.append(('update_vec', None))

    monkeypatch.setattr(control, '_search_for_keff', fake_search)
    monkeypatch.setattr(control, '_update_vec', fake_update_vec)
    return control, calls


def test_materials_updated_before_search(control_and_calls):
    """Compositions must be pushed to openmc.lib before the search runs.

    The keff search is executed at the beginning of a depletion step, before
    the transport operator is called. Without an explicit update, both
    openmc.lib.materials and the operator's AtomNumber still hold the previous
    call's compositions, and _update_vec() overwrites the depleted vector with
    them -- freezing nuclide densities at their initial values.
    """
    control, calls = control_and_calls

    n = [np.array([1.0, 2.0, 3.0]), np.array([4.0, 5.0, 6.0])]
    root = control.run(n)

    assert root == 0.5
    assert [name for name, _ in calls] == [
        'update_materials', 'search', 'update_vec']

    # The vector handed to the operator must be the current composition
    recorded = calls[0][1]
    assert len(recorded) == len(n)
    for actual, expected in zip(recorded, n):
        np.testing.assert_array_equal(actual, expected)


def test_bracket_validation():
    operator = MockOperator([])
    with pytest.raises(ValueError, match='exactly 2 elements'):
        _KeffSearchControl(operator, lambda x: None, 0.0, 1.0, [0.0])
    with pytest.raises(ValueError, match=r'bracket\[0\] must be'):
        _KeffSearchControl(operator, lambda x: None, 0.0, 1.0, [2.0, 1.0])
