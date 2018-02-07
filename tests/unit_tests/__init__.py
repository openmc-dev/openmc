import numpy as np
import pytest


def assert_unbounded(obj):
    """Assert that a region/cell is unbounded."""
    ll, ur = obj.bounding_box
    assert ll == pytest.approx((-np.inf, -np.inf, -np.inf))
    assert ur == pytest.approx((np.inf, np.inf, np.inf))
