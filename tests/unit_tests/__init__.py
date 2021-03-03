import shutil

import numpy as np
import pytest


# Check if NJOY is available
needs_njoy = pytest.mark.skipif(shutil.which('njoy') is None,
                                reason="NJOY not installed")


def assert_unbounded(obj):
    """Assert that a region/cell is unbounded."""
    ll, ur = obj.bounding_box
    assert ll == pytest.approx((-np.inf, -np.inf, -np.inf))
    assert ur == pytest.approx((np.inf, np.inf, np.inf))
