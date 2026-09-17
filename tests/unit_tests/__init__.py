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


def assert_sample_mean(samples, expected_mean):
    # Calculate sample standard deviation
    std_dev = samples.std() / np.sqrt(samples.size - 1)

    # Means should agree within 4 sigma 99.993% of the time. Note that this is
    # expected to fail about 1 out of 16,000 times
    assert np.abs(expected_mean - samples.mean()) < 4*std_dev
