"""Standard deviation of a tally whose realizations are nearly identical."""

import warnings

import numpy as np
import pytest

import openmc


def _tally_with_sums(sum_, sum_sq, n):
    """Build a Tally carrying given accumulated sums, bypassing statepoint I/O."""
    tally = openmc.Tally()
    tally._sp_filename = 'unused.h5'
    tally._results_read = True
    tally._num_realizations = n
    tally._sum = np.array(sum_, dtype=float).reshape(1, 1, 1)
    tally._sum_sq = np.array(sum_sq, dtype=float).reshape(1, 1, 1)
    tally._shape = (1, 1, 1)
    return tally


def test_std_dev_negative_variance_from_roundoff():
    """Rounding can drive sum_sq/n below mean**2; std_dev must not be NaN.

    A tally whose batches all score the same value has zero variance, but it is
    computed as a difference of accumulated sums and can land just below zero.
    """
    # 15 realizations of the same value; the accumulated sums land such that
    # sum_sq/n is one ulp below mean**2. Whether a given pair of sums cancels
    # this way depends on the exact order of operations, so the inputs are
    # pinned rather than recomputed, and the premise is asserted below.
    n = 15
    sum_ = 142569.55519193487
    sum_sq = 1355071871.175077
    mean = np.array(sum_).reshape(1, 1, 1) / n
    assert (np.array(sum_sq).reshape(1, 1, 1) / n - mean ** 2).ravel()[0] < 0.0

    tally = _tally_with_sums(sum_, sum_sq, n)
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        std_dev = tally.std_dev
    assert not np.isnan(std_dev).any()
    assert std_dev.ravel()[0] == 0.0


def test_std_dev_single_realization():
    """One realization reports an undefined uncertainty, and says why.

    The standard deviation stays NaN, as it was when 1/(n - 1) divided by zero,
    but the user gets a statement of the problem instead of a numpy warning
    about an invalid value.
    """
    tally = _tally_with_sums(3.0, 9.0, 1)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        std_dev = tally.std_dev

    assert np.isnan(std_dev.ravel()[0])
    messages = [str(w.message) for w in caught]
    assert any('single realization' in m for m in messages), messages
    assert not any(issubclass(w.category, RuntimeWarning) for w in caught)


def test_std_dev_single_realization_zero_mean():
    """Bins that never scored stay at zero, as before."""
    tally = _tally_with_sums(0.0, 0.0, 1)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        std_dev = tally.std_dev
    assert std_dev.ravel()[0] == 0.0


def test_std_dev_unchanged_for_ordinary_tally():
    """Normal statistics are untouched by the clamp."""
    n = 5
    values = np.array([1.0, 1.2, 0.9, 1.1, 1.3])
    sum_, sum_sq = values.sum(), (values ** 2).sum()
    tally = _tally_with_sums(sum_, sum_sq, n)
    expected = np.sqrt((sum_sq / n - (sum_ / n) ** 2) / (n - 1))
    assert tally.std_dev.ravel()[0] == pytest.approx(expected, rel=1e-12)
