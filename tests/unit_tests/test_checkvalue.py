from numbers import Real

import numpy as np
import pytest

from openmc.checkvalue import check_iterable_type


@pytest.mark.parametrize('dtype', (np.float16, np.float32, np.float64))
@pytest.mark.parametrize('expected_type', (Real, float))
def test_check_iterable_type_float_array(dtype, expected_type):
    values = np.ones((2, 3, 4), dtype=dtype)
    check_iterable_type(
        'values', values, expected_type, min_depth=1, max_depth=3)


def test_check_iterable_type_float_array_depth():
    values = np.ones((2, 3))
    with pytest.raises(TypeError, match='maximum depth'):
        check_iterable_type('values', values, Real, max_depth=1)


def test_check_iterable_type_nonfloat_array():
    values = np.ones(3, dtype=np.complex128)
    with pytest.raises(TypeError, match='Items must be of type'):
        check_iterable_type('values', values, Real)
