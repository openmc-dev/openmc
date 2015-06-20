import numpy as np


def is_integer(val):
    return isinstance(val, (int, np.int32, np.int64))


def is_float(val):
    return isinstance(val, (float, np.float32, np.float64))


def is_string(val):
    return isinstance(val, (str, np.str))
