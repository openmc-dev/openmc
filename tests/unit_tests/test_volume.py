import numpy as np
import pytest

import openmc


def test_infinity_handling():
    surf1 = openmc.Sphere(boundary_type="vacuum")
    reg1 = -surf1
    cell1 = openmc.Cell(region=reg1)

    lower_left = [-2, -np.inf, -2]
    upper_right = [np.inf, 2, 2]

    msg = "Infinite value found in lower_left or upper_right. Could not compute volume."

    with pytest.raises(ValueError, match=msg):
        openmc.VolumeCalculation([cell1], 100, lower_left, upper_right)
