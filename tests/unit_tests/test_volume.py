import numpy as np
import pytest

import openmc


def test_infinity_handling():
    surf1 = openmc.Sphere(boundary_type="vacuum")
    cell1 = openmc.Cell(region=-surf1)

    lower_left = (-2, -np.inf, -2)
    upper_right = (np.inf, 2, 2)

    with pytest.raises(ValueError, match="must be finite"):
        openmc.VolumeCalculation([cell1], 100, lower_left, upper_right)
