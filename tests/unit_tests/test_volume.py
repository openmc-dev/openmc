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


@pytest.mark.parametrize('cls', [openmc.Cell, openmc.Material, openmc.Universe])
def test_invalid_id(run_in_tmpdir, cls):
    m = openmc.Material()
    m.add_nuclide('U235', 0.02)
    sph = openmc.Sphere(boundary_type='vacuum')
    cell = openmc.Cell(fill=m, region=-sph)
    model = openmc.Model(geometry=openmc.Geometry([cell]))

    # Apply volume calculation with unused domains
    model.settings.volume_calculations = openmc.VolumeCalculation(
        [cls()], 10000, *model.geometry.bounding_box)

    with pytest.raises(RuntimeError):
        model.calculate_volumes()
