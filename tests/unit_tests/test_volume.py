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


@pytest.mark.parametrize("cls", [openmc.Cell, openmc.Material, openmc.Universe])
def test_invalid_id(run_in_tmpdir, cls):
    m = openmc.Material()
    m.add_nuclide("U235", 0.02)
    sph = openmc.Sphere(boundary_type="vacuum")
    cell = openmc.Cell(fill=m, region=-sph)
    model = openmc.Model(geometry=openmc.Geometry([cell]))

    # Apply volume calculation with unused domains
    model.settings.volume_calculations = openmc.VolumeCalculation(
        [cls()], 10000, *model.geometry.bounding_box
    )

    with pytest.raises(RuntimeError):
        model.calculate_volumes()


def test_no_bcs(run_in_tmpdir):
    """Ensure that a model without boundary conditions can be used in a volume calculation"""
    model = openmc.examples.pwr_pin_cell()
    for surface in model.geometry.get_all_surfaces().values():
        surface.boundary_type = "transmission"

    bbox = openmc.BoundingBox([-1.0] * 3, [1.0] * 3)
    cells = list(model.geometry.get_all_cells().values())
    vc = openmc.VolumeCalculation(
        cells, samples=10, lower_left=bbox[0], upper_right=bbox[1]
    )

    model.settings.volume_calculations = [vc]
    model.calculate_volumes()
