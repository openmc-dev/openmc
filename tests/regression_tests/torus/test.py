import openmc
import pytest

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def model():
    model = openmc.Model()
    fuel = openmc.Material()
    fuel.set_density("g/cm3", 12.0)
    fuel.add_nuclide("U235", 1.0)
    al = openmc.Material()
    al.set_density("g/cm3", 1.0)
    al.add_nuclide("H1", 1.0)
    model.materials.extend([fuel, al])

    # 🍩🍩🍩
    zt = openmc.ZTorus(a=3, b=1.5, c=1)
    xt = openmc.XTorus(x0=6, a=3, b=1.5, c=1)
    yt = openmc.YTorus(x0=6, a=6, b=1, c=0.75)
    box = openmc.model.RectangularParallelepiped(
        -5, 14, -5, 5, -8, 8, boundary_type="vacuum"
    )

    xt_cell = openmc.Cell(fill=fuel, region=-xt)
    yt_cell = openmc.Cell(fill=fuel, region=-yt)
    zt_cell = openmc.Cell(fill=fuel, region=-zt)
    outer_cell = openmc.Cell(fill=al, region=-box & +xt & +yt & +zt)
    model.geometry = openmc.Geometry([xt_cell, yt_cell, zt_cell, outer_cell])

    model.settings.particles = 1000
    model.settings.batches = 10
    model.settings.inactive = 5
    return model


def test_torus(model):
    harness = PyAPITestHarness("statepoint.10.h5", model)
    harness.main()
