import openmc
import pytest
import numpy as np
 
from openmc.surface import ImplicitSurface
from openmc.implicit import X, Y, Z, Cos, Sin, Cached
from tests.testing_harness import PyAPITestHarness
 
 
@pytest.fixture()
def implicit_fuel_graded_nocache_model():
    # Material
    material = openmc.Material()
    material.add_nuclide('U235',  1.0)
    material.set_density('g/cm3', 16.0)
 
    # Box
    x0 = openmc.XPlane(-16., boundary_type="vacuum")
    x1 = openmc.XPlane(+16., boundary_type="vacuum")
    y0 = openmc.YPlane(-16., boundary_type="vacuum")
    y1 = openmc.YPlane(+16., boundary_type="vacuum")
    z0 = openmc.ZPlane(-16., boundary_type="vacuum")
    z1 = openmc.ZPlane(+16., boundary_type="vacuum")
    box = +x0 & -x1 & +y0 & -y1 & +z0 & -z1
    # Fuel grading
    invPitch = (Z() ** 2 + 8.) / 64.
    x = 2 * np.pi * X() * invPitch
    y = 2 * np.pi * Y() * invPitch
    z = 2 * np.pi * Z() * invPitch
    func = Cos(x) + Cos(y) + Cos(z)
    impl = ImplicitSurface(function=func)

    fuel_cell = openmc.Cell(region=-impl & box, fill=material)
    void_cell = openmc.Cell(region=+impl & box)
    geometry  = openmc.Geometry([fuel_cell, void_cell])
 
    # Settings
    settings = openmc.Settings()
    settings.particles = 1000
    settings.batches   = 20
    settings.inactive  = 5
 
    model = openmc.Model(settings=settings, geometry=geometry)
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point((0., 0., 0.)),
    )
    return model
 
 
def test_implicit_fuel_graded_nocache(implicit_fuel_graded_nocache_model):
    harness = PyAPITestHarness('statepoint.20.h5', model=implicit_fuel_graded_nocache_model)
    harness.main()