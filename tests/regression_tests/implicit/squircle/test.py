import openmc
import pytest
import numpy as np
 
from openmc.surface import ImplicitSurface
from openmc.implicit import X, Y, Z, Cos, Sqrt, Cached
from tests.testing_harness import PyAPITestHarness
 
 
@pytest.fixture()
def implicit_sphere_model():
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
    # Squircle
    L = 8.0
    scale = 0.5/L
    xp, yp, zp = Cached(scale * X()),  Cached(scale * Y()),  Cached(scale * Z())
    x = 2 * np.pi * xp * Sqrt(1 - 0.5 * yp**2 - 0.5 * zp**2 + yp**2*zp**2/3)
    y = 2 * np.pi * yp * Sqrt(1 - 0.5 * xp**2 - 0.5 * zp**2 + xp**2*zp**2/3)
    z = 2 * np.pi * zp * Sqrt(1 - 0.5 * xp**2 - 0.5 * yp**2 + xp**2*yp**2/3)
    func = Cos(x) + Cos(y) + Cos(z)
    impl = ImplicitSurface(function=func)
    sphere = openmc.Sphere(r=2*L)

    fuel_cell = openmc.Cell(region=-impl & -sphere, fill=material)
    void_cell = openmc.Cell(region=+impl & -sphere)
    outer_cell = openmc.Cell(region=box & +sphere)
    geometry  = openmc.Geometry([fuel_cell, void_cell, outer_cell])
 
    # Settings
    settings = openmc.Settings()
    settings.particles = 1000
    settings.batches   = 20
    settings.inactive  = 5
    settings.implicit = {"name": "fast", "atol":1e-10, "ftol":1e-10, "margin":5e-9}
 
    model = openmc.Model(settings=settings, geometry=geometry)
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point((0., 0., 0.)),
    )
    return model
 
 
def test_implicit_sphere(implicit_sphere_model):
    harness = PyAPITestHarness('statepoint.20.h5', model=implicit_sphere_model)
    harness.main()