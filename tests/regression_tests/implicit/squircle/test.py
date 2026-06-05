import openmc
import pytest
import numpy as np
 
from openmc.surface import ImplicitSurface
from openmc.implicit import X, Y, Z, Cos, Sqrt
from tests.testing_harness import PyAPITestHarness
 
 
@pytest.fixture()
def implicit_sphere_model():
    # Material
    material = openmc.Material()
    material.add_nuclide('U235',  1.0)
    material.set_density('g/cm3', 16.0)
 
    # Box
    x0 = openmc.XPlane(-1., boundary_type="reflective")
    x1 = openmc.XPlane(+1., boundary_type="reflective")
    y0 = openmc.YPlane(-1., boundary_type="reflective")
    y1 = openmc.YPlane(+1., boundary_type="reflective")
    z0 = openmc.ZPlane(-1., boundary_type="reflective")
    z1 = openmc.ZPlane(+1., boundary_type="reflective")
    box = +x0 & -x1 & +y0 & -y1 & +z0 & -z1
    # Squircle
    x = X() * Sqrt(1 - Y()**2/2 - Z()**2/2 + Y()**2*Z()**2/3)
    y = Y() * Sqrt(1 - X()**2/2 - Z()**2/2 + X()**2*Z()**2/3)
    z = Z() * Sqrt(1 - X()**2/2 - Y()**2/2 + X()**2*Y()**2/3)
    L = 0.5
    xp, yp, zp = 2*np.pi*x/L, 2*np.pi*y/L, 2*np.pi*z/L
    func = Cos(xp) + Cos(yp) + Cos(zp)
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