import openmc
import pytest
 
from openmc.surface import ImplicitSurface
from openmc.implicit import X, Y, Z
from tests.testing_harness import PyAPITestHarness
 
 
@pytest.fixture()
def implicit_sphere_model():
    # Material
    material = openmc.Material()
    material.add_nuclide('U235', 1.0)
    material.set_density('g/cm3', 16.0)
 
    # Geometry: implicit sphere enclosed in an analytic vacuum sphere.
    # The outer analytic sphere provides the finite region required by
    # the two-pass Region::distance algorithm.
    R     = 10.0
    outer = openmc.Sphere(r=R * 1.1, boundary_type='vacuum')
    impl  = ImplicitSurface(function=X()**2 + Y()**2 + Z()**2, isovalue=R**2)
 
    fuel_cell = openmc.Cell(region=-impl & -outer, fill=material)
    void_cell = openmc.Cell(region=+impl & -outer)
    geometry  = openmc.Geometry([fuel_cell, void_cell])
 
    # Settings
    settings = openmc.Settings()
    settings.particles = 1000
    settings.batches   = 20
    settings.inactive  = 5
 
    model = openmc.Model(settings=settings, geometry=geometry)
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point((0., 0., 0.)),
        constraints={'fissionable': True},
    )
    return model
 
 
def test_implicit_sphere(implicit_sphere_model):
    harness = PyAPITestHarness('statepoint.20.h5', model=implicit_sphere_model)
    harness.main()