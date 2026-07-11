import openmc
import pytest
 
from openmc.surface import ImplicitSurface
from openmc.implicit import X, Y, Z
from tests.testing_harness import PyAPITestHarness
 
 
@pytest.fixture()
def implicit_primitive_model():
    # Material
    material = openmc.Material()
    material.add_nuclide('U235',  5.0)
    material.add_nuclide('U238', 95.0)
    material.set_density('g/cm3', 16.0)
 
    # Primitive
    x0 = openmc.XPlane(-16.0, boundary_type="vacuum")
    x1 = openmc.XPlane(+16.0, boundary_type="vacuum")
    y0 = openmc.YPlane(-16.0, boundary_type="vacuum")
    y1 = openmc.YPlane(+16.0, boundary_type="vacuum")
    z0 = openmc.ZPlane(-16.0, boundary_type="vacuum")
    z1 = openmc.ZPlane(+16.0, boundary_type="vacuum")
    box = +x0 & -x1 & +y0 & -y1 & +z0 & -z1
    
    impl  = openmc.TPMS.from_pitch_isovalue("primitive", 4., 0.)
 
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
 
 
def test_implicit_primitive(implicit_primitive_model):
    harness = PyAPITestHarness('statepoint.20.h5', model=implicit_primitive_model)
    harness.main()