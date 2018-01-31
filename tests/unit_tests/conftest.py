import openmc
import pytest


@pytest.fixture
def run_in_tmpdir(tmpdir):
    orig = tmpdir.chdir()
    try:
        yield
    finally:
        orig.chdir()


@pytest.fixture(scope='module')
def uo2():
    m = openmc.Material(material_id=100, name='UO2')
    m.add_nuclide('U235', 1.0)
    m.add_nuclide('O16', 2.0)
    m.set_density('g/cm3', 10.0)
    m.depletable = True
    return m


@pytest.fixture(scope='module')
def sphere_model():
    model = openmc.model.Model()

    m = openmc.Material()
    m.add_nuclide('U235', 1.0)
    m.set_density('g/cm3', 1.0)
    model.materials.append(m)

    sph = openmc.Sphere(boundary_type='vacuum')
    c = openmc.Cell(fill=m, region=-sph)
    model.geometry.root_universe = openmc.Universe(cells=[c])

    model.settings.particles = 100
    model.settings.batches = 10
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.Source(space=openmc.stats.Point())
    return model
