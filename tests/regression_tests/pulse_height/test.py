import numpy as np
import openmc

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def sphere_model():
    model = openmc.model.Model()
    # Define materials
    NaI = openmc.Material(material_id=1)
    NaI.set_density('g/cc', 3.7)
    NaI.add_element('Na', 1.0)
    NaI.add_element('I', 1.0)

    model.materials = openmc.Materials([NaI])

    # Define geometry: two spheres in each other
    s1 = openmc.Sphere(r=1, surface_id=1)
    s2 = openmc.Sphere(r=2, surface_id=2, boundary_type='vacuum')
    inner_sphere = openmc.Cell(cell_id=1, name='inner sphere')
    inner_sphere.region = -s1
    inner_sphere.fill = NaI
    outer_sphere = openmc.Cell(cell_id=2, name='outer sphere')
    outer_sphere.region = +s1 & -s2
    root = openmc.Universe(universe_id=0, name='root universe')
    root.add_cell(inner_sphere)
    root.add_cell(outer_sphere)
    model.geometry = openmc.Geometry(root)

    # Define settings
    settings_file = openmc.Settings()
    settings_file.run_mode = 'fixed source'
    settings_file.batches = 5
    settings_file.particles = 100
    settings_file.photon_transport = True
    settings_file.source = openmc.source.Source(space=openmc.stats.Point(),
                                                energy=openmc.stats.Discrete([1e6],[1]),
                                                particle='photon')
    model.settings = settings_file

    # Define tallies
    tallies = openmc.Tallies()
    tally = openmc.Tally(name="pht tally")
    tally.scores = ['pulse-height']
    cell_filter = openmc.CellFilter(inner_sphere)
    energy_filter = openmc.EnergyFilter(np.linspace(0, 1_000_000, 101))
    tally.filters = [cell_filter, energy_filter]
    tallies.append(tally)
    model.tallies = tallies

    return model


def test_periodic(box_model):
    harness = PyAPITestHarness('statepoint.5.h5', sphere_model)
    harness.main()