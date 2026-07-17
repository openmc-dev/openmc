import numpy as np
import pytest
import openmc


@pytest.fixture
def model():
    openmc.reset_auto_ids()
    model = openmc.Model()

    # Define materials
    NaI = openmc.Material()
    NaI.set_density('g/cm3', 3.7)
    NaI.add_element('Na', 1.0)
    NaI.add_element('I', 1.0)

    # Define geometry: two spheres in each other
    s1 = openmc.Sphere(r=1)
    s2 = openmc.Sphere(r=2, boundary_type='vacuum')
    inner_sphere = openmc.Cell(name='inner sphere', fill=NaI, region=-s1)
    outer_sphere = openmc.Cell(name='outer sphere', fill=NaI, region=+s1 & -s2)
    model.geometry = openmc.Geometry([inner_sphere, outer_sphere])

    # Define settings
    model.settings.run_mode = 'fixed source'
    model.settings.batches = 1
    model.settings.particles = 10000
    model.settings.photon_transport = True
    model.settings.source = openmc.IndependentSource(
        energy=openmc.stats.delta_function(1e6),
        particle='photon'
    )

    # Define tallies
    energy_filter = openmc.EnergyFilter([1e3, 1e7])    

    tally1 = openmc.Tally()
    tally1.scores = ['pulse-height']
    cell_filter1 = openmc.CellFilter([inner_sphere, outer_sphere])
    tally1.filters = [cell_filter1, energy_filter]

    tally2 = openmc.Tally()
    tally2.scores = ['pulse-height']
    cell_filter2 = openmc.CellFilter([outer_sphere, inner_sphere])
    tally2.filters = [cell_filter2, energy_filter]

    model.tallies = [tally1, tally2]
    return model


def test_pulse_height(model, run_in_tmpdir):
    sp_path = model.run()
    sp = openmc.StatePoint(sp_path)
    t1 = sp.tallies[1].mean.squeeze()
    t2 = sp.tallies[2].mean.squeeze()
    
    np.testing.assert_array_equal(t1, t2[::-1])


