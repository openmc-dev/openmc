import openmc
import pytest


@pytest.fixture
def th232_model():
    # URR boundaries for Th232
    e_min, e_max = 4000.0, 100000.0

    model = openmc.model.Model()
    th232 = openmc.Material()
    th232.add_nuclide('Th232', 1.0)

    surf = openmc.Sphere(r=100.0, boundary_type='reflective')
    cell = openmc.Cell(fill=th232, region=-surf)
    model.geometry = openmc.Geometry([cell])

    model.settings.particles = 100
    model.settings.batches = 10
    model.settings.run_mode = 'fixed source'
    energies = openmc.stats.Uniform(e_min, e_max)
    model.settings.source = openmc.IndependentSource(energy=energies)

    tally = openmc.Tally(name='rates')
    tally.filters = [openmc.EnergyFilter([e_min, e_max])]
    tally.scores = ['(n,gamma)', 'absorption', 'fission']
    model.tallies.append(tally)
    return model


def test_urr_capture(run_in_tmpdir, th232_model):
    # Export and run model
    th232_model.export_to_xml()
    openmc.run()

    # Get reaction rates from tally
    with openmc.StatePoint('statepoint.10.h5') as sp:
        t = sp.get_tally(name='rates')
        ngamma, absorption, fission = t.mean.flatten()

    # In URR, the (n,gamma) rate should be equal to absorption - fission
    assert ngamma == pytest.approx(absorption - fission)
