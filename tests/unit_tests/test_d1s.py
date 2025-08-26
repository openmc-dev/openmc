from pathlib import Path
from math import exp

import numpy as np
import pytest
import openmc
import openmc.deplete
from openmc.deplete import d1s


CHAIN_PATH = Path(__file__).parents[1] / "chain_ni.xml"


@pytest.fixture
def model():
    """Simple model with natural Ni"""
    mat = openmc.Material()
    mat.add_element('Ni', 1.0)
    geom = openmc.Geometry([openmc.Cell(fill=mat)])
    return openmc.Model(geometry=geom)


def test_get_radionuclides(model):
    # Check that radionuclides are correct and are unstable
    chain = openmc.deplete.Chain.from_xml(CHAIN_PATH)
    nuclides = d1s.get_radionuclides(model, chain)
    assert sorted(nuclides) == [
        'Co58', 'Co60', 'Co61', 'Co62', 'Co64',
        'Fe55', 'Fe59', 'Fe61', 'Ni57', 'Ni59', 'Ni63', 'Ni65'
    ]
    for nuc in nuclides:
        assert openmc.data.half_life(nuc) is not None


@pytest.mark.parametrize("nuclide", ['Co60', 'Ni63', 'H3', 'Na24', 'K40'])
def test_time_correction_factors(nuclide):
    # Irradiation schedule turning unit neutron source on and off
    timesteps = [1.0, 1.0, 1.0]
    source_rates = [1.0, 0.0, 1.0]

    # Compute expected solution
    decay_rate = openmc.data.decay_constant(nuclide)
    g = exp(-decay_rate)
    expected = [0.0, (1 - g), (1 - g)*g, (1 - g)*(1 + g*g)]

    # Test against expected solution
    tcf = d1s.time_correction_factors([nuclide], timesteps, source_rates)
    assert tcf[nuclide] == pytest.approx(expected)

    # Make sure all values at first timestep and onward are positive (K40 case
    # has very small decay constant that stresses this)
    assert np.all(tcf[nuclide][1:] > 0.0)

    # Timesteps as a tuple
    timesteps = [(1.0, 's'), (1.0, 's'), (1.0, 's')]
    tcf = d1s.time_correction_factors([nuclide], timesteps, source_rates)
    assert tcf[nuclide] == pytest.approx(expected)

    # Test changing units
    timesteps = [1.0/60.0, 1.0/60.0, 1.0/60.0]
    tcf = d1s.time_correction_factors([nuclide], timesteps, source_rates,
                                      timestep_units='min')
    assert tcf[nuclide] == pytest.approx(expected)


def test_prepare_tallies(model):
    tally = openmc.Tally()
    tally.filters = [openmc.ParticleFilter('photon')]
    tally.scores = ['flux']
    model.tallies = [tally]

    # Check that prepare_tallies adds a ParentNuclideFilter
    nuclides = ['Co58', 'Co60', 'Fe55']
    d1s.prepare_tallies(model, nuclides, chain_file=CHAIN_PATH)
    assert tally.contains_filter(openmc.ParentNuclideFilter)
    assert list(tally.filters[-1].bins) == nuclides

    # Get rid of parent nuclide filter
    tally.filters.pop()

    # With no nuclides specified, filter should use get_radionuclides
    radionuclides = d1s.get_radionuclides(model, CHAIN_PATH)
    d1s.prepare_tallies(model, chain_file=CHAIN_PATH)
    assert tally.contains_filter(openmc.ParentNuclideFilter)
    assert sorted(tally.filters[-1].bins) == sorted(radionuclides)

    assert len(tally.filters) == 2
    # calling prepare_tallies twice should not add another ParentNuclideFilter
    d1s.prepare_tallies(model, chain_file=CHAIN_PATH)
    assert len(tally.filters) == 2


def test_apply_time_correction(run_in_tmpdir):
    # Make simple sphere model with elemental Ni
    mat = openmc.Material()
    mat.add_element('Ni', 1.0)
    sphere = openmc.Sphere(r=10.0, boundary_type='vacuum')
    cell = openmc.Cell(fill=mat, region=-sphere)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])
    model.settings.run_mode = 'fixed source'
    model.settings.batches = 3
    model.settings.particles = 10
    model.settings.photon_transport = True
    model.settings.use_decay_photons = True
    particle_filter = openmc.ParticleFilter('photon')
    tally = openmc.Tally()
    tally.filters = [particle_filter]
    tally.scores = ['flux']
    model.tallies = [tally]

    # Prepare tallies for D1S and compute time correction factors
    nuclides = d1s.prepare_tallies(model, chain_file=CHAIN_PATH)
    factors = d1s.time_correction_factors(nuclides, [1.0e10], [1.0])

    # Run OpenMC and get tally result
    with openmc.config.patch('chain_file', CHAIN_PATH):
        output_path = model.run()
    with openmc.StatePoint(output_path) as sp:
        tally = sp.tallies[tally.id]
        flux = tally.mean.flatten()

    # Apply TCF and make sure results are consistent
    result = d1s.apply_time_correction(tally, factors, sum_nuclides=False)
    tcf = np.array([factors[nuc][-1] for nuc in nuclides])
    assert result.mean.flatten() == pytest.approx(tcf * flux)

    # Make sure summed results match a manual sum
    result_summed = d1s.apply_time_correction(tally, factors)
    assert result_summed.mean.flatten()[0] == pytest.approx(result.mean.sum())

    # Make sure various tally methods work
    result.get_values()
    result_summed.get_values()
    result.get_reshaped_data()
    result_summed.get_reshaped_data()
    result.get_pandas_dataframe()
    result_summed.get_pandas_dataframe()
