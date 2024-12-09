from pathlib import Path
from math import exp, log

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
    nuclides = d1s.get_radionuclides(model, CHAIN_PATH)
    assert sorted(nuclides) == [
        'Co58', 'Co60', 'Co61', 'Co62', 'Co64',
        'Fe55', 'Fe59', 'Fe61', 'Ni57', 'Ni59', 'Ni63', 'Ni65'
    ]
    for nuc in nuclides:
        assert openmc.data.half_life(nuc) is not None


def test_time_correction_factors():
    nuclides = ['Co60', 'Ni63', 'H3', 'Na24', 'K40']

    # Irradiation schedule turning unit neutron source on and off
    timesteps = [1.0, 1.0, 1.0]
    source_rates = [1.0, 0.0, 1.0]

    # Compute expected solutin
    decay_rate = {x: log(2.0) / openmc.data.half_life(x) for x in nuclides}
    g = {nuc: exp(-r) for nuc, r in decay_rate.items()}
    expected = {}
    for nuc in nuclides:
        gi = g[nuc]
        expected[nuc] = [0.0, (1 - gi), (1 - gi)*gi, (1 - gi)*(1 + gi**2)]

    tcf = d1s.time_correction_factors(nuclides, timesteps, source_rates)
    for nuc, actual in tcf.items():
        assert actual == pytest.approx(expected[nuc])

    # Timesteps as a tuple
    timesteps = [(1.0, 's'), (1.0, 's'), (1.0, 's')]
    tcf = d1s.time_correction_factors(nuclides, timesteps, source_rates)
    for nuc, actual in tcf.items():
        assert actual == pytest.approx(expected[nuc])

    # Test changing units
    timesteps = [1.0/60.0, 1.0/60.0, 1.0/60.0]
    tcf = d1s.time_correction_factors(nuclides, timesteps, source_rates,
                                      timestep_units='min')
    for nuc, actual in tcf.items():
        assert actual == pytest.approx(expected[nuc]), nuc


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


def test_apply_time_correction():
    pass
