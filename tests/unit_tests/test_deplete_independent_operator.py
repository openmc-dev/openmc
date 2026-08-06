"""Basic unit tests for openmc.deplete.IndependentOperator instantiation

"""

from pathlib import Path

import numpy as np
import pytest

from openmc import Material
from openmc.deplete import IndependentOperator, MicroXS, Chain
from openmc.deplete.independent_operator import _neutrons_emitted

CHAIN_PATH = Path(__file__).parents[1] / "chain_simple.xml"
ONE_GROUP_XS = Path(__file__).parents[1] / "micro_xs_simple.csv"


def test_operator_init():
    """The test uses a temporary dummy chain. This file will be removed
    at the end of the test, and only contains a depletion_chain node."""
    volume = 1
    nuclides = {'U234': 8.922411359424315e+18,
                'U235': 9.98240191860822e+20,
                'U238': 2.2192386373095893e+22,
                'U236': 4.5724195495061115e+18,
                'O16': 4.639065406771322e+22,
                'O17': 1.7588724018066158e+19}
    flux = 1.0
    micro_xs = MicroXS.from_csv(ONE_GROUP_XS)
    chain = Chain.from_xml(CHAIN_PATH)
    IndependentOperator.from_nuclides(
        volume, nuclides, flux, micro_xs, chain, nuc_units='atom/cm3')

    fuel = Material(name="uo2")
    fuel.add_element("U", 1, percent_type="ao", enrichment=4.25)
    fuel.add_element("O", 2)
    fuel.set_density("g/cc", 10.4)
    fuel.depletable = True
    fuel.volume = 1
    materials = [fuel]
    fluxes = [1.0]
    micros = [micro_xs]
    IndependentOperator(materials, fluxes, micros, CHAIN_PATH)


def test_error_handling():
    micro_xs = MicroXS.from_csv(ONE_GROUP_XS)
    fuel = Material(name="oxygen")
    fuel.add_element("O", 2)
    fuel.set_density("g/cc", 1)
    fuel.depletable = True
    fuel.volume = 1
    materials = [fuel]
    fluxes = [1.0, 2.0]
    micros = [micro_xs]
    with pytest.raises(ValueError, match=r"The length of fluxes \(2\)"):
        IndependentOperator(materials, fluxes, micros, CHAIN_PATH)


def _uranium_material():
    fuel = Material(name="u metal")
    fuel.add_nuclide("U235", 0.05)
    fuel.add_nuclide("U238", 0.95)
    fuel.set_density("g/cc", 19.0)
    fuel.depletable = True
    fuel.volume = 1.0
    return fuel


def test_neutrons_emitted():
    assert _neutrons_emitted('fission') == 0
    assert _neutrons_emitted('(n,gamma)') == 0
    assert _neutrons_emitted('(n,p)') == 0
    assert _neutrons_emitted('(n,a)') == 0
    assert _neutrons_emitted('(n,3He)') == 0
    assert _neutrons_emitted('(n,2a)') == 0
    assert _neutrons_emitted('(n,np)') == 1
    assert _neutrons_emitted('(n,nd2a)') == 1
    assert _neutrons_emitted('(n,2n)') == 2
    assert _neutrons_emitted('(n,2nd)') == 2
    assert _neutrons_emitted('(n,3n)') == 3
    assert _neutrons_emitted('(n,3np)') == 3
    assert _neutrons_emitted('(n,4n)') == 4


def test_calculate_kinf():
    # Only U235 has nonzero cross sections so that the k-infinity estimate
    # does not depend on the material composition
    nuclides = ['U235', 'U238']
    reactions = ['fission', 'nu-fission', '(n,gamma)', '(n,2n)']
    data = np.array([
        [[50.0], [120.0], [10.0], [2.0]],
        [[0.0], [0.0], [0.0], [0.0]],
    ])
    micro_xs = MicroXS(data, nuclides, reactions)

    op = IndependentOperator(
        [_uranium_material()], [np.array([1.0])], [micro_xs], CHAIN_PATH,
        normalization_mode='source-rate', calculate_kinf=True)
    vec = op.initial_condition()

    # Production is nu-fission; loss is fission + (n,gamma) - (n,2n), since
    # (n,2n) produces one net neutron and is not counted as production
    expected = 120.0 / (50.0 + 10.0 - 2.0)
    result = op(vec, 1.0)
    assert result.k.n == pytest.approx(expected)
    assert result.k.s == 0.0

    # A decay step (zero source rate) reports the same estimate
    result = op(vec, 0.0)
    assert result.k.n == pytest.approx(expected)


def test_calculate_kinf_multigroup():
    # Two-group cross sections and flux, U235 only so that the estimate does
    # not depend on the material composition
    nuclides = ['U235']
    reactions = ['fission', 'nu-fission', '(n,gamma)']
    data = np.array([[[10.0, 50.0], [25.0, 120.0], [2.0, 10.0]]])
    micro_xs = MicroXS(data, nuclides, reactions)
    flux = np.array([0.75, 0.25])

    fuel = Material(name="u metal")
    fuel.add_nuclide("U235", 1.0)
    fuel.set_density("g/cc", 19.0)
    fuel.depletable = True
    fuel.volume = 1.0

    op = IndependentOperator(
        [fuel], [flux], [micro_xs], CHAIN_PATH,
        normalization_mode='source-rate', calculate_kinf=True)
    vec = op.initial_condition()

    production = 25.0*0.75 + 120.0*0.25
    loss = (10.0 + 2.0)*0.75 + (50.0 + 10.0)*0.25
    result = op(vec, 1.0)
    assert result.k.n == pytest.approx(production / loss)


def test_calculate_kinf_errors():
    # MicroXS without nu-fission data cannot be used to estimate k-infinity
    micro_xs = MicroXS.from_csv(ONE_GROUP_XS)
    with pytest.raises(ValueError, match="nu-fission"):
        IndependentOperator([_uranium_material()], [1.0], [micro_xs],
                            CHAIN_PATH, calculate_kinf=True)

    # keff and calculate_kinf are mutually exclusive
    data = np.zeros((1, 2, 1))
    micro_xs = MicroXS(data, ['U235'], ['fission', 'nu-fission'])
    with pytest.raises(ValueError, match="mutually exclusive"):
        IndependentOperator([_uranium_material()], [1.0], [micro_xs],
                            CHAIN_PATH, keff=(1.0, 0.0), calculate_kinf=True)
