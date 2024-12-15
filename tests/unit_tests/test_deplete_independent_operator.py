"""Basic unit tests for openmc.deplete.IndependentOperator instantiation

"""

from pathlib import Path

import pytest

from openmc import Material
from openmc.deplete import IndependentOperator, MicroXS

CHAIN_PATH = Path(__file__).parents[1] / "chain_simple.xml"
ONE_GROUP_XS = Path(__file__).parents[1] / "micro_xs_simple.csv"


def test_operator_init():
    """The test uses a temporary dummy chain. This file will be removed
    at the end of the test, and only contains a depletion_chain node."""
    volume = 1
    nuclides = {
        "U234": 8.922411359424315e18,
        "U235": 9.98240191860822e20,
        "U238": 2.2192386373095893e22,
        "U236": 4.5724195495061115e18,
        "O16": 4.639065406771322e22,
        "O17": 1.7588724018066158e19,
    }
    flux = 1.0
    micro_xs = MicroXS.from_csv(ONE_GROUP_XS)
    IndependentOperator.from_nuclides(
        volume, nuclides, flux, micro_xs, CHAIN_PATH, nuc_units="atom/cm3"
    )

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
