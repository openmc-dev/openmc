"""Basic unit tests for openmc.deplete.IndependentOperator instantiation

"""

from pathlib import Path

from openmc.deplete import IndependentOperator, MicroXS
from openmc import Material, Materials

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
    IndependentOperator.from_nuclides(
        volume, nuclides, flux, micro_xs, CHAIN_PATH, nuc_units='atom/cm3')

    fuel = Material(name="uo2")
    fuel.add_element("U", 1, percent_type="ao", enrichment=4.25)
    fuel.add_element("O", 2)
    fuel.set_density("g/cc", 10.4)
    fuel.depletable = True
    fuel.volume = 1
    materials = Materials([fuel])
    fluxes = [1.0]
    micros = [micro_xs]
    IndependentOperator(materials, fluxes, micros, CHAIN_PATH)
