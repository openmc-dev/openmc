"""Test the FissionYieldHelpers"""

from collections import namedtuple

import pytest

from openmc.deplete.nuclide import Nuclide, FissionYieldDistribution
from openmc.deplete.helpers import ConstantFissionYieldHelper


@pytest.fixture(scope="module")
def nuclide_bundle():
    u5yield_dict = {
        0.0253: {"Xe135": 7.85e-4, "Gd155": 4.08e-12, "Sm149": 1.71e-12},
        5.0e5: {"Xe135": 7.85e-4, "Sm149": 1.71e-12},
        1.40e7: {"Xe135": 4.54e-3, "Gd155": 5.83e-8}}
    u235 = Nuclide()
    u235.name = "U235"
    u235.yield_energies = (0.0253, 5.0e5, 1.40e7)
    u235.yield_data = FissionYieldDistribution.from_dict(u5yield_dict)

    u8yield_dict = {5.00e5: {"Xe135": 1.12e-3, "Gd155": 1.32e-12}}
    u238 = Nuclide()
    u238.name = "U238"
    u238.yield_energies = (5.00e5, )
    u238.yield_data = FissionYieldDistribution.from_dict(u8yield_dict)

    xe135 = Nuclide()
    xe135.name = "Xe135"

    NuclideBundle = namedtuple("NuclideBundle", "u235 u238 xe135")
    return NuclideBundle(u235, u238, xe135)
@pytest.mark.parametrize("input_energy, u5_yield_energy", (
    (0.0253, 0.0253), (0.01, 0.0253), (4e5, 5e5)))
def test_constant_helper(nuclide_bundle, input_energy, u5_yield_energy):
    helper = ConstantFissionYieldHelper(
        nuclide_bundle, energy=input_energy)
    assert helper.energy == input_energy
    assert helper.constant_yields == {
        "U235": nuclide_bundle.u235.yield_data[u5_yield_energy],
        "U238": nuclide_bundle.u238.yield_data[5.00e5]}  # only epithermal
    assert helper.constant_yields == helper.weighted_yields(1)
