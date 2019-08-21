"""Test the FissionYieldHelpers"""

from collections import namedtuple
from unittest.mock import Mock

import pytest
import numpy
from openmc.deplete.nuclide import Nuclide, FissionYieldDistribution
from openmc.deplete.helpers import (
    FissionYieldCutoffHelper, ConstantFissionYieldHelper,
    AveragedFissionYieldHelper)


MATERIALS = ["1", "2"]


@pytest.fixture(scope="module")
def nuclide_bundle():
    u5yield_dict = {
        0.0253: {"Xe135": 7.85e-4, "Gd155": 4.08e-12, "Sm149": 1.71e-12},
        5.0e5: {"Xe135": 7.85e-4, "Sm149": 1.71e-12},
        1.40e7: {"Xe135": 4.54e-3, "Gd155": 5.83e-8}}
    u235 = Nuclide("U235")
    u235.yield_data = FissionYieldDistribution(u5yield_dict)

    u8yield_dict = {5.00e5: {"Xe135": 1.12e-3, "Gd155": 1.32e-12}}
    u238 = Nuclide("U238")
    u238.yield_data = FissionYieldDistribution(u8yield_dict)

    xe135 = Nuclide("Xe135")

    NuclideBundle = namedtuple("NuclideBundle", "u235 u238 xe135")
    return NuclideBundle(u235, u238, xe135)


@pytest.mark.parametrize(
    "input_energy, u5_yield_energy",
    ((0.0253, 0.0253), (0.01, 0.0253), (4e5, 5e5)))
def test_constant_helper(nuclide_bundle, input_energy, u5_yield_energy):
    helper = ConstantFissionYieldHelper(nuclide_bundle, energy=input_energy)
    assert helper.energy == input_energy
    assert helper.constant_yields == {
        "U235": nuclide_bundle.u235.yield_data[u5_yield_energy],
        "U238": nuclide_bundle.u238.yield_data[5.00e5]}  # only epithermal
    assert helper.constant_yields == helper.weighted_yields(1)


def test_cutoff_construction(nuclide_bundle):
    # defaults
    helper = FissionYieldCutoffHelper(nuclide_bundle, 1)
    assert helper.constant_yields == {
        "U238": nuclide_bundle.u238.yield_data[5.0e5]}
    assert helper.thermal_yields == {
        "U235": nuclide_bundle.u235.yield_data[0.0253]}
    assert helper.fast_yields == {"U235": nuclide_bundle.u235.yield_data[5e5]}
    # use 14 MeV yields
    helper = FissionYieldCutoffHelper(nuclide_bundle, 1, fast_energy=14e6)
    assert helper.constant_yields == {
        "U238": nuclide_bundle.u238.yield_data[5.0e5]}
    assert helper.thermal_yields == {
        "U235": nuclide_bundle.u235.yield_data[0.0253]}
    assert helper.fast_yields == {"U235": nuclide_bundle.u235.yield_data[14e6]}
    # specify missing thermal yields -> use 0.0253
    helper = FissionYieldCutoffHelper(nuclide_bundle, 1, thermal_energy=1)
    assert helper.thermal_yields == {
        "U235": nuclide_bundle.u235.yield_data[0.0253]}
    assert helper.fast_yields == {"U235": nuclide_bundle.u235.yield_data[5e5]}
    # request missing fast yields -> use epithermal
    helper = FissionYieldCutoffHelper(nuclide_bundle, 1, fast_energy=1e4)
    assert helper.thermal_yields == {
        "U235": nuclide_bundle.u235.yield_data[0.0253]}
    assert helper.fast_yields == {"U235": nuclide_bundle.u235.yield_data[5e5]}
    # test failures in cutoff: super low, super high
    with pytest.raises(ValueError, match="replacement fission yields"):
        FissionYieldCutoffHelper(
            nuclide_bundle, 1, thermal_energy=0.001, cutoff=0.002)
    with pytest.raises(ValueError, match="replacement fission yields"):
        FissionYieldCutoffHelper(
            nuclide_bundle, 1, cutoff=15e6, fast_energy=17e6)


@pytest.mark.parametrize("key", ("cutoff", "thermal_energy", "fast_energy"))
def test_cutoff_failure(key):
    with pytest.raises(TypeError, match=key):
        FissionYieldCutoffHelper(None, None, **{key: None})
    with pytest.raises(ValueError, match=key):
        FissionYieldCutoffHelper(None, None, **{key: -1})


class ProxyMixin:
    """Mixing that overloads the tally generation"""
    def generate_tallies(self, materials, mat_indexes):
        self._fission_rate_tally = Mock()
        self._local_indexes = numpy.asarray(mat_indexes)


class CutoffProxy(ProxyMixin, FissionYieldCutoffHelper):
    """Proxy that supplies a set of tallies"""


# emulate some split between fast and thermal U235 fissions
@pytest.mark.parametrize("therm_frac", (0.5, 0.2, 0.8))
def test_cutoff_helper(nuclide_bundle, therm_frac):
    n_bmats = len(MATERIALS)
    proxy = CutoffProxy(nuclide_bundle, n_bmats)
    proxy.generate_tallies(MATERIALS, [0])
    non_zero_nucs = [n.name for n in nuclide_bundle]
    tally_nucs = proxy.update_tally_nuclides(non_zero_nucs)
    assert tally_nucs == ("U235",)
    # Emulate building tallies
    # material x energy, tallied_nuclides, 3
    proxy_flux = 1e6
    tally_data = numpy.empty((n_bmats * 2, 1, 3))
    tally_data[0, 0, 1] = therm_frac * proxy_flux
    tally_data[1, 0, 1] = (1 - therm_frac) * proxy_flux
    proxy._fission_rate_tally.results = tally_data

    proxy.unpack()
    # expected results of shape (n_mats, 2, n_tnucs)
    expected_results = numpy.empty((1, 2, 1))
    expected_results[:, 0] = therm_frac
    expected_results[:, 1] = 1 - therm_frac
    assert proxy.results == pytest.approx(expected_results)

    actual_yields = proxy.weighted_yields(0)
    assert actual_yields["U238"] == nuclide_bundle.u238.yield_data[5e5]
    assert actual_yields["U235"] == (
        proxy.thermal_yields["U235"] * therm_frac
        + proxy.fast_yields["U235"] * (1 - therm_frac))


class AverageProxy(ProxyMixin, AveragedFissionYieldHelper):
    """Proxy for generating mock set of tallies"""
    def generate_tallies(self, materials, mat_indexes):
        super().generate_tallies(materials, mat_indexes)
        self._weighted_tally = Mock()


@pytest.mark.parametrize("avg_energy", (0.01, 100, 15e6))
def test_averaged_helper(nuclide_bundle, avg_energy):
    proxy = AverageProxy(nuclide_bundle)
    proxy.generate_tallies(MATERIALS, [0])
    tallied_nucs = proxy.update_tally_nuclides(
        [n.name for n in nuclide_bundle])
    assert tallied_nucs == ("U235", )
    # enforce some average energy
    proxy_flux = 1e16
    fission_results = numpy.ones((len(MATERIALS), 1, 3)) * proxy_flux
    weighted_results = fission_results * avg_energy
    proxy._fission_rate_tally.results = fission_results
    proxy._weighted_tally.results = weighted_results
    proxy.unpack()
    expected_results = numpy.ones((1, 1)) * avg_energy
    assert proxy.results == pytest.approx(expected_results)

    actual_yields = proxy.weighted_yields(0)
    # constant U238 => no interpolation
    assert actual_yields["U238"] == nuclide_bundle.u238.yield_data[5e5]
    # construct expected yields
    if avg_energy < 0.0253:  # take thermal U235 yields
        exp_u235_yields = nuclide_bundle.u235.yield_data[0.0253]
    elif avg_energy > 14e6:  # take fastest U235 yields
        exp_u235_yields = nuclide_bundle.u235.yield_data[14e6]
    else:  # reconstruct between thermal and epithermal
        thermal = nuclide_bundle.u235.yield_data[0.0253]
        epithermal = nuclide_bundle.u235.yield_data[5e5]
        split = (avg_energy - 0.0253) / (5e5 - 0.0253)
        exp_u235_yields = thermal * (1 - split) + epithermal * split
    assert actual_yields["U235"] == exp_u235_yields
