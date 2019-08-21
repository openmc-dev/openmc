"""Test the FissionYieldHelpers"""

import os
from collections import namedtuple
from unittest.mock import Mock
import bisect

import pytest
import numpy
from openmc import capi
from openmc.deplete.nuclide import Nuclide, FissionYieldDistribution
from openmc.deplete.helpers import (
    FissionYieldCutoffHelper, ConstantFissionYieldHelper,
    AveragedFissionYieldHelper)


@pytest.fixture(scope="module")
def materials(tmpdir_factory):
    """Use C API to construct realistic materials for testing tallies"""
    tmpdir = tmpdir_factory.mktemp("capi")
    orig = tmpdir.chdir()
    # Create proxy xml files to please openmc
    with open("geometry.xml", "w") as stream:
        stream.write("""
<?xml version='1.0' encoding='utf-8'?>
<geometry>
  <cell id="1" material="1" name="fuel" region="1 -2 3 -4" universe="1" />
  <surface boundary="reflective" coeffs="-0.63" id="1" type="x-plane" />
  <surface boundary="reflective" coeffs="0.63" id="2" type="x-plane" />
  <surface boundary="reflective" coeffs="-0.63" id="3" type="y-plane" />
  <surface boundary="reflective" coeffs="0.63" id="4" type="y-plane" />
</geometry>
""")
    with open("settings.xml", "w") as stream:
        stream.write("""
<?xml version='1.0' encoding='utf-8'?>
<settings>
  <run_mode>eigenvalue</run_mode>
  <particles>100</particles>
  <batches>10</batches>
  <inactive>0</inactive>
  <verbosity>1</verbosity>
  <source strength="1.0">
    <space type="point">
      <parameters>0.0 0.0 0.0</parameters>
    </space>
  </source>
</settings>
""")
    with open("materials.xml", "w") as stream:
        stream.write("""
<?xml version='1.0' encoding='utf-8'?>
<materials>
  <material depletable="true" id="1" name="U" volume="71.67537585">
    <density units="g/cc" value="10.4" />
    <nuclide ao="1.0" name="U235" />
    <nuclide ao="1.0" name="U238" />
    <nuclide ao="1.0" name="Xe135" />
    <nuclide ao="1.0" name="Pu239" />
  </material>
</materials>""")
    try:
        with capi.run_in_memory():
            yield [capi.Material(), capi.Material()]
    finally:
        print(os.path.abspath(os.curdir))
        os.remove(tmpdir / "settings.xml")
        os.remove(tmpdir / "geometry.xml")
        os.remove(tmpdir / "materials.xml")
        os.remove(tmpdir / "summary.h5")
        orig.chdir()
        os.rmdir(tmpdir)


def proxy_tally_data(tally, fill=None):
    """Construct an empty matrix built from a C tally

    The shape of tally.results will be
    ``(n_bins, n_nuc * n_scores, 3)``
    """
    n_nucs = max(len(tally.nuclides), 1)
    n_scores = max(len(tally.scores), 1)
    n_bins = 1
    for tfilter in tally.filters:
        if not hasattr(tfilter, "bins"):
            continue
        this_bins = len(tfilter.bins)
        if isinstance(tfilter, capi.EnergyFilter):
            this_bins -= 1
        n_bins *= max(this_bins, 1)
    data = numpy.empty((n_bins, n_nucs * n_scores, 3))
    if fill is not None:
        data.fill(fill)
    return data


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

    pu239 = Nuclide("Pu239")
    pu239.yield_data = FissionYieldDistribution({
        0.0253: {"Xe135": 3.141e-3, "Sm149": 8.19e-10, "Gd155": 1.66e-9},
        5.0e5: {"Xe135": 6.14e-3, "Sm149": 9.429e-10, "Gd155": 5.24e-9},
        2e6: {"Xe135": 6.15e-3, "Sm149": 9.42e-10, "Gd155": 5.29e-9}})

    NuclideBundle = namedtuple("NuclideBundle", "u235 u238 xe135 pu239")
    return NuclideBundle(u235, u238, xe135, pu239)


@pytest.mark.parametrize(
    "input_energy, yield_energy",
    ((0.0253, 0.0253), (0.01, 0.0253), (4e5, 5e5)))
def test_constant_helper(nuclide_bundle, input_energy, yield_energy):
    helper = ConstantFissionYieldHelper(nuclide_bundle, energy=input_energy)
    assert helper.energy == input_energy
    assert helper.constant_yields == {
        "U235": nuclide_bundle.u235.yield_data[yield_energy],
        "U238": nuclide_bundle.u238.yield_data[5.00e5],  # only epithermal
        "Pu239": nuclide_bundle.pu239.yield_data[yield_energy]}
    assert helper.constant_yields == helper.weighted_yields(1)


def test_cutoff_construction(nuclide_bundle):
    u235 = nuclide_bundle.u235
    u238 = nuclide_bundle.u238
    pu239 = nuclide_bundle.pu239

    # defaults
    helper = FissionYieldCutoffHelper(nuclide_bundle, 1)
    assert helper.constant_yields == {
        "U238": u238.yield_data[5.0e5]}
    assert helper.thermal_yields == {
        "U235": u235.yield_data[0.0253],
        "Pu239": pu239.yield_data[0.0253]}
    assert helper.fast_yields == {
        "U235": u235.yield_data[5e5],
        "Pu239": pu239.yield_data[5e5]}

    # use 14 MeV yields
    helper = FissionYieldCutoffHelper(nuclide_bundle, 1, fast_energy=14e6)
    assert helper.constant_yields == {
        "U238": u238.yield_data[5.0e5]}
    assert helper.thermal_yields == {
        "U235": u235.yield_data[0.0253],
        "Pu239": pu239.yield_data[0.0253]}
    assert helper.fast_yields == {
        "U235": u235.yield_data[14e6],
        "Pu239": pu239.yield_data[2e6]}

    # specify missing thermal yields -> use 0.0253
    helper = FissionYieldCutoffHelper(nuclide_bundle, 1, thermal_energy=1)
    assert helper.thermal_yields == {
        "U235": u235.yield_data[0.0253],
        "Pu239": pu239.yield_data[0.0253]}
    assert helper.fast_yields == {
        "U235": u235.yield_data[5e5],
        "Pu239": pu239.yield_data[5e5]}

    # request missing fast yields -> use epithermal
    helper = FissionYieldCutoffHelper(nuclide_bundle, 1, fast_energy=1e4)
    assert helper.thermal_yields == {
        "U235": u235.yield_data[0.0253],
        "Pu239": pu239.yield_data[0.0253]}
    assert helper.fast_yields == {
        "U235": u235.yield_data[5e5],
        "Pu239": pu239.yield_data[5e5]}

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


# emulate some split between fast and thermal U235 fissions
@pytest.mark.parametrize("therm_frac", (0.5, 0.2, 0.8))
def test_cutoff_helper(materials, nuclide_bundle, therm_frac):
    helper = FissionYieldCutoffHelper(nuclide_bundle, len(materials))
    helper.generate_tallies(materials, [0])

    non_zero_nucs = [n.name for n in nuclide_bundle]
    tally_nucs = helper.update_tally_nuclides(non_zero_nucs)
    assert tally_nucs == ("Pu239", "U235",)

    # Check tallies
    fission_tally = helper._fission_rate_tally
    assert fission_tally is not None
    filters = fission_tally.filters
    assert len(filters) == 2
    assert isinstance(filters[0], capi.MaterialFilter)
    assert len(filters[0].bins) == len(materials)
    assert isinstance(filters[1], capi.EnergyFilter)
    # lower, cutoff, and upper energy
    assert len(filters[1].bins) == 3

    # Emulate building tallies
    # material x energy, tallied_nuclides, 3
    tally_data = proxy_tally_data(fission_tally)
    helper._fission_rate_tally = Mock()
    helper_flux = 1e6
    tally_data[0, :, 1] = therm_frac * helper_flux
    tally_data[1, :, 1] = (1 - therm_frac) * helper_flux
    helper._fission_rate_tally.results = tally_data

    helper.unpack()
    # expected results of shape (n_mats, 2, n_tnucs)
    expected_results = numpy.empty((1, 2, len(tally_nucs)))
    expected_results[:, 0] = therm_frac
    expected_results[:, 1] = 1 - therm_frac
    assert helper.results == pytest.approx(expected_results)

    actual_yields = helper.weighted_yields(0)
    assert actual_yields["U238"] == nuclide_bundle.u238.yield_data[5e5]
    for nuc in tally_nucs:
        assert actual_yields[nuc] == (
            helper.thermal_yields[nuc] * therm_frac
            + helper.fast_yields[nuc] * (1 - therm_frac))


@pytest.mark.parametrize("avg_energy", (0.01, 6e5, 15e6))
def test_averaged_helper(materials, nuclide_bundle, avg_energy):
    helper = AveragedFissionYieldHelper(nuclide_bundle)
    helper.generate_tallies(materials, [0])
    tallied_nucs = helper.update_tally_nuclides(
        [n.name for n in nuclide_bundle])
    assert tallied_nucs == ("Pu239", "U235")

    # check generated tallies
    fission_tally = helper._fission_rate_tally
    assert fission_tally is not None
    fission_filters = fission_tally.filters
    assert len(fission_filters) == 2
    assert isinstance(fission_filters[0], capi.MaterialFilter)
    assert len(fission_filters[0].bins) == len(materials)
    assert isinstance(fission_filters[1], capi.EnergyFilter)
    assert len(fission_filters[1].bins) == 2
    assert fission_tally.scores == ["fission"]
    assert fission_tally.nuclides == list(tallied_nucs)

    weighted_tally = helper._weighted_tally
    assert weighted_tally is not None
    weighted_filters = weighted_tally.filters
    assert len(weighted_filters) == 2
    assert isinstance(weighted_filters[0], capi.MaterialFilter)
    assert len(weighted_filters[0].bins) == len(materials)
    assert isinstance(weighted_filters[1], capi.EnergyFunctionFilter)
    assert len(weighted_filters[1].energy) == 2
    assert len(weighted_filters[1].y) == 2
    assert weighted_tally.scores == ["fission"]
    assert weighted_tally.nuclides == list(tallied_nucs)

    helper_flux = 1e16
    fission_results = proxy_tally_data(fission_tally, helper_flux)
    weighted_results = proxy_tally_data(
        weighted_tally, helper_flux * avg_energy)

    helper._fission_rate_tally = Mock()
    helper._weighted_tally = Mock()
    helper._fission_rate_tally.results = fission_results
    helper._weighted_tally.results = weighted_results

    helper.unpack()
    expected_results = numpy.ones((1, len(tallied_nucs))) * avg_energy
    assert helper.results == pytest.approx(expected_results)

    actual_yields = helper.weighted_yields(0)
    # constant U238 => no interpolation
    assert actual_yields["U238"] == nuclide_bundle.u238.yield_data[5e5]
    # construct expected yields
    exp_u235_yields = interp_average_yields(nuclide_bundle.u235, avg_energy)
    assert actual_yields["U235"] == exp_u235_yields
    exp_pu239_yields = interp_average_yields(nuclide_bundle.pu239, avg_energy)
    assert actual_yields["Pu239"] == exp_pu239_yields


def interp_average_yields(nuc, avg_energy):
    """Construct a set of yields by interpolation between neighbors"""
    energies = nuc.yield_energies
    yields = nuc.yield_data
    if avg_energy < energies[0]:
        return yields[energies[0]]
    if avg_energy > energies[-1]:
        return yields[energies[-1]]
    thermal_ix = bisect.bisect_left(energies, avg_energy)
    thermal_E, fast_E = energies[thermal_ix - 1:thermal_ix + 1]
    assert thermal_E < avg_energy < fast_E
    split = (avg_energy - thermal_E)/(fast_E - thermal_E)
    return yields[thermal_E]*(1 - split) + yields[fast_E]*split
