"""Test the Operator helpers"""

from unittest.mock import Mock

import pytest
import numpy
from numpy.testing import assert_array_equal

from openmc.deplete.nuclide import Nuclide, FissionYieldDistribution
from openmc.deplete.helpers import FissionYieldHelper


class FissionYieldHelperProxy(FissionYieldHelper):

    def __init__(self, chain_nuclides, n_bmats):
        super().__init__(chain_nuclides, n_bmats)
        self._reaction_tally = Mock()
        self.local_indexes = numpy.array([0])

    def generate_tallies(self, *args, **kwargs):
        # Avoid calls to the C-API
        pass


def test_fission_yield_helper():
    """Test the collection of fission yield data using approximated tallies
    """
    u5yield_dict = {
        0.0253: {"Xe135": 7.85e-4, "Gd155": 4.08e-12, "Sm149": 1.71e-12},
        1.40e7: {"Xe135": 4.54e-3, "Gd155": 5.83e-8}}
    u235 = Nuclide()
    u235.name = "U235"
    u235.yield_energies = (0.0253, 1.40e7)
    u235.yield_data = FissionYieldDistribution.from_dict(u5yield_dict)

    u8yield_dict = {5.00e5: {"Xe135": 1.12e-3, "Gd155": 1.32e-12}}
    u238 = Nuclide()
    u238.name = "U238"
    u238.yield_energies = (5.00e5, )
    u238.yield_data = FissionYieldDistribution.from_dict(u8yield_dict)

    xe135 = Nuclide()
    xe135.name = "Xe135"

    n_bmats = 2
    helper = FissionYieldHelperProxy([u235, u238, xe135], n_bmats)

    assert helper.energy_bounds == (0, 0.0253, 5.00e5, 1.40e7)

    # test that tally must be created with nuclides with yields
    with pytest.raises(ValueError, match="No overlap"):
        helper.set_fissionable_nuclides(["Xe135", ])

    act_nucs = helper.set_fissionable_nuclides(["U235", "U238", "Xe135"])
    assert act_nucs == ("U235", "U238")

    # Emulate getting tally data from transport run
    # Test as if this Helper is responsible for one of two materials
    # Tally results ordered [n_mat * n_ene, n_fiss_nuc, 3]

    tally_res = numpy.zeros((3 * n_bmats, 2, 3))
    u5_fiss_rates = numpy.array([1.0, 1.5, 2.0])
    u8_fiss_rates = numpy.array([0.0, 1.0, 1.0])
    tally_res[:3, 0, 1] = u5_fiss_rates
    tally_res[:3, 1, 1] = u8_fiss_rates

    helper._reaction_tally.results = tally_res
    helper.unpack()  # compute yield fractions

    # Compare fraction fission rate from helper.reset
    exp_results = numpy.empty((1, 3, 2))
    # Fraction of fission events in each energy range
    u5_vec = u5_fiss_rates / u5_fiss_rates.sum()
    u8_vec = u8_fiss_rates / u8_fiss_rates.sum()
    exp_results[..., 0] = u5_vec
    exp_results[..., 1] = u8_vec
    assert_array_equal(helper.results, exp_results)

    # Compute and compare the library of fission yields
    exp_lib = {
        "U235": {
            "Xe135": (u5yield_dict[0.0253]["Xe135"] * u5_vec[0]
                      + u5yield_dict[1.40e7]["Xe135"] * u5_vec[2]),
            "Gd155": (u5yield_dict[0.0253]["Gd155"] * u5_vec[0]
                      + u5yield_dict[1.40e7]["Gd155"] * u5_vec[2]),
            "Sm149": u5yield_dict[0.0253]["Sm149"] * u5_vec[0],
        },
        "U238": {
            "Xe135": u8yield_dict[5.00e5]["Xe135"] * u8_vec[1],
            "Gd155": u8yield_dict[5.00e5]["Gd155"] * u8_vec[1],
        }
    }

    assert len(helper.libraries) == 0
    helper.compute_yields(0)
    assert len(helper.libraries) == 1
    act_library = helper.libraries[0]
    for parent, sub_yields in exp_lib.items():
        assert act_library[parent] == pytest.approx(sub_yields)
