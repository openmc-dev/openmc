"""Tests for the openmc.deplete.ReactionRates class."""

import numpy as np
from openmc.deplete import ReactionRates


def test_get_set():
    """Tests the get/set methods."""

    local_mats = ["10000", "10001"]
    nuclides = ["U238", "U235"]
    reactions = ["fission", "(n,gamma)"]

    rates = ReactionRates(local_mats, nuclides, reactions)
    assert rates.shape == (2, 2, 2)
    assert np.all(rates == 0.0)

    rates.set("10000", "U238", "fission", 1.0)
    rates.set("10001", "U238", "fission", 2.0)
    rates.set("10000", "U235", "fission", 3.0)
    rates.set("10001", "U235", "fission", 4.0)
    rates.set("10000", "U238", "(n,gamma)", 5.0)
    rates.set("10001", "U238", "(n,gamma)", 6.0)
    rates.set("10000", "U235", "(n,gamma)", 7.0)
    rates.set("10001", "U235", "(n,gamma)", 8.0)

    # String indexing
    assert rates.get("10000", "U238", "fission") == 1.0
    assert rates.get("10001", "U238", "fission") == 2.0
    assert rates.get("10000", "U235", "fission") == 3.0
    assert rates.get("10001", "U235", "fission") == 4.0
    assert rates.get("10000", "U238", "(n,gamma)") == 5.0
    assert rates.get("10001", "U238", "(n,gamma)") == 6.0
    assert rates.get("10000", "U235", "(n,gamma)") == 7.0
    assert rates.get("10001", "U235", "(n,gamma)") == 8.0

    # Int indexing
    assert rates[0, 0, 0] == 1.0
    assert rates[1, 0, 0] == 2.0
    assert rates[0, 1, 0] == 3.0
    assert rates[1, 1, 0] == 4.0
    assert rates[0, 0, 1] == 5.0
    assert rates[1, 0, 1] == 6.0
    assert rates[0, 1, 1] == 7.0
    assert rates[1, 1, 1] == 8.0

    rates[0, 0, 0] = 5.0

    assert rates[0, 0, 0] == 5.0
    assert rates.get("10000", "U238", "fission") == 5.0


def test_properties():
    """Test number of materials property."""
    local_mats = ["10000", "10001"]
    nuclides = ["U238", "U235", "Gd157"]
    reactions = ["fission", "(n,gamma)", "(n,2n)", "(n,3n)"]

    rates = ReactionRates(local_mats, nuclides, reactions)

    assert rates.n_mat == 2
    assert rates.n_nuc == 3
    assert rates.n_react == 4
