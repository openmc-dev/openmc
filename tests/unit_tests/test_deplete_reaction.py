"""Tests for the openmc.deplete.ReactionRates class."""

from openmc.deplete import reaction_rates


def test_indexing():
    """Tests the __getitem__ and __setitem__ routines simultaneously."""

    mat_to_ind = {"10000" : 0, "10001" : 1}
    nuc_to_ind = {"U238" : 0, "U235" : 1}
    react_to_ind = {"fission" : 0, "(n,gamma)" : 1}

    rates = reaction_rates.ReactionRates(mat_to_ind, nuc_to_ind, react_to_ind)

    rates["10000", "U238", "fission"] = 1.0
    rates["10001", "U238", "fission"] = 2.0
    rates["10000", "U235", "fission"] = 3.0
    rates["10001", "U235", "fission"] = 4.0
    rates["10000", "U238", "(n,gamma)"] = 5.0
    rates["10001", "U238", "(n,gamma)"] = 6.0
    rates["10000", "U235", "(n,gamma)"] = 7.0
    rates["10001", "U235", "(n,gamma)"] = 8.0

    # String indexing
    assert rates["10000", "U238", "fission"] == 1.0
    assert rates["10001", "U238", "fission"] == 2.0
    assert rates["10000", "U235", "fission"] == 3.0
    assert rates["10001", "U235", "fission"] == 4.0
    assert rates["10000", "U238", "(n,gamma)"] == 5.0
    assert rates["10001", "U238", "(n,gamma)"] == 6.0
    assert rates["10000", "U235", "(n,gamma)"] == 7.0
    assert rates["10001", "U235", "(n,gamma)"] == 8.0

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
    assert rates["10000", "U238", "fission"] == 5.0


def test_n_mat():
    """Test number of materials property."""
    mat_to_ind = {"10000" : 0, "10001" : 1}
    nuc_to_ind = {"U238" : 0, "U235" : 1, "Gd157" : 2}
    react_to_ind = {"fission" : 0, "(n,gamma)" : 1, "(n,2n)" : 2, "(n,3n)" : 3}

    rates = reaction_rates.ReactionRates(mat_to_ind, nuc_to_ind, react_to_ind)

    assert rates.n_mat == 2


def test_n_nuc():
    """Test number of nuclides property."""
    mat_to_ind = {"10000" : 0, "10001" : 1}
    nuc_to_ind = {"U238" : 0, "U235" : 1, "Gd157" : 2}
    react_to_ind = {"fission" : 0, "(n,gamma)" : 1, "(n,2n)" : 2, "(n,3n)" : 3}

    rates = reaction_rates.ReactionRates(mat_to_ind, nuc_to_ind, react_to_ind)

    assert rates.n_nuc == 3


def test_n_react():
    """ Test number of reactions property. """
    mat_to_ind = {"10000" : 0, "10001" : 1}
    nuc_to_ind = {"U238" : 0, "U235" : 1, "Gd157" : 2}
    react_to_ind = {"fission" : 0, "(n,gamma)" : 1, "(n,2n)" : 2, "(n,3n)" : 3}

    rates = reaction_rates.ReactionRates(mat_to_ind, nuc_to_ind, react_to_ind)

    assert rates.n_react == 4
