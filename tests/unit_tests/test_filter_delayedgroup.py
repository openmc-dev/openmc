"""Test that delayed group filter bins follow the groups they name rather than
their position in the filter's list of bins."""

import pytest
import openmc


GROUPS_ASCENDING = [1, 2, 3, 4, 5, 6]
GROUPS_SHUFFLED = [4, 1, 6, 3, 5, 2]
GROUPS_SUBSET = [2, 4]

SCORES = ["delayed-nu-fission", "ifp-beta-numerator"]


@pytest.fixture(scope="module")
def geometry():
    openmc.reset_auto_ids()
    material = openmc.Material()
    material.add_nuclide("U235", 1.0)
    material.set_density("g/cm3", 16.0)
    sphere = openmc.Sphere(r=10.0, boundary_type="vacuum")
    cell = openmc.Cell(region=-sphere, fill=material)
    return openmc.Geometry([cell])


def build_model(geometry, score):
    """Model tallying one score four ways: three delayed group filters that
    differ only in how the bins are listed, plus an unfiltered reference."""
    model = openmc.Model(geometry=geometry)
    model.settings.particles = 5000
    model.settings.batches = 20
    model.settings.inactive = 5
    model.settings.ifp_n_generation = 5

    for name, groups in [
        ("ascending", GROUPS_ASCENDING),
        ("shuffled", GROUPS_SHUFFLED),
        ("subset", GROUPS_SUBSET),
    ]:
        tally = openmc.Tally(name=name)
        tally.scores = [score]
        tally.filters = [openmc.DelayedGroupFilter(groups)]
        model.tallies.append(tally)

    tally = openmc.Tally(name="unfiltered")
    tally.scores = [score]
    model.tallies.append(tally)

    return model


def means_by_group(statepoint, name):
    """Map delayed group number to tallied mean.

    The group order is read back from the filter stored in the statepoint, so
    a mismatch between the reported bin labels and the values fails here
    instead of being hidden by assuming an order.
    """
    tally = statepoint.get_tally(name=name)
    delayed_group_filter = tally.find_filter(openmc.DelayedGroupFilter)
    means = tally.mean.ravel()
    assert len(means) == len(delayed_group_filter.bins)
    return dict(zip(delayed_group_filter.bins, means))


@pytest.mark.parametrize("score", SCORES)
def test_delayed_group_bins(run_in_tmpdir, geometry, score):
    """Delayed group labels determine values independently of bin order.

    All tallies come from a single run, so they see the same histories and the
    same scoring events. Corresponding groups must therefore agree within
    statistics, and all selected groups must sum to the unfiltered result.
    """
    model = build_model(geometry, score)
    sp_file = model.run()
    with openmc.StatePoint(sp_file) as sp:
        ascending = means_by_group(sp, "ascending")
        shuffled = means_by_group(sp, "shuffled")
        subset = means_by_group(sp, "subset")
        unfiltered = sp.get_tally(name="unfiltered").mean.ravel()[0]

    # Guard against a vacuous pass if the run were too short to populate bins.
    assert all(value > 0.0 for value in ascending.values())

    # Listing groups out of order must not move the corresponding results.
    assert sorted(shuffled) == sorted(ascending)
    for group in GROUPS_ASCENDING:
        assert shuffled[group] == pytest.approx(ascending[group], rel=1e-6)

    # A partial selection must reproduce those groups from the full list. A
    # group-to-bin index assumption would write out of bounds for group 4.
    assert sorted(subset) == sorted(GROUPS_SUBSET)
    for group in GROUPS_SUBSET:
        assert subset[group] == pytest.approx(ascending[group], rel=1e-6)

    # Every delayed neutron lands in exactly one bin and none are lost.
    assert sum(ascending.values()) == pytest.approx(unfiltered, rel=1e-6)
