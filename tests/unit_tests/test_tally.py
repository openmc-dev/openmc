import openmc


def test_tally_init_args():
    """Test that Tally constructor kwargs are applied correctly."""
    filter = openmc.EnergyFilter([0.0, 1.0, 20.0e6])
    tally = openmc.Tally(
        name='my tally',
        scores=['flux', 'fission'],
        filters=[filter],
        nuclides=['U235'],
        estimator='tracklength',
    )

    assert tally.name == 'my tally'
    assert tally.scores == ['flux', 'fission']
    assert tally.filters == [filter]
    assert tally.nuclides == ['U235']
    assert tally.estimator == 'tracklength'
