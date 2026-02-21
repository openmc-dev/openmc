"""Basic unit tests for openmc.deplete.IndependentOperator instantiation

Modifies and resets environment variable OPENMC_CROSS_SECTIONS
to a custom file with new depletion_chain node
"""

from os import remove
from pathlib import Path
from unittest.mock import patch

import pytest
import openmc
from openmc.deplete import MicroXS, get_microxs_and_flux
import numpy as np

ONE_GROUP_XS = Path(__file__).parents[1] / "micro_xs_simple.csv"
CHAIN_FILE = Path(__file__).parents[1] / "chain_simple.xml"


def test_from_array():
    nuclides = [
        'U234',
        'U235',
        'U238',
        'U236',
        'O16',
        'O17',
        'I135',
        'Xe135',
        'Xe136',
        'Cs135',
        'Gd157',
        'Gd156']
    reactions = ['fission', '(n,gamma)']
    # These values are placeholders and are not at all
    # physically meaningful.
    data = np.array([[0.1, 0.],
                     [0.1, 0.],
                     [0.9, 0.],
                     [0.4, 0.],
                     [0., 0.],
                     [0., 0.],
                     [0., 0.1],
                     [0., 0.9],
                     [0., 0.],
                     [0., 0.],
                     [0., 0.1],
                     [0., 0.1]])
    data.shape = (12, 2, 1)

    MicroXS(data, nuclides, reactions)
    with pytest.raises(ValueError, match='Data array must be 3D'):
        MicroXS(data[:, 0], nuclides, reactions)


def test_csv():
    ref_xs = MicroXS.from_csv(ONE_GROUP_XS)
    ref_xs.to_csv('temp_xs.csv')
    temp_xs = MicroXS.from_csv('temp_xs.csv')
    assert np.all(ref_xs.data == temp_xs.data)
    remove('temp_xs.csv')


def test_from_multigroup_flux():
    energies = [0., 6.25e-1, 5.53e3, 8.21e5, 2.e7]
    flux = [1.1e-7, 1.2e-6, 1.3e-5, 1.4e-4]
    chain_file = Path(__file__).parents[1] / 'chain_simple.xml'
    kwargs = {'multigroup_flux': flux, 'chain_file': chain_file}

    # test with energy group structure from string
    microxs = MicroXS.from_multigroup_flux(energies='CASMO-4', **kwargs)
    assert isinstance(microxs, MicroXS)

    # test with energy group structure as floats
    microxs = MicroXS.from_multigroup_flux(energies=energies, **kwargs)
    assert isinstance(microxs, MicroXS)

    # test with nuclides provided
    microxs = MicroXS.from_multigroup_flux(
        energies=energies, nuclides=['Gd157', 'H1'], **kwargs
    )
    assert isinstance(microxs, MicroXS)
    assert microxs.nuclides == ['Gd157', 'H1']

    # test with reactions provided
    microxs = MicroXS.from_multigroup_flux(
        energies=energies, reactions=['fission', '(n,2n)'], **kwargs
    )
    assert isinstance(microxs, MicroXS)
    assert microxs.reactions == ['fission', '(n,2n)']


def test_multigroup_flux_same():
    chain_file = Path(__file__).parents[1] / 'chain_simple.xml'

    # Generate micro XS based on 4-group flux
    energies = [0., 6.25e-1, 5.53e3, 8.21e5, 2.e7]
    flux_per_ev = [0.3, 0.3, 1.0, 1.0]
    flux = flux_per_ev * np.diff(energies)
    flux_sum = flux.sum()
    microxs_4g = MicroXS.from_multigroup_flux(
        energies=energies, multigroup_flux=flux, chain_file=chain_file)

    # from_multigroup_flux should not modify the flux
    assert flux.sum() == flux_sum

    # Generate micro XS based on 2-group flux, where the boundaries line up with
    # the 4 group flux and have the same flux per eV across the full energy
    # range
    energies = [0., 5.53e3, 2.0e7]
    flux_per_ev = [0.3, 1.0]
    flux = flux_per_ev * np.diff(energies)
    microxs_2g = MicroXS.from_multigroup_flux(
        energies=energies, multigroup_flux=flux, chain_file=chain_file)

    assert microxs_4g.data == pytest.approx(microxs_2g.data)


def test_microxs_zero_flux():
    chain_file = Path(__file__).parents[1] / 'chain_simple.xml'

    # Generate micro XS based on zero flux
    energies = [0., 6.25e-1, 5.53e3, 8.21e5, 2.e7]
    flux = [0.0, 0.0, 0.0, 0.0]
    microxs = MicroXS.from_multigroup_flux(
        energies=energies, multigroup_flux=flux, chain_file=chain_file)

    # All microscopic cross sections should be zero
    assert np.all(microxs.data == 0.0)


def test_hybrid_tally_setup():
    """In hybrid mode a 1-group RR tally is added alongside the flux tally."""
    # Create a simple model with one material and a few nuclides for testing
    model = openmc.Model()
    mat = openmc.Material(components={'U235': 1.0, 'O16': 2.0})
    sphere = openmc.Sphere(r=10.0, boundary_type='vacuum')
    cell = openmc.Cell(region=-sphere, fill=mat)
    model.geometry = openmc.Geometry([cell])
    model.settings.batches = 2
    model.settings.particles = 10

    # Define 2-group energy structure for the test
    energies = [0., 0.625, 2.0e7]

    # Function to replace Model.run and capture the tallies that were created
    captured = {}
    def capture_run(**kwargs):
        captured['tallies'] = list(model.tallies)
        raise StopIteration

    # Call get_microxs_and_flux but replace Model.run with a function that
    # captures the tallies and raises StopIteration to exit early
    with patch.object(model, 'run', side_effect=capture_run):
        with pytest.raises(StopIteration):
            get_microxs_and_flux(
                model, [mat],
                nuclides=['U235', 'O16'],
                reactions=['fission', '(n,gamma)'],
                energies=energies,
                reaction_rate_mode='flux',
                reaction_rate_opts={'nuclides': ['U235'], 'reactions': ['fission']},
                chain_file=CHAIN_FILE,
            )

    # Check that both tallies were created with the expected properties
    tally_names = [t.name for t in captured['tallies']]
    assert 'MicroXS flux' in tally_names
    assert 'MicroXS RR' in tally_names

    # Check that the RR tally has the expected nuclides and reactions
    rr = next(t for t in captured['tallies'] if t.name == 'MicroXS RR')
    assert rr.nuclides == ['U235']
    assert rr.scores == ['fission']

    # RR tally must use a 1-group energy filter spanning the full energy range
    ef = next(f for f in rr.filters if isinstance(f, openmc.EnergyFilter))
    assert len(ef.values) == 2
    assert ef.values[0] == pytest.approx(energies[0])
    assert ef.values[-1] == pytest.approx(energies[-1])

# ---------------------------------------------------------------------------
# Tests for MicroXS.merge()
# ---------------------------------------------------------------------------

def _make_microxs(nuclides, reactions, values, groups=1):
    """Helper: build a MicroXS from a flat list of values (nuclide-major order)."""
    data = np.array(values, dtype=float).reshape(
        len(nuclides), len(reactions), groups)
    return MicroXS(data, nuclides, reactions)


def test_merge_disjoint():
    """Merging two MicroXS with no overlapping nuclides or reactions."""
    m1 = _make_microxs(['U235', 'U238'], ['fission', '(n,gamma)'],
                        [1., 2., 3., 4.])
    m2 = _make_microxs(['Pu239'], ['(n,2n)'], [5.])

    merged = m1.merge(m2)

    assert merged.nuclides == ['U235', 'U238', 'Pu239']
    assert merged.reactions == ['fission', '(n,gamma)', '(n,2n)']
    assert merged.data.shape == (3, 3, 1)

    # Self data preserved
    assert merged['U235', 'fission'] == pytest.approx([1.])
    assert merged['U238', '(n,gamma)'] == pytest.approx([4.])
    # New data from other
    assert merged['Pu239', '(n,2n)'] == pytest.approx([5.])
    # Cross-terms that had no data should be zero
    assert merged['U235', '(n,2n)'] == pytest.approx([0.])
    assert merged['Pu239', 'fission'] == pytest.approx([0.])


def test_merge_prefer_other():
    """prefer='other': other overwrites shared entries, adds new ones."""
    # m1: U235 and U238, reactions fission and (n,gamma)
    m1 = _make_microxs(['U235', 'U238'], ['fission', '(n,gamma)'],
                        [1., 2., 3., 4.])
    # m2: only U235, reactions fission (overlap) and (n,2n) (new)
    m2 = _make_microxs(['U235'], ['fission', '(n,2n)'], [9., 5.])

    merged = m1.merge(m2)

    # Nuclide/reaction sets
    assert set(merged.nuclides) == {'U235', 'U238'}
    assert set(merged.reactions) == {'fission', '(n,gamma)', '(n,2n)'}

    # U235/fission overwritten by m2
    assert merged['U235', 'fission'] == pytest.approx([9.])
    # U235/(n,gamma) untouched
    assert merged['U235', '(n,gamma)'] == pytest.approx([2.])
    # New reaction added from m2
    assert merged['U235', '(n,2n)'] == pytest.approx([5.])
    # U238 data preserved
    assert merged['U238', 'fission'] == pytest.approx([3.])
    assert merged['U238', '(n,gamma)'] == pytest.approx([4.])
    # U238/(n,2n) not in either → zero
    assert merged['U238', '(n,2n)'] == pytest.approx([0.])


def test_merge_prefer_self():
    """prefer='self': shared pairs keep self's value; new entries use other's value."""
    m1 = _make_microxs(['U235', 'U238'], ['fission', '(n,gamma)'],
                        [1., 2., 3., 4.])
    m2 = _make_microxs(['U235'], ['fission', '(n,2n)'], [9., 5.])

    merged = m1.merge(m2, prefer='self')

    # U235/fission: self wins
    assert merged['U235', 'fission'] == pytest.approx([1.])
    # U235/(n,2n): other used
    assert merged['U235', '(n,2n)'] == pytest.approx([5.])
