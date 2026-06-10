"""Unit tests for the continuous-energy group cross section / flux-collapse path:
``openmc.lib.Nuclide.group_xs``, ``_SparseXSTable``, ``_build_xs_table_ce`` and
``_collapse_fluxes`` in :mod:`openmc.deplete.microxs`.

Requires nuclear data (``OPENMC_CROSS_SECTIONS``) and the compiled library.
"""

from pathlib import Path

import numpy as np
import pytest

import openmc
import openmc.lib
from openmc.data import REACTION_MT
from openmc.mgxs import GROUP_STRUCTURES
from openmc.deplete import MicroXS
from openmc.deplete.microxs import _build_xs_table_ce, _collapse_fluxes
from openmc.deplete.coupled_operator import (
    _find_cross_sections, _get_nuclides_with_data)

CHAIN_FILE = Path(__file__).parents[1] / "chain_simple.xml"
GROUP_STRUCTURE = "CASMO-40"
TEMPERATURE = 293.6


def _energies():
    return np.asarray(GROUP_STRUCTURES[GROUP_STRUCTURE], dtype=float)


def _reference_micros(nuclides, reactions, mts, energies, fluxes, temperature):
    """Old per-domain reference: collapse_rate(mt, T, E, flux/sum) per (nuc, mt)."""
    out = []
    with openmc.lib.TemporarySession():
        lib_nucs = {n: openmc.lib.load_nuclide(n) for n in nuclides}
        for flux in fluxes:
            flux = np.asarray(flux, dtype=float)
            arr = np.zeros((len(nuclides), len(reactions), 1))
            flux_sum = flux.sum()
            if flux_sum != 0.0:
                phi = flux / flux_sum
                for i, nuc in enumerate(nuclides):
                    for j, mt in enumerate(mts):
                        arr[i, j, 0] = lib_nucs[nuc].collapse_rate(
                            mt, temperature, energies, phi)
            out.append(arr)
    return out


def test_group_xs_equivalence():
    """group_xs reproduces collapse_rate: sum_g(phi_g * group_xs_g) == collapse_rate.

    A single dot product is under-determined, so equivalence is checked over
    several independent fluxes (random + the flat-per-eV bin-width flux). A scale
    guard catches a forgotten division by the group width (the undivided integral
    would be ~1e5-1e7 eV-b, far outside the barns range).
    """
    energies = _energies()
    n_groups = energies.size - 1
    dE = np.diff(energies)
    rng = np.random.default_rng(0)

    # (nuclide, MT): include a threshold reaction (U238 (n,2n)) plus non-threshold
    cases = [
        ('U238', REACTION_MT['(n,2n)']),
        ('U235', REACTION_MT['fission']),
        ('O16', REACTION_MT['(n,gamma)']),
    ]
    with openmc.lib.TemporarySession():
        for name, mt in cases:
            nuc = openmc.lib.load_nuclide(name)
            group_xs = nuc.group_xs(mt, TEMPERATURE, energies)
            assert group_xs.shape == (n_groups,)

            for flux in (rng.random(n_groups), rng.random(n_groups), dE):
                assert np.dot(flux, group_xs) == pytest.approx(
                    nuc.collapse_rate(mt, TEMPERATURE, energies, flux), rel=1e-10)

            # Per-group AVERAGE cross section, in barns -- not the raw integral.
            assert group_xs.max() < 1e6

            if mt == REACTION_MT['(n,2n)']:
                n2n = group_xs

    # The threshold reaction must leave below-threshold groups at zero.
    assert (n2n == 0.0).any() and (n2n != 0.0).any()


def test_group_xs_structure_below_threshold():
    """A structure entirely below a reaction threshold collapses to all zeros."""
    energies = np.array([1.0e-5, 1.0, 1.0e2])   # U238 (n,2n) threshold ~6.2 MeV
    flux = np.array([1.0, 2.0])
    with openmc.lib.TemporarySession():
        nuc = openmc.lib.load_nuclide('U238')
        mt = REACTION_MT['(n,2n)']
        assert np.all(nuc.group_xs(mt, TEMPERATURE, energies) == 0.0)
        assert nuc.collapse_rate(mt, TEMPERATURE, energies, flux) == 0.0


def test_ce_collapse_matches_collapse_rate():
    """_build_xs_table_ce + _collapse_fluxes (and the delegated from_multigroup_flux)
    reproduce the old per-domain collapse_rate result, with a threshold reaction
    and a zero-flux domain.
    """
    nuclides = ['U235', 'U238', 'O16']
    reactions = ['fission', '(n,gamma)', '(n,2n)']  # (n,2n) is a threshold reaction
    mts = [REACTION_MT[r] for r in reactions]
    energies = _energies()
    n_groups = energies.size - 1
    nuclides_with_data = _get_nuclides_with_data(_find_cross_sections(model=None))

    rng = np.random.default_rng(1)
    fluxes = [rng.random(n_groups), rng.random(n_groups), np.zeros(n_groups)]

    table = _build_xs_table_ce(
        nuclides, reactions, energies, TEMPERATURE, nuclides_with_data)
    new = _collapse_fluxes(table, fluxes, nuclides, reactions)
    ref = _reference_micros(nuclides, reactions, mts, energies, fluxes, TEMPERATURE)

    for micro, ref_data in zip(new, ref):
        assert micro.data == pytest.approx(ref_data, rel=1e-9, abs=1e-12)
    # Zero-flux domain -> all-zero MicroXS of the right shape.
    assert new[-1].data.shape == (len(nuclides), len(reactions), 1)
    assert np.all(new[-1].data == 0.0)

    # from_multigroup_flux delegates to the same engine -> identical result.
    micro = MicroXS.from_multigroup_flux(
        energies=energies, multigroup_flux=fluxes[0], chain_file=CHAIN_FILE,
        nuclides=nuclides, reactions=reactions)
    assert micro.data == pytest.approx(ref[0], rel=1e-9, abs=1e-12)


def test_collapse_fluxes_guards():
    """_collapse_fluxes rejects non-finite / negative flux and zero-fills zero flux."""
    nuclides = ['U235']
    reactions = ['fission']
    energies = _energies()
    n_groups = energies.size - 1
    nuclides_with_data = _get_nuclides_with_data(_find_cross_sections(model=None))
    table = _build_xs_table_ce(
        nuclides, reactions, energies, TEMPERATURE, nuclides_with_data)

    nan_flux = np.ones(n_groups)
    nan_flux[0] = np.nan
    with pytest.raises(ValueError):
        _collapse_fluxes(table, [nan_flux], nuclides, reactions)

    neg_flux = np.ones(n_groups)
    neg_flux[0] = -1.0
    with pytest.raises(ValueError):
        _collapse_fluxes(table, [neg_flux], nuclides, reactions)

    zero = _collapse_fluxes(table, [np.zeros(n_groups)], nuclides, reactions)[0]
    assert zero.data.shape == (1, 1, 1)
    assert np.all(zero.data == 0.0)


def test_from_multigroup_flux_batch():
    """A 2-D batch builds one shared table and matches the per-flux singular
    result and the old per-domain collapse_rate reference, including a threshold
    reaction and a zero-flux row, without mutating the input.
    """
    nuclides = ['U235', 'U238', 'O16']
    reactions = ['fission', '(n,gamma)', '(n,2n)']  # (n,2n) is a threshold reaction
    mts = [REACTION_MT[r] for r in reactions]
    energies = _energies()
    n_groups = energies.size - 1

    rng = np.random.default_rng(2)
    fluxes = np.vstack([rng.random(n_groups), rng.random(n_groups),
                        np.zeros(n_groups)])
    flux_copy = fluxes.copy()

    batch = MicroXS.from_multigroup_flux(
        energies=energies, multigroup_flux=fluxes, chain_file=CHAIN_FILE,
        nuclides=nuclides, reactions=reactions)
    assert isinstance(batch, list)
    assert len(batch) == len(fluxes)
    assert all(isinstance(m, MicroXS) for m in batch)

    # The batch call must not mutate its input
    assert np.array_equal(fluxes, flux_copy)

    # Each batch element equals the per-flux singular call ...
    for row, micro in zip(fluxes, batch):
        single = MicroXS.from_multigroup_flux(
            energies=energies, multigroup_flux=row, chain_file=CHAIN_FILE,
            nuclides=nuclides, reactions=reactions)
        assert isinstance(single, MicroXS)
        assert micro.data == pytest.approx(single.data, rel=1e-9, abs=1e-12)

    # ... and the old per-domain collapse_rate reference
    ref = _reference_micros(nuclides, reactions, mts, energies, fluxes, TEMPERATURE)
    for micro, ref_data in zip(batch, ref):
        assert micro.data == pytest.approx(ref_data, rel=1e-9, abs=1e-12)

    # Zero-flux row -> all-zero MicroXS
    assert np.all(batch[-1].data == 0.0)


def test_from_multigroup_flux_one_row_batch():
    """A 1-row 2-D batch returns a 1-element list, not an unwrapped MicroXS."""
    energies = _energies()
    n_groups = energies.size - 1
    flux = np.random.default_rng(3).random((1, n_groups))

    result = MicroXS.from_multigroup_flux(
        energies=energies, multigroup_flux=flux, chain_file=CHAIN_FILE,
        nuclides=['U235'], reactions=['fission'])
    assert isinstance(result, list)
    assert len(result) == 1
    assert isinstance(result[0], MicroXS)


def test_from_multigroup_flux_invalid_shape():
    """ndim 0/3 and ragged fluxes raise ValueError."""
    energies = _energies()
    n_groups = energies.size - 1
    kwargs = dict(energies=energies, chain_file=CHAIN_FILE,
                  nuclides=['U235'], reactions=['fission'])

    with pytest.raises(ValueError):
        MicroXS.from_multigroup_flux(multigroup_flux=1.0, **kwargs)

    with pytest.raises(ValueError):
        MicroXS.from_multigroup_flux(
            multigroup_flux=np.ones((2, 1, n_groups)), **kwargs)

    with pytest.raises(ValueError):
        MicroXS.from_multigroup_flux(
            multigroup_flux=[np.ones(n_groups), np.ones(n_groups - 1)], **kwargs)


def test_from_multigroup_flux_explicit_chain_no_config(monkeypatch):
    """An explicit chain_file works even when openmc.config has no chain_file."""
    monkeypatch.delitem(openmc.config, 'chain_file', raising=False)
    energies = _energies()
    n_groups = energies.size - 1
    flux = np.random.default_rng(4).random(n_groups)

    # With no nuclides/reactions the chain must be resolved from chain_file
    micro = MicroXS.from_multigroup_flux(
        energies=energies, multigroup_flux=flux, chain_file=CHAIN_FILE)
    assert isinstance(micro, MicroXS)
    assert len(micro.nuclides) > 0
