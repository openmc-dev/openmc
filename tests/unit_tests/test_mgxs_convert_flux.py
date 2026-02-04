"""Tests for openmc.mgxs.convert_flux_groups function."""

import numpy as np
import pytest
from pytest import approx

import openmc.mgxs


def test_coarse_to_fine():
    """Test coarse to fine conversion with flux conservation."""
    source = openmc.mgxs.EnergyGroups([1.0, 10.0, 100.0])
    target = openmc.mgxs.EnergyGroups([1.0, 5.0, 10.0, 50.0, 100.0])
    flux_source = np.array([1e8, 2e8])

    flux_target = openmc.mgxs.convert_flux_groups(flux_source, source, target)

    # Check conservation
    assert np.sum(flux_target) == approx(np.sum(flux_source))
    assert len(flux_target) == 4
    assert np.all(flux_target >= 0)


def test_fine_to_coarse():
    """Test fine to coarse conversion (reverse direction)."""
    source = openmc.mgxs.EnergyGroups([1.0, 5.0, 10.0, 50.0, 100.0])
    target = openmc.mgxs.EnergyGroups([1.0, 10.0, 100.0])
    flux_source = np.array([1e7, 2e7, 3e7, 4e7])

    flux_target = openmc.mgxs.convert_flux_groups(flux_source, source, target)

    assert np.sum(flux_target) == approx(np.sum(flux_source))
    assert len(flux_target) == 2


def test_lethargy_distribution():
    """Test that flux is distributed by lethargy, not linear energy."""
    # Single group from 1 to 100 eV
    source = openmc.mgxs.EnergyGroups([1.0, 100.0])
    # Split into two groups: 1-10 eV and 10-100 eV
    target = openmc.mgxs.EnergyGroups([1.0, 10.0, 100.0])
    flux_source = np.array([1e8])

    flux_target = openmc.mgxs.convert_flux_groups(flux_source, source, target)

    # Each target group spans one decade (ln(10) lethargy each)
    # So flux should be split 50/50 by lethargy
    assert flux_target[0] == approx(5e7)
    assert flux_target[1] == approx(5e7)


def test_fns_ccfe709_to_ukaea1102():
    """Test CCFE-709 to UKAEA-1102 conversion with real FNS flux spectrum."""
    from pathlib import Path
    flux_file = Path(__file__).parent.parent / 'fns_flux_709.npy'
    fns_flux_709 = np.load(flux_file)

    flux_1102 = openmc.mgxs.convert_flux_groups(fns_flux_709, 'CCFE-709', 'UKAEA-1102')

    assert len(flux_1102) == 1102
    assert np.sum(flux_1102) == approx(np.sum(fns_flux_709), rel=1e-10)
    assert np.all(flux_1102 >= 0)
