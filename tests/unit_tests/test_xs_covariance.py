"""Unit tests for openmc.data.xs_covariance_njoy

These tests use synthetic data only (no NJOY executable or ENDF files
are required).  They test compute_covariance_factor on positive-definite, semi-definite, 
and zero matrices. The HDF5 round-trip for NeutronXSCovariances (with and without raw 
covariance storage). Finally, its testing the integration with export_to_hdf5
 and from_hdf5 via the mg_covariance property.
"""

from pathlib import Path

import h5py
import numpy as np
import pytest

from openmc.data.xs_covariance_njoy import (
    CovFactorResult,
    NeutronXSCovariances,
    compute_covariance_factor,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_posdef_matrix(n, seed=42):
    """Return a random (n x n) symmetric positive-definite matrix."""
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((n, n))
    return A @ A.T + 0.1 * np.eye(n)


def _make_semidef_matrix(n, rank, seed=42):
    """Return a random (n x n) symmetric PSD matrix with given rank."""
    rng = np.random.default_rng(seed)
    B = rng.standard_normal((n, rank))
    return B @ B.T


def _make_mock_covariances(G=10, seed=42):
    """Build a NeutronXSCovariances object with synthetic data."""
    energy = np.logspace(-5, 7, G + 1)
    cov_mat = _make_posdef_matrix(G, seed=seed)
    result = compute_covariance_factor(cov_mat)

    reactions = {
        2: {
            "MAT": 2631, "MF": 33, "MT": 2,
            "ZA": 26056.0, "AWR": 55.45,
            "COVS": {2: cov_mat},
        },
    }
    factors = {2: {2: result}}

    return NeutronXSCovariances(
        name="Fe56",
        energy_grid_ev=energy,
        reactions=reactions,
        mat=2631,
        temperature_k=293.6,
        factor_results=factors,
    )


# ---------------------------------------------------------------------------
# Tests — compute_covariance_factor
# ---------------------------------------------------------------------------

def test_factor_posdef_cholesky():
    """Positive-definite input should use the Cholesky path."""
    Sigma = _make_posdef_matrix(5)
    result = compute_covariance_factor(Sigma)

    assert result.method == "cholesky"
    assert result.effective_rank == 5
    assert result.full_size == 5
    assert result.L.shape == (5, 5)
    np.testing.assert_allclose(result.L @ result.L.T, Sigma, atol=1e-10)


def test_factor_semidef_eigen_qr():
    """Rank-deficient input should fall back to eigen_qr."""
    Sigma = _make_semidef_matrix(8, rank=3)
    result = compute_covariance_factor(Sigma)

    assert result.method == "eigen_qr"
    assert result.effective_rank == 3
    assert result.full_size == 8
    assert result.L.shape == (8, 3)
    np.testing.assert_allclose(result.L @ result.L.T, Sigma, atol=1e-8)


def test_factor_zero_matrix():
    """All-zero input should return rank-0 with the zero_matrix method."""
    result = compute_covariance_factor(np.zeros((6, 6)))

    assert result.method == "zero_matrix"
    assert result.effective_rank == 0
    assert result.full_size == 6
    assert result.L.shape == (6, 0)


def test_factor_partial_zero_variance():
    """Matrix with some zero-diagonal rows should still factorise the
    non-zero sub-block correctly."""
    Sigma = np.zeros((5, 5))
    small = _make_posdef_matrix(3, seed=99)
    Sigma[1:4, 1:4] = small

    result = compute_covariance_factor(Sigma)

    assert result.effective_rank == 3
    assert result.L.shape[0] == 5
    np.testing.assert_allclose(result.L @ result.L.T, Sigma, atol=1e-10)


def test_factor_symmetry_enforced():
    """Slightly asymmetric input should not crash — symmetry is forced."""
    Sigma = _make_posdef_matrix(4)
    Sigma[0, 1] += 1e-12

    result = compute_covariance_factor(Sigma)
    assert result.effective_rank == 4


def test_cov_factor_result_fields():
    """CovFactorResult should have exactly the four documented fields."""
    r = CovFactorResult(
        L=np.eye(2), effective_rank=2, full_size=2, method="test",
    )
    assert r.L.shape == (2, 2)
    assert r.effective_rank == 2
    assert r.full_size == 2
    assert r.method == "test"


# ---------------------------------------------------------------------------
# Tests — HDF5 round-trip (standalone file)
# ---------------------------------------------------------------------------

def test_hdf5_roundtrip_with_raw_covariance(tmp_path):
    """Full round-trip storing both raw covariance and L factors."""
    cov = _make_mock_covariances(G=10)
    h5_path = tmp_path / "cov_test.h5"

    cov.to_hdf5(h5_path, store_raw_covariance=True)
    cov_read = NeutronXSCovariances.from_hdf5(h5_path)

    # Metadata
    assert cov_read.mat == cov.mat
    assert cov_read.temperature_k == cov.temperature_k
    np.testing.assert_allclose(cov_read.energy_grid_ev, cov.energy_grid_ev)

    # Reaction keys preserved
    assert set(cov_read.reactions.keys()) == set(cov.reactions.keys())

    # Raw covariance matrices match
    for mt in cov.reactions:
        for mt1 in cov.reactions[mt]["COVS"]:
            original = cov.reactions[mt]["COVS"][mt1]
            loaded = cov_read.reactions[mt]["COVS"][mt1]
            np.testing.assert_allclose(loaded, original, atol=1e-12)

    # L factors match
    for mt in cov.factor_results:
        for mt1 in cov.factor_results[mt]:
            L_orig = cov.factor_results[mt][mt1].L
            L_read = cov_read.factor_results[mt][mt1].L
            np.testing.assert_allclose(L_read, L_orig, atol=1e-12)


def test_hdf5_roundtrip_factors_only(tmp_path):
    """Round-trip with store_raw_covariance=False.

    The raw covariance should be reconstructed as L @ L.T on read.
    """
    cov = _make_mock_covariances(G=8)
    h5_path = tmp_path / "cov_factors_only.h5"

    cov.to_hdf5(h5_path, store_raw_covariance=False)
    cov_read = NeutronXSCovariances.from_hdf5(h5_path)

    # Verify "reactions" group is absent in the HDF5
    with h5py.File(h5_path, "r") as f:
        assert "reactions" not in f["mf33"]

    # COVS should be reconstructed from L factors
    for mt in cov.reactions:
        for mt1 in cov.reactions[mt]["COVS"]:
            original = cov.reactions[mt]["COVS"][mt1]
            reconstructed = cov_read.reactions[mt]["COVS"][mt1]
            np.testing.assert_allclose(reconstructed, original, atol=1e-8)


def test_hdf5_schema_attributes(tmp_path):
    """Verify the expected HDF5 attributes and datasets exist."""
    cov = _make_mock_covariances(G=5)
    h5_path = tmp_path / "schema_check.h5"
    cov.to_hdf5(h5_path)

    with h5py.File(h5_path, "r") as f:
        mf33 = f["mf33"]
        assert mf33.attrs["format"] == b"openmc.mf33.v1"
        assert mf33.attrs["relative"] == 1
        assert mf33.attrs["mat"] == 2631
        assert "energy_grid_ev" in mf33
        assert "mts" in mf33
        assert "factors" in mf33

        # Energy grid has G+1 entries
        assert mf33["energy_grid_ev"].shape == (6,)


def test_hdf5_multiple_reactions(tmp_path):
    """Round-trip with multiple MT values."""
    G = 6
    energy = np.logspace(-5, 7, G + 1)
    cov_2 = _make_posdef_matrix(G, seed=10)
    cov_102 = _make_posdef_matrix(G, seed=20)

    reactions = {
        2: {"MAT": 2631, "MF": 33, "MT": 2,
            "ZA": 26056.0, "AWR": 55.45,
            "COVS": {2: cov_2}},
        102: {"MAT": 2631, "MF": 33, "MT": 102,
              "ZA": 26056.0, "AWR": 55.45,
              "COVS": {102: cov_102}},
    }
    cov_obj = NeutronXSCovariances(
        name="Fe56", energy_grid_ev=energy,
        reactions=reactions, mat=2631, temperature_k=293.6,
    )
    h5_path = tmp_path / "multi_mt.h5"
    cov_obj.to_hdf5(h5_path)
    cov_read = NeutronXSCovariances.from_hdf5(h5_path)

    assert set(cov_read.reactions.keys()) == {2, 102}
    np.testing.assert_allclose(
        cov_read.reactions[2]["COVS"][2], cov_2, atol=1e-12,
    )
    np.testing.assert_allclose(
        cov_read.reactions[102]["COVS"][102], cov_102, atol=1e-12,
    )


# ---------------------------------------------------------------------------
# Tests — write_mf33_group / _read_mf33_group (embedded in another HDF5)
# ---------------------------------------------------------------------------

def test_write_read_mf33_into_existing_group(tmp_path):
    """Write mf33 into a pre-existing HDF5 group, then read it back."""
    cov = _make_mock_covariances(G=5)
    h5_path = tmp_path / "embedded.h5"

    with h5py.File(h5_path, "w") as f:
        nuc = f.create_group("Fe56")
        cov_root = nuc.create_group("covariance")
        cov.write_mf33_group(cov_root, store_raw_covariance=True)

    with h5py.File(h5_path, "r") as f:
        mf33_group = f["Fe56"]["covariance"]["mf33"]
        cov_read = NeutronXSCovariances._read_mf33_group(
            mf33_group, name="Fe56",
        )

    assert cov_read.name == "Fe56"
    assert cov_read.mat == 2631
    assert 2 in cov_read.reactions


def test_overwrite_existing_mf33(tmp_path):
    """Calling write_mf33_group twice should replace, not duplicate."""
    cov = _make_mock_covariances(G=4)
    h5_path = tmp_path / "overwrite.h5"

    with h5py.File(h5_path, "w") as f:
        root = f.create_group("root")
        cov.write_mf33_group(root)
        cov.write_mf33_group(root)

    with h5py.File(h5_path, "r") as f:
        assert "mf33" in f["root"]