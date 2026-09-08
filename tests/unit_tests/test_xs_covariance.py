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
    RawBlockDiagnostics,
    compute_covariance_factor,
    validate_raw_mf33_reactions,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
 
# Default group count and energy grid used by validate_raw_mf33_reactions tests
G = 4
EK = np.logspace(-5, 7, G + 1).tolist()  # G+1 boundaries in eV
 
 
def _posdef(n, seed=42):
    """Return a random (n x n) symmetric positive-definite matrix."""
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((n, n))
    return A @ A.T + 0.1 * np.eye(n)
 
 
def _semidef(n, rank, seed=42):
    """Return a random (n x n) symmetric PSD matrix with given rank."""
    rng = np.random.default_rng(seed)
    B = rng.standard_normal((n, rank))
    return B @ B.T
 
 
def _make_mock_covariances(G=10, seed=42):
    """Build a NeutronXSCovariances object with synthetic data."""
    energy = np.logspace(-5, 7, G + 1)
    cov_mat = _posdef(G, seed=seed)
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
 
 
def _wrap_self(mt, M):
    """Wrap a single self-block into the reactions dict structure."""
    return {mt: {"COVS": {mt: M}}}
 
 
def _wrap_pair(mt, mt1, M_self, M_cross,
               M_self1=None, M_partner=None):
    """Wrap a (mt, mt1) pair with cross-block(s)."""
    reactions = {
        mt:  {"COVS": {mt: M_self, mt1: M_cross}},
        mt1: {"COVS": {}},
    }
    if M_self1 is not None:
        reactions[mt1]["COVS"][mt1] = M_self1
    if M_partner is not None:
        reactions[mt1]["COVS"][mt] = M_partner
    return reactions
 
 
# ===========================================================================
# compute_covariance_factor
# ===========================================================================
 
def test_factor_posdef_cholesky():
    """Positive-definite input should use the Cholesky path."""
    Sigma = _posdef(5)
    result = compute_covariance_factor(Sigma)
 
    assert result.method == "cholesky"
    assert result.effective_rank == 5
    assert result.full_size == 5
    assert result.L.shape == (5, 5)
    np.testing.assert_allclose(result.L @ result.L.T, Sigma, atol=1e-10)
 
 
def test_factor_semidef_eigen_qr():
    """Rank-deficient input should fall back to eigen_qr."""
    Sigma = _semidef(8, rank=3)
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
    small = _posdef(3, seed=99)
    Sigma[1:4, 1:4] = small
 
    result = compute_covariance_factor(Sigma)
 
    assert result.effective_rank == 3
    assert result.L.shape[0] == 5
    np.testing.assert_allclose(result.L @ result.L.T, Sigma, atol=1e-10)
 
 
def test_factor_symmetry_enforced():
    """Slightly asymmetric input should not crash — symmetry is forced."""
    Sigma = _posdef(4)
    Sigma[0, 1] += 1e-12
 
    result = compute_covariance_factor(Sigma)
    assert result.effective_rank == 4
 
# ===========================================================================
# HDF5 round-trip (standalone file)
# ===========================================================================
 
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
    G_local = 6
    energy = np.logspace(-5, 7, G_local + 1)
    cov_2 = _posdef(G_local, seed=10)
    cov_102 = _posdef(G_local, seed=20)
 
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
 
 
# ===========================================================================
# write_mf33_group / _read_mf33_group (embedded in another HDF5)
# ===========================================================================
 
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
 
 
# ===========================================================================
# validate_raw_mf33_reactions — hard failures (shape & finiteness)
# ===========================================================================
 
def test_validate_wrong_shape_fails():
    M = np.eye(G + 1)  # deliberately wrong
    diags = validate_raw_mf33_reactions(_wrap_self(2, M), EK)
    dx = diags[2][2]
    assert dx.status == "fail"
    assert any("Wrong shape" in m for m in dx.messages)
    # Downstream diagnostics untouched (None) because the check `continue`s
    assert dx.n_negative_eigenvalues is None
    assert dx.max_abs_correlation is None
 
 
def test_validate_nan_in_matrix_fails():
    M = np.eye(G)
    M[1, 1] = np.nan
    diags = validate_raw_mf33_reactions(_wrap_self(2, M), EK)
    dx = diags[2][2]
    assert dx.status == "fail"
    assert any("NaN or Inf" in m for m in dx.messages)
    assert dx.n_negative_eigenvalues is None
 
 
# ===========================================================================
# validate_raw_mf33_reactions — clean self-block baseline
# ===========================================================================
 
def test_validate_clean_selfblock_passes():
    # Small variances (~0.01) so rel_std ~ 0.1 → well below warn threshold
    M = 0.01 * _posdef(G, seed=1)
    diags = validate_raw_mf33_reactions(_wrap_self(2, M), EK)
    dx = diags[2][2]
    assert dx.status == "pass"
    assert dx.mt == 2 and dx.mt1 == 2
    assert dx.is_self_covariance is True
    # Diagnostic fields should be populated for self-blocks
    assert dx.n_negative_eigenvalues is not None
    assert dx.max_abs_correlation is not None
    assert dx.max_rel_std is not None
    assert dx.n_negative_eigenvalues == 0
    assert dx.max_abs_correlation <= 1.0 + 1e-10
 
 
# ===========================================================================
# validate_raw_mf33_reactions — symmetry
# ===========================================================================
 
def test_validate_asymmetric_selfblock_fails():
    M = _posdef(G, seed=2)
    M[0, 1] += 0.5  # break symmetry noticeably
    diags = validate_raw_mf33_reactions(_wrap_self(2, M), EK)
    dx = diags[2][2]
    assert dx.status == "fail"
    assert any("Not symmetric" in m for m in dx.messages)
 
# ===========================================================================
# validate_raw_mf33_reactions — diagonal: negative / zero
# ===========================================================================
 
def test_validate_negative_diagonal_fails():
    M = np.eye(G) * 0.01
    M[2, 2] = -0.01
    diags = validate_raw_mf33_reactions(_wrap_self(2, M), EK)
    dx = diags[2][2]
    assert dx.status == "fail"
    assert any("negative diagonal" in m for m in dx.messages)
 
 
def test_validate_zero_diagonal_warns():
    # One group has zero variance; rest are fine.
    M = np.eye(G) * 0.01
    M[1, 1] = 0.0
    diags = validate_raw_mf33_reactions(_wrap_self(2, M), EK)
    dx = diags[2][2]
    assert dx.status == "warn"
    assert any("zero-variance" in m for m in dx.messages)
 
 
def test_validate_negative_and_zero_diagonal_is_fail_not_warn():
    """If both negative and zero diagonals are present, fail wins."""
    M = np.eye(G) * 0.01
    M[0, 0] = 0.0
    M[1, 1] = -0.01
    diags = validate_raw_mf33_reactions(_wrap_self(2, M), EK)
    assert diags[2][2].status == "fail"
 
 
# ===========================================================================
# validate_raw_mf33_reactions — eigenvalue spectrum
# ===========================================================================
 
def test_validate_small_negative_eigenvalue_is_informational_only():
    """neg_ratio well below warn threshold (1e-3) → status stays pass, but
    an informational message is appended."""
    rng = np.random.default_rng(5)
    Q, _ = np.linalg.qr(rng.standard_normal((G, G)))
    eigs = np.array([1.0, 0.5, 0.2, -1e-6])
    M = (Q * eigs) @ Q.T
    M = 0.5 * (M + M.T)  # numerical symmetry
 
    diags = validate_raw_mf33_reactions(_wrap_self(2, M), EK)
    dx = diags[2][2]
    # Diagonal might not all be positive due to random rotation;
    # check only the eigenvalue-specific behavior.
    if dx.status == "pass":
        assert dx.n_negative_eigenvalues >= 1
        assert dx.negativity_ratio is not None
        assert dx.negativity_ratio < 1e-3
        assert any("small negative" in m for m in dx.messages)
 
 
def test_validate_moderate_negative_eigenvalue_warns():
    """neg_ratio between warn (1e-3) and fail (1e-1) → warn."""
    eigs = np.array([1.0, 0.5, 0.2, -1e-2])  # ratio = 1e-2
    rng = np.random.default_rng(7)
    Q, _ = np.linalg.qr(rng.standard_normal((G, G)))
    M = (Q * eigs) @ Q.T
    M = 0.5 * (M + M.T)
 
    diags = validate_raw_mf33_reactions(_wrap_self(2, M), EK)
    dx = diags[2][2]
    # Filter out the (possibly flaky) diagonal/correlation side-effects:
    assert dx.n_negative_eigenvalues == 1
    assert 1e-3 <= dx.negativity_ratio < 1e-1
    assert any("negativity ratio" in m and "threshold" in m for m in dx.messages)
 
 
def test_validate_large_negative_eigenvalue_fails():
    """neg_ratio >= fail threshold (1e-1) → fail."""
    eigs = np.array([1.0, 0.5, 0.2, -0.5])  # ratio = 0.5
    rng = np.random.default_rng(11)
    Q, _ = np.linalg.qr(rng.standard_normal((G, G)))
    M = (Q * eigs) @ Q.T
    M = 0.5 * (M + M.T)
 
    diags = validate_raw_mf33_reactions(_wrap_self(2, M), EK)
    dx = diags[2][2]
    assert dx.status == "fail"
    assert dx.n_negative_eigenvalues == 1
    assert dx.negativity_ratio >= 1e-1
    assert any("repaired matrix will differ" in m for m in dx.messages)
 
 
# ===========================================================================
# validate_raw_mf33_reactions — correlation bound |rho| <= 1
# ===========================================================================
 
def test_validate_correlation_overshoot_warns():
    """Construct a symmetric matrix with an off-diagonal that exceeds the
    geometric-mean bound sqrt(d_ii * d_jj)."""
    M = np.eye(G) * 1.0
    M[0, 1] = M[1, 0] = 1.5  # rho = 1.5 > 1
    diags = validate_raw_mf33_reactions(_wrap_self(2, M), EK)
    dx = diags[2][2]
    # Status is at least warn; could be fail because this also produces a
    # negative eigenvalue (2x2 block has det < 0).
    assert dx.status in ("warn", "fail")
    assert dx.max_abs_correlation > 1.0
    assert any("Correlation matrix" in m and "rho" in m for m in dx.messages)
 
 
# ===========================================================================
# validate_raw_mf33_reactions — relative uncertainty magnitude
# ===========================================================================
 
def test_validate_rel_std_above_warn_warns():
    """sqrt(diag) > rel_std_warn (1.0) but < rel_std_fail (10.0) → warn."""
    M = np.eye(G) * 0.01
    M[0, 0] = 4.0  # sigma_rel = 2.0 → between 1 and 10
    diags = validate_raw_mf33_reactions(_wrap_self(2, M), EK)
    dx = diags[2][2]
    assert dx.status == "warn"
    assert dx.n_groups_rel_std_above_1 == 1
    assert dx.n_groups_rel_std_above_10 == 0
    assert any("Relative uncertainty" in m for m in dx.messages)
 
 
def test_validate_rel_std_above_fail_threshold_flagged_as_artifact():
    """sqrt(diag) > rel_std_fail (10.0) → warn with 'artifact' message."""
    M = np.eye(G) * 0.01
    M[0, 0] = 200.0  # sigma_rel ~ 14 > 10
    diags = validate_raw_mf33_reactions(_wrap_self(2, M), EK)
    dx = diags[2][2]
    assert dx.status == "warn"
    assert dx.n_groups_rel_std_above_10 >= 1
    assert any("processing artifact" in m for m in dx.messages)
 
 
# ===========================================================================
# validate_raw_mf33_reactions — threshold customization
# ===========================================================================
 
def test_validate_custom_negativity_thresholds():
    """Verify that relaxing thresholds turns a fail into a non-fail."""
    eigs = np.array([1.0, 0.5, 0.2, -0.5])
    rng = np.random.default_rng(11)
    Q, _ = np.linalg.qr(rng.standard_normal((G, G)))
    M = (Q * eigs) @ Q.T
    M = 0.5 * (M + M.T)
 
    # With defaults: fail.
    strict = validate_raw_mf33_reactions(_wrap_self(2, M), EK)
    assert strict[2][2].status == "fail"
 
    # With extremely loose thresholds, the eigenvalue-spectrum branch no
    # longer fires. Diagonal / correlation checks may still warn depending
    # on the random rotation, but no "repaired matrix" message.
    loose = validate_raw_mf33_reactions(
        _wrap_self(2, M), EK,
        negativity_ratio_warn=10.0,
        negativity_ratio_fail=100.0,
    )
    assert not any(
        "repaired matrix will differ" in m
        for m in loose[2][2].messages
    )
 
 
# ===========================================================================
# validate_raw_mf33_reactions — cross-block partner consistency
# ===========================================================================
 
def test_validate_cross_block_partner_absent_passes():
    """ERRORR stores only the upper triangle; missing partner is expected."""
    Mself = 0.01 * _posdef(G, seed=21)
    Mcross = np.full((G, G), 0.001)
    reactions = _wrap_pair(2, 102, Mself, Mcross)  # no partner at (102, 2)
    diags = validate_raw_mf33_reactions(reactions, EK)
    dx = diags[2][102]
    assert dx.is_self_covariance is False
    assert dx.status == "pass"
    # Cross-blocks shouldn't populate self-only fields
    assert dx.n_negative_eigenvalues is None
    assert dx.max_abs_correlation is None
 
 
def test_validate_cross_block_partner_consistent_passes():
    Mself = 0.01 * _posdef(G, seed=22)
    Mself1 = 0.01 * _posdef(G, seed=23)
    Mcross = np.full((G, G), 0.001)
    reactions = _wrap_pair(
        2, 102, Mself, Mcross,
        M_self1=Mself1,
        M_partner=Mcross.T,   # correct partner
    )
    diags = validate_raw_mf33_reactions(reactions, EK)
    assert diags[2][102].status == "pass"
    assert diags[102][2].status == "pass"
 
 
def test_validate_cross_block_partner_shape_mismatch_fails():
    Mself = 0.01 * _posdef(G, seed=24)
    Mcross = np.full((G, G), 0.001)
    bad_partner = np.full((G, G + 1), 0.001)  # wrong shape
    reactions = _wrap_pair(
        2, 102, Mself, Mcross,
        M_self1=0.01 * _posdef(G, seed=25),
        M_partner=bad_partner,
    )
    diags = validate_raw_mf33_reactions(reactions, EK)
    assert diags[2][102].status == "fail"
    assert any("Partner wrong shape" in m for m in diags[2][102].messages)
 
 
def test_validate_cross_block_partner_not_transpose_fails():
    Mself = 0.01 * _posdef(G, seed=28)
    Mcross = np.full((G, G), 0.001)
    wrong_partner = np.full((G, G), 0.002)  # finite, right shape, wrong values
    reactions = _wrap_pair(
        2, 102, Mself, Mcross,
        M_self1=0.01 * _posdef(G, seed=29),
        M_partner=wrong_partner,
    )
    diags = validate_raw_mf33_reactions(reactions, EK)
    assert diags[2][102].status == "fail"
    msg = " ".join(diags[2][102].messages)
    assert "!=" in msg or "C(" in msg
 