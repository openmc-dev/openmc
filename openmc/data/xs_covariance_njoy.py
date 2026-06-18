"""openmc.data.xs_covariance

Tiny model + HDF5 I/O for multigroup **cross section** covariances (ERRORR MF=33).

This is intentionally narrow:
- neutron MF=33 only
- explicit energy grid (group edges in eV)
- **relative** covariances (what ERRORR typically produces with IRELCO=1)

Public API
----------
NeutronXSCovariances.from_endf(...)
NeutronXSCovariances.to_hdf5(...)
NeutronXSCovariances.from_hdf5(...)

Everything else is internal parsing.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import scipy.linalg as la
from .mf33_njoy import generate_errorr_mf33

log = logging.getLogger(__name__)

# -----------------------------------------------------------------------------
# ERRORR tape33 (ENDF-6 text) parser for MF=33
# -----------------------------------------------------------------------------

from collections import namedtuple

_CONT = namedtuple("CONT", "C1 C2 L1 L2 N1 N2")
_LIST = namedtuple("LIST", "C1 C2 L1 L2 NPL N2 B")


def _endf_int(field: str) -> int:
    field = field.strip()
    if not field:
        return 0
    try:
        return int(field)
    except ValueError:
        return int(float(field))


def _endf_float(field: str) -> float:
    s = field.strip()
    if not s:
        return 0.0
    if "e" in s.lower():
        return float(s)
    # ENDF-style: mantissa + exponent sign+digits, no 'e'
    for i in range(len(s) - 1, 0, -1):
        if s[i] in "+-" and s[i - 1].isdigit():
            mant = s[:i]
            exp = s[i:]
            return float(f"{mant}e{exp}")
    return float(s)


def _parse_endf_ids(line: str) -> Tuple[int, int, int]:
    if len(line) < 75:
        return (0, 0, 0)
    return (int(line[66:70]), int(line[70:72]), int(line[72:75]))


def _read_cont_line(line: str) -> _CONT:
    return _CONT(
        _endf_float(line[0:11]),
        _endf_float(line[11:22]),
        _endf_int(line[22:33]),
        _endf_int(line[33:44]),
        _endf_int(line[44:55]),
        _endf_int(line[55:66]),
    )


def _read_list_record(lines: List[str], ipos: int) -> Tuple[_LIST, int]:
    C = _read_cont_line(lines[ipos])
    npl = int(C.N1)
    ipos += 1
    vals: List[float] = []
    nlines = int(np.ceil(npl / 6.0))
    for _ in range(nlines):
        ln = lines[ipos]
        for j in range(6):
            a = 11 * j
            b = a + 11
            if a >= 66:
                break
            vals.append(_endf_float(ln[a:b]))
        ipos += 1
    vals = vals[:npl]
    return _LIST(C.C1, C.C2, C.L1, C.L2, npl, C.N2, vals), ipos


def _parse_mf33_section_lines(section_lines: List[str], mat: int, mt: int) -> Dict[str, Any]:
    i = 0
    C0 = _read_cont_line(section_lines[i])
    i += 1
    out: Dict[str, Any] = {"MAT": mat, "MF": 33, "MT": mt, "ZA": C0.C1, "AWR": C0.C2}

    reaction_pairs: Dict[int, np.ndarray] = {}
    for _ in range(int(C0.N2)):  # number of reaction pairs
        C = _read_cont_line(section_lines[i])
        i += 1
        mt1 = int(C.L2)
        ng = int(C.N2)
        M = np.zeros((ng, ng))
        while True:
            L, i = _read_list_record(section_lines, i)
            ngcol = int(L.L1)
            grow = int(L.N2)
            gcol = int(L.L2)
            # LIST rows are 1-based
            M[grow - 1, gcol - 1 : gcol - 1 + ngcol] = np.asarray(L.B, dtype=float)
            if grow >= ng:
                break
        reaction_pairs[mt1] = M

    out["COVS"] = reaction_pairs
    return out


def parse_errorr_mf33_text(tape33_text: str, mat: Optional[int] = None) -> Dict[int, Dict[str, Any]]:
    lines = tape33_text.splitlines(True)
    reactions: Dict[int, Dict[str, Any]] = {}

    current: List[str] = []
    current_key: Optional[Tuple[int, int, int]] = None

    def flush():
        nonlocal current, current_key
        if current_key is None or not current:
            current = []
            current_key = None
            return
        m, mf, mt = current_key
        if mf == 33 and mt > 0 and (mat is None or m == mat):
            reactions[mt] = _parse_mf33_section_lines(current, m, mt)
        current = []
        current_key = None

    for ln in lines:
        m, mf, mt = _parse_endf_ids(ln)

        if m == 0 and mf == 0 and mt == 0:
            flush()
            continue
        if mt == 0:
            flush()
            continue

        key = (m, mf, mt)
        if current_key is None:
            current_key = key
        elif key != current_key:
            flush()
            current_key = key

        current.append(ln)

    flush()
    return reactions

# -----------------------------------------------------------------------------
# Validating raw covariance data from NJOY
# -----------------------------------------------------------------------------
@dataclass
class RawBlockDiagnostics:
    """Stage-1 QA diagnostics for one raw MF=33 block."""
    mt: int
    mt1: int
    is_self_covariance: bool
    status: str = "pass"   # "pass", "warn", "fail"
    messages: List[str] = field(default_factory=list)

    # --- eigenvalue spectrum diagnostics (self-blocks only) ---
    n_negative_eigenvalues: Optional[int] = None
    negativity_ratio: Optional[float] = None      # |λ_min_neg| / λ_max_pos
    negative_spectral_fraction: Optional[float] = None  # Σ|λ_neg| / Σ|λ_all|

    # --- correlation matrix bounds (self-blocks only) ---
    max_abs_correlation: Optional[float] = None

    # --- relative uncertainty magnitude (self-blocks only) ---
    max_rel_std: Optional[float] = None
    n_groups_rel_std_above_1: Optional[int] = None
    n_groups_rel_std_above_10: Optional[int] = None
 
    def _flag(self, level: str, msg: str):
        """Set status to *level* (unless already 'fail') and append *msg*."""
        if level == "fail" or self.status != "fail":
            self.status = level
        self.messages.append(msg)
 
 
def _rel_frob(A: np.ndarray, B: np.ndarray) -> float:
    """||A-B||_F / max(||A||_F, ||B||_F, 1)."""
    nA, nB = np.linalg.norm(A, "fro"), np.linalg.norm(B, "fro")
    return float(np.linalg.norm(A - B, "fro") / max(nA, nB, 1.0))
 
 
def validate_raw_mf33_reactions(
    reactions: Dict[int, Dict[str, Any]],
    energy_grid_ev: Sequence[float],
    *,
    atol: float = 1e-14,
    rtol: float = 1e-10,
    warn_only_missing_partner: bool = True,
    negativity_ratio_warn: float = 1e-3,
    negativity_ratio_fail: float = 1e-1,
    rel_std_warn: float = 1.0,
    rel_std_fail: float = 10.0,
) -> Dict[int, Dict[int, RawBlockDiagnostics]]:
    """Stage-1 integrity checks on raw MF=33 matrices from ERRORR tape33.
 
    Checks per block: shape == (G,G), finiteness, and—for self blocks—
    symmetry, non-negative diagonal, eigenvalue spectrum, correlation
    matrix bounds, and relative uncertainty magnitude; for cross blocks,
    transpose-partner consistency C(mt,mt1) ≈ C(mt1,mt)^T.

    Parameters
    ----------
    negativity_ratio_warn, negativity_ratio_fail
        Thresholds on |λ_min_neg| / λ_max_pos for the eigenvalue spectrum
        check.  Following the convention that 1e-3 is a warning and 1e-1
        signals the approximate (repaired) matrix will differ substantially.
    rel_std_warn, rel_std_fail
        Thresholds on group-wise relative standard deviations
        (sqrt of diagonal, since IRELCO=1).  Values above rel_std_warn
        (default 1.0 = 100%) are flagged as warnings; above rel_std_fail
        (default 10.0 = 1000%) as likely processing artifacts.
 
    Returns nested dict ``diagnostics[mt][mt1] = RawBlockDiagnostics``.
    """
    G = len(energy_grid_ev) - 1
    shape_exp = (G, G)
    diags: Dict[int, Dict[int, RawBlockDiagnostics]] = {}
 
    for mt, sec in reactions.items():
        diags[mt] = {}
        for mt1, M in sec.get("COVS", {}).items():
            A = np.asarray(M, dtype=np.float64)
            is_self = int(mt) == int(mt1)
            dx = RawBlockDiagnostics(int(mt), int(mt1), is_self)
 
            # 1) Shape / finiteness (hard failures) -------------------------
            if A.shape != shape_exp:
                dx._flag("fail", f"Wrong shape {A.shape}, expected {shape_exp}.")
            if not np.all(np.isfinite(A)):
                dx._flag("fail", "Matrix contains NaN or Inf values.")
            if dx.status == "fail":
                diags[mt][mt1] = dx
                continue
 
            # 2) Self-covariance: symmetry + diagonal -----------------------
            if is_self:
                res = _rel_frob(A, A.T)
                if not np.allclose(A, A.T, atol=atol, rtol=rtol):
                    dx._flag("fail", f"Not symmetric (residual={res:.3e}).")
                d = np.diag(A)
                n_neg = int(np.count_nonzero(d < -atol))
                n_zero = int(np.count_nonzero(np.isclose(d, 0.0, atol=atol)))
                if n_neg > 0:
                    dx._flag("fail", f"{n_neg} negative diagonal entries; "
                             f"min(diag)={float(d.min()):.6e}.")
                elif n_zero > 0:
                    dx._flag("warn", f"{n_zero} zero-variance diagonal groups.")

                # 2a) Eigenvalue spectrum diagnostics -----------------------
                #     Symmetrise before computing eigenvalues so that eigh
                #     is applicable even when ERRORR left tiny asymmetries.
                A_sym = 0.5 * (A + A.T)
                eigenvalues = la.eigvalsh(A_sym)

                neg_mask = eigenvalues < -atol
                pos_mask = eigenvalues > atol
                n_neg_eig = int(neg_mask.sum())
                dx.n_negative_eigenvalues = n_neg_eig

                lam_max_pos = float(eigenvalues[pos_mask].max()) if pos_mask.any() else 0.0
                abs_neg = np.abs(eigenvalues[neg_mask]) if neg_mask.any() else np.empty(0)
                lam_min_neg_abs = float(abs_neg.max()) if abs_neg.size > 0 else 0.0

                if lam_max_pos > 0.0 and lam_min_neg_abs > 0.0:
                    dx.negativity_ratio = lam_min_neg_abs / lam_max_pos
                else:
                    dx.negativity_ratio = 0.0

                total_spectral_weight = float(np.abs(eigenvalues).sum())
                if total_spectral_weight > 0.0:
                    dx.negative_spectral_fraction = float(abs_neg.sum()) / total_spectral_weight
                else:
                    dx.negative_spectral_fraction = 0.0

                if n_neg_eig > 0:
                    if dx.negativity_ratio >= negativity_ratio_fail:
                        dx._flag(
                            "fail",
                            f"Eigenvalue spectrum: {n_neg_eig} negative eigenvalue(s), "
                            f"negativity ratio={dx.negativity_ratio:.3e} "
                            f"(>= {negativity_ratio_fail:.0e} threshold); "
                            f"repaired matrix will differ substantially from original."
                        )
                    elif dx.negativity_ratio >= negativity_ratio_warn:
                        dx._flag(
                            "warn",
                            f"Eigenvalue spectrum: {n_neg_eig} negative eigenvalue(s), "
                            f"negativity ratio={dx.negativity_ratio:.3e} "
                            f"(>= {negativity_ratio_warn:.0e} threshold)."
                        )
                    else:
                        # Small negative eigenvalues — informational only
                        dx.messages.append(
                            f"Eigenvalue spectrum: {n_neg_eig} small negative "
                            f"eigenvalue(s), negativity ratio={dx.negativity_ratio:.3e}."
                        )

                # 2b) Correlation matrix bounds check -----------------------
                #     Convert self-covariance to correlation and verify |ρ|≤1.
                d_safe = d.copy()
                d_safe[d_safe <= 0.0] = np.inf  # skip zero-variance groups
                inv_std = 1.0 / np.sqrt(d_safe)
                corr = A_sym * np.outer(inv_std, inv_std)
                # Zero out rows/cols that had zero variance (they are undefined)
                zero_var_mask = d <= 0.0
                corr[zero_var_mask, :] = 0.0
                corr[:, zero_var_mask] = 0.0
                np.fill_diagonal(corr, 1.0)  # diagonal is ρ=1 by definition

                max_abs_rho = float(np.max(np.abs(
                    corr[np.triu_indices_from(corr, k=1)]
                ))) if G > 1 else 0.0
                dx.max_abs_correlation = max_abs_rho

                if max_abs_rho > 1.0 + rtol:
                    n_violating = int(np.count_nonzero(
                        np.abs(corr[np.triu_indices_from(corr, k=1)]) > 1.0 + rtol
                    ))
                    dx._flag(
                        "warn",
                        f"Correlation matrix has {n_violating} off-diagonal "
                        f"element(s) with |rho| > 1 (max |rho|={max_abs_rho:.6f}); "
                        f"unphysical covariance entries."
                    )

                # 2c) Relative uncertainty magnitude check ------------------
                #     For IRELCO=1, sqrt(diag) gives group-wise relative sigma.
                rel_std = np.sqrt(np.maximum(d, 0.0))
                dx.max_rel_std = float(rel_std.max()) if rel_std.size > 0 else 0.0
                dx.n_groups_rel_std_above_1 = int(np.count_nonzero(rel_std > rel_std_warn))
                dx.n_groups_rel_std_above_10 = int(np.count_nonzero(rel_std > rel_std_fail))

                if dx.n_groups_rel_std_above_10 > 0:
                    dx._flag(
                        "warn",
                        f"Relative uncertainty: {dx.n_groups_rel_std_above_10} group(s) "
                        f"with sigma_rel > {rel_std_fail} "
                        f"(max={dx.max_rel_std:.3f}); "
                        f"likely ERRORR processing artifact."
                    )
                elif dx.n_groups_rel_std_above_1 > 0:
                    dx._flag(
                        "warn",
                        f"Relative uncertainty: {dx.n_groups_rel_std_above_1} group(s) "
                        f"with sigma_rel > {rel_std_warn} "
                        f"(max={dx.max_rel_std:.3f})."
                    )
 
            # 3) Cross-covariance: transpose partner ------------------------
            else:
                B_raw = reactions.get(int(mt1), {}).get("COVS", {}).get(int(mt))
                if B_raw is None:
                    pass  # Expected: ERRORR stores only the upper-triangle block
                else:
                    B = np.asarray(B_raw, dtype=np.float64)
                    if B.shape != shape_exp:
                        dx._flag("fail", f"Partner wrong shape {B.shape}.")
                    elif not np.all(np.isfinite(B)):
                        dx._flag("fail", "Partner contains NaN/Inf.")
                    elif not np.allclose(A, B.T, atol=atol, rtol=rtol):
                        dx._flag("fail",
                                f"C({mt},{mt1}) != C({mt1},{mt})^T "
                                f"(residual={_rel_frob(A, B.T):.3e}).")
 
            diags[mt][mt1] = dx
 
    return diags

# -----------------------------------------------------------------------------
# Covariance factor computation (eigendecomposition + QR)
# -----------------------------------------------------------------------------

@dataclass
class CovFactorResult:
    """Container for a covariance factorization after diagnostics are done.

    Attributes
    ----------
    L : np.ndarray
        Lower-triangular factor, shape (G, r) where r = effective_rank.
        Satisfies L @ L.T ≈ original covariance (after any regularization).
    effective_rank : int
        Number of positive eigenvalues retained.
    full_size : int
        Original matrix dimension G.
    method : str
        Which path was taken: "cholesky", "eigen_qr", or "zero_matrix".
    """
    L: np.ndarray
    effective_rank: int
    full_size: int
    method: str = "eigen_qr"
  

def _enforce_symmetry(A: np.ndarray) -> np.ndarray:
    """Force exact symmetry: A_sym = (A + A^T) / 2."""
    return 0.5 * (A + A.T)


def _strip_zero_variance(A: np.ndarray):
    """
    Remove rows/columns with zero diagonal (zero variance).
    """
    diag = np.diag(A)
    nz = np.flatnonzero(diag > 0.0)
    if len(nz) == 0:
        return np.empty((0, 0), dtype=A.dtype), nz
    return A[np.ix_(nz, nz)], nz

def compute_covariance_factor(
    cov_matrix: np.ndarray,
    tol: float = 1e-10,

) -> CovFactorResult:
    """
    Compute a lower-triangular factor L such that L @ L.T ≈ Sigma.
    """
    G = cov_matrix.shape[0]

    # --- 0. Enforce exact symmetry ---
    A = _enforce_symmetry(np.asarray(cov_matrix, dtype=np.float64))

    # --- 1. Strip zero-variance rows/columns ---
    Ar, nz = _strip_zero_variance(A)

    if Ar.size == 0:
        return CovFactorResult(
            L=np.zeros((G, 0), dtype=np.float64),
            effective_rank=0,
            full_size=G,
            method="zero_matrix",
        )

    Gr = Ar.shape[0]

    # --- 2a. Fast path: standard Cholesky ---
    try:
        Lr = la.cholesky(Ar, lower=True)

        L_full = np.zeros((G, Gr), dtype=np.float64)
        L_full[nz, :] = Lr

        return CovFactorResult(
            L=L_full,
            effective_rank=Gr,
            full_size=G,
            method="cholesky"
        )
    except la.LinAlgError:
        pass

    # --- 2b. Fallback: eigendecomposition + QR ---
    eigenvalues, V = la.eigh(Ar)

    max_eig = eigenvalues.max() if eigenvalues.max() > 0 else 1.0
    eigenvalues[eigenvalues / max_eig < tol] = 0.0
    eigenvalues[eigenvalues < 0.0] = 0.0

    pos_mask = eigenvalues > 0.0
    r = int(pos_mask.sum())

    if r == 0:
        return CovFactorResult(
            L=np.zeros((G, 0), dtype=np.float64),
            effective_rank=0,
            full_size=G,
            method="eigen_qr",
        )

    V_pos = V[:, pos_mask]
    sqrt_lam = np.sqrt(eigenvalues[pos_mask])
    Athin = V_pos * sqrt_lam[np.newaxis, :]

    _, R = la.qr(Athin.T, mode="economic")
    Lr = R.T

    L_full = np.zeros((G, r), dtype=np.float64)
    L_full[nz, :] = Lr

    return CovFactorResult(
        L=L_full,
        effective_rank=r,
        full_size=G,
        method="eigen_qr",
    )

# -----------------------------------------------------------------------------
# Data model + HDF5 I/O
# -----------------------------------------------------------------------------

@dataclass
class NeutronXSCovariances:
    """Neutron multigroup **relative** covariance data for cross sections (MF=33)."""

    name: str
    energy_grid_ev: np.ndarray              # shape (G+1,)
    reactions: Dict[int, Dict[str, Any]]    # MT -> parsed dict from ERRORR
    mat: Optional[int] = None
    temperature_k: Optional[float] = None   # what we processed at (single T)
    factor_results: Optional[Dict[int, Dict[int, CovFactorResult]]] = None 

    @classmethod
    def from_endf(
        cls,
        endf_path: str | Path,
        energy_grid_ev: Sequence[float],
        *,
        njoy_exec: Optional[str] = None,
        mat: Optional[int] = None,
        temperature: float = 293.6,
        name: Optional[str] = None,
        compute_factors: bool = True,
        eig_tol: float = 1e-10,
    ) -> "NeutronXSCovariances":
        res = generate_errorr_mf33(
            endf_path,
            energy_grid_ev,
            njoy_exec=njoy_exec,
            mat=mat,
            temperature=temperature,
        )
        mat_used = int(res["mat"])
        ek = np.asarray(res["ek"], dtype=float)

        reactions = parse_errorr_mf33_text(res["tape33"], mat=mat_used)

        # --- Parsed block inventory ---
        n_self, n_cross = 0, 0
        cross_pairs: List[Tuple[int, int]] = []
        for mt_key, sec in reactions.items():
            for mt1_key in sec.get("COVS", {}):
                if int(mt_key) == int(mt1_key):
                    n_self += 1
                else:
                    n_cross += 1
                    cross_pairs.append((int(mt_key), int(mt1_key)))
        log.info(
            "ERRORR tape33 parsed: %d MT section(s), "
            "%d self-covariance block(s), %d cross-covariance block(s)",
            len(reactions), n_self, n_cross,
        )
        if cross_pairs:
            for mt_a, mt1_a in cross_pairs:
                log.info("  cross-covariance: MT %d <-> MT %d", mt_a, mt1_a)
        else:
            log.info(
                "  No explicit cross-covariance blocks in ERRORR output.  "
                "Cross-reaction correlations (if any) are handled implicitly "
                "via NC-type derivation relations in the evaluation."
            )

        raw_validation = validate_raw_mf33_reactions(reactions, ek, atol=1e-14, rtol=1e-10)
        
        n_warnings = 0
        for mt, blocks in raw_validation.items():
            for mt1, qa in blocks.items():
                if qa.status == "fail":
                    raise ValueError(
                        f"Stage-1 MF=33 validation failed for MT {mt} -> MT1 {mt1}: "
                        + "; ".join(qa.messages)
                    )
                elif qa.status == "warn":
                    n_warnings += 1
                    log.info(
                        "Stage-1 MF=33 warning for MT %s -> MT1 %s: %s",
                        mt, mt1, "; ".join(qa.messages)
                    )
        # Pre-compute covariance factors for all sub-blocks
        n_eigen_qr = 0
        fcache: Optional[Dict[int, Dict[int, CovFactorResult]]] = None
        if compute_factors:
            fcache = {}
            for mt_key, sec in reactions.items():
                fcache[mt_key] = {}
                for mt1, M in sec.get("COVS", {}).items():
                    if int(mt_key) != int(mt1):
                        log.info(
                            "  MT %s -> MT1 %s: cross-block (%d x %d), "
                            "stored raw (no factorization)",
                            mt_key, mt1,
                            M.shape[0], M.shape[1] if M.ndim > 1 else 0,
                        )
                        continue  # cross-block: store raw only
                    result = compute_covariance_factor(
                        np.asarray(M, dtype=np.float64),
                        tol=eig_tol,
                    )
                    fcache[mt_key][mt1] = result
                    if result.method == "eigen_qr":
                        n_eigen_qr += 1
                    log.info(
                        "  MT %s -> MT1 %s: method=%-9s  rank=%d/%d  ",
                        mt_key, mt1,
                        result.method,
                        result.effective_rank,
                        result.full_size
                    )

        # --- One-line summary ---
        log.info(
            "%s: %d self-block(s), %d cross-block(s), "
            "%d eigen_qr fallback(s), %d warning(s)",
            name if name is not None else Path(endf_path).stem,
            n_self, n_cross, n_eigen_qr, n_warnings,
        )

        if name is None:
            name = Path(endf_path).stem

        return cls(
            name=str(name),
            energy_grid_ev=ek,
            reactions=reactions,
            mat=mat_used,
            temperature_k=float(temperature),
            factor_results=fcache,
        )
    
    # -----------------------------------------------------------------
    # HDF5 writing
    # -----------------------------------------------------------------

    def to_hdf5(self, filename: str | Path, store_raw_covariance: bool = True) -> None:
        """
        Write covariances to a standalone HDF5 file.
        """
        import h5py

        filename = Path(filename)
        filename.parent.mkdir(parents=True, exist_ok=True)

        with h5py.File(filename, "w") as f:
            self.write_mf33_group(f, store_raw_covariance=store_raw_covariance)

    def write_mf33_group(self, h5_group, store_raw_covariance: bool = True,
                         eig_tol: float = 1e-10) -> None:
        """
        Write the "mf33" sub-tree into an already-open HDF5 group.
        """
        if "mf33" in h5_group:
            del h5_group["mf33"]
        mf33 = h5_group.create_group("mf33")

        mf33.attrs["format"] = np.bytes_("openmc.mf33.v1")
        mf33.attrs["source"] = np.bytes_("njoy errorr")
        mf33.attrs["relative"] = 1  # int flag – portable across languages
        mf33.attrs["store_raw_covariance"] = int(store_raw_covariance)
        #mf33.attrs["factorization"] = np.bytes_("eigen_qr_thin")
        if self.mat is not None:
            mf33.attrs["mat"] = int(self.mat)
        if self.temperature_k is not None:
            mf33.attrs["temperature_k"] = float(self.temperature_k)

        mf33.create_dataset(
            "energy_grid_ev",
            data=np.asarray(self.energy_grid_ev, dtype=np.float64),
        )
        nonempty_mts = sorted(mt for mt, sec in self.reactions.items() if sec.get("COVS", {}))
        
        mf33.create_dataset(
            "mts",
            data=np.asarray(nonempty_mts, dtype=np.int32),
        )

        if store_raw_covariance:
            greact = mf33.create_group("reactions")
        gfact = mf33.create_group("factors")

        cached = self.factor_results or {}

        for mt, sec in self.reactions.items():
            if not sec.get("COVS", {}):
                continue 
            mt_str = str(int(mt))

            if store_raw_covariance:
                gmt = greact.create_group(mt_str)
                for attr_name in ("ZA", "AWR"):
                    if attr_name in sec:
                        gmt.attrs[attr_name] = float(sec[attr_name])

            covs: Dict[int, np.ndarray] = sec.get("COVS", {})
            for mt1, M in covs.items():
                M_arr = np.asarray(M, dtype=np.float64)
                ds_name = str(int(mt1))
                is_self = int(mt) == int(mt1)

                # ---- raw covariance ----
                # Self-blocks: optional. Cross-blocks: always stored (no factor).
                if store_raw_covariance or not is_self:
                    if not store_raw_covariance:
                        # Ensure the reactions group exists for cross-blocks
                        if "reactions" not in mf33:
                            greact = mf33.create_group("reactions")
                        if mt_str not in greact:
                            gmt = greact.create_group(mt_str)
                            for attr_name in ("ZA", "AWR"):
                                if attr_name in sec:
                                    gmt.attrs[attr_name] = float(sec[attr_name])
                        else:
                            gmt = greact[mt_str]
                    gmt.create_dataset(
                        ds_name,
                        data=M_arr,
                        compression="gzip",
                        shuffle=True,
                    )

                # ---- Triangular factor (self-blocks only) ----
                if is_self:
                    result = cached.get(mt, {}).get(mt1, None)
                    if result is None:
                        result = compute_covariance_factor(
                            M_arr, tol=eig_tol,)

                    if mt_str not in gfact:
                        gmt_fact = gfact.create_group(mt_str)
                        for attr_name in ("ZA", "AWR"):
                            if attr_name in sec:
                                gmt_fact.attrs[attr_name] = float(sec[attr_name])
                    else:
                        gmt_fact = gfact[mt_str]

                    gmt_fact.create_dataset(
                        ds_name,
                        data=result.L,
                        compression="gzip",
                        shuffle=True,
                    )

    # -----------------------------------------------------------------
    # HDF5 reading
    # -----------------------------------------------------------------
 
    @classmethod
    def _read_mf33_group(cls, mf33_group, name: str = "unknown") -> "NeutronXSCovariances":
        """Reconstruct a :class:`NeutronXSCovariances` from an open HDF5
        ``mf33`` group (the inverse of :meth:`write_mf33_group`).
 
        Parameters
        ----------
        mf33_group : h5py.Group
            The ``mf33`` group inside the HDF5 file.
        name : str
            Nuclide name to attach to the returned object.
 
        Returns
        -------
        NeutronXSCovariances
        """
        energy_grid_ev = mf33_group["energy_grid_ev"][()]
        mt_list = mf33_group["mts"][()]
 
        mat = int(mf33_group.attrs["mat"]) if "mat" in mf33_group.attrs else None
        temperature_k = (
            float(mf33_group.attrs["temperature_k"])
            if "temperature_k" in mf33_group.attrs
            else None
        )
 
        has_raw = "reactions" in mf33_group
        has_factors = "factors" in mf33_group
 
        reactions: Dict[int, Dict[str, Any]] = {}
        factor_results: Dict[int, Dict[int, CovFactorResult]] = {}
 
        for mt in mt_list:
            mt_str = str(int(mt))
            sec: Dict[str, Any] = {"MAT": mat, "MF": 33, "MT": int(mt)}
 
            # --- Read raw covariance matrices (optional) ---
            covs: Dict[int, np.ndarray] = {}
            if has_raw and mt_str in mf33_group["reactions"]:
                gmt = mf33_group["reactions"][mt_str]
                for attr_name in ("ZA", "AWR"):
                    if attr_name in gmt.attrs:
                        sec[attr_name] = float(gmt.attrs[attr_name])
                for ds_name in gmt:
                    mt1 = int(ds_name)
                    covs[mt1] = gmt[ds_name][()]
            sec["COVS"] = covs
 
            # --- Read triangular factors ---
            mt_factors: Dict[int, CovFactorResult] = {}
            if has_factors and mt_str in mf33_group["factors"]:
                gmt_fact = mf33_group["factors"][mt_str]
                for attr_name in ("ZA", "AWR"):
                    if attr_name in gmt_fact.attrs:
                        sec.setdefault(attr_name, float(gmt_fact.attrs[attr_name]))
                for ds_name in gmt_fact:
                    mt1 = int(ds_name)
                    L = gmt_fact[ds_name][()]
                    G = L.shape[0]
                    r = L.shape[1] if L.ndim == 2 else 0
                    mt_factors[mt1] = CovFactorResult(
                        L=L,
                        effective_rank=r,
                        full_size=G,
                        method="from_hdf5",
                    )
 
                    # If raw covariance was not stored, synthesise a
                    # minimal COVS entry so downstream code that iterates
                    # over reactions["COVS"] still works.
                    if mt1 not in covs:
                        covs[mt1] = L @ L.T
 
            reactions[int(mt)] = sec
            factor_results[int(mt)] = mt_factors
 
        return cls(
            name=name,
            energy_grid_ev=energy_grid_ev,
            reactions=reactions,
            mat=mat,
            temperature_k=temperature_k,
            factor_results=factor_results if factor_results else None,
        )
 
    @classmethod
    def from_hdf5(cls, filename: "str | Path") -> "NeutronXSCovariances":
        """Read covariances from a standalone HDF5 file written by
        :meth:`to_hdf5`.
 
        Parameters
        ----------
        filename : str or Path
            Path to the HDF5 file.
 
        Returns
        -------
        NeutronXSCovariances
        """
        import h5py
 
        filename = Path(filename)
        with h5py.File(filename, "r") as f:
            if "mf33" in f:
                mf33 = f["mf33"]
            else:
                raise KeyError(
                    f"HDF5 file '{filename}' does not contain an 'mf33' group."
                )
            name = filename.stem
            return cls._read_mf33_group(mf33, name=name)