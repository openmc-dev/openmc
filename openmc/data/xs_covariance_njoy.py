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
# Covariance factor computation (eigendecomposition + QR)
# -----------------------------------------------------------------------------

@dataclass
class CovFactorResult:
    """Container for a covariance factorization and its diagnostics.

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
            condition_number=np.inf,
            nonzero_variance_indices=nz,
        )

    Gr = Ar.shape[0]

    # --- 2a. Fast path: standard Cholesky ---
    try:
        Lr = la.cholesky(Ar, lower=True)
        eig_pos = np.diag(Lr) ** 2
        cond = float(eig_pos.max() / eig_pos.min()) if eig_pos.min() > 0 else np.inf

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

        # Pre-compute covariance factors for all sub-blocks
        fcache: Optional[Dict[int, Dict[int, CovFactorResult]]] = None
        if compute_factors:
            fcache = {}
            for mt_key, sec in reactions.items():
                fcache[mt_key] = {}
                for mt1, M in sec.get("COVS", {}).items():
                    result = compute_covariance_factor(
                        np.asarray(M, dtype=np.float64),
                        tol=eig_tol,
                    )
                    fcache[mt_key][mt1] = result
                    log.info(
                        "  MT %s -> MT1 %s: method=%-9s  rank=%d/%d  "
                        "cond=%.2e  neg_mass=%.2e%%",
                        mt_key, mt1,
                        result.method,
                        result.effective_rank,
                        result.full_size,
                        result.condition_number,
                        result.negative_eig_mass_pct,
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
        mf33.attrs["factorization"] = np.bytes_("eigen_qr_thin")
        if self.mat is not None:
            mf33.attrs["mat"] = int(self.mat)
        if self.temperature_k is not None:
            mf33.attrs["temperature_k"] = float(self.temperature_k)

        mf33.create_dataset(
            "energy_grid_ev",
            data=np.asarray(self.energy_grid_ev, dtype=np.float64),
        )
        mf33.create_dataset(
            "mts",
            data=np.asarray(sorted(self.reactions.keys()), dtype=np.int32),
        )

        if store_raw_covariance:
            greact = mf33.create_group("reactions")
        gfact = mf33.create_group("factors")

        cached = self.factor_results or {}

        for mt, sec in self.reactions.items():
            mt_str = str(int(mt))

            if store_raw_covariance:
                gmt = greact.create_group(mt_str)
                for attr_name in ("ZA", "AWR"):
                    if attr_name in sec:
                        gmt.attrs[attr_name] = float(sec[attr_name])

            gmt_fact = gfact.create_group(mt_str)
            for attr_name in ("ZA", "AWR"):
                if attr_name in sec:
                    gmt_fact.attrs[attr_name] = float(sec[attr_name])

            covs: Dict[int, np.ndarray] = sec.get("COVS", {})
            for mt1, M in covs.items():
                M_arr = np.asarray(M, dtype=np.float64)
                ds_name = str(int(mt1))

                # ---- raw covariance (optional) ----
                if store_raw_covariance:
                    gmt.create_dataset(
                        ds_name,
                        data=M_arr,
                        compression="gzip",
                        shuffle=True,
                    )

                # ---- Triangular factor (always stored) ----
                result = cached.get(mt, {}).get(mt1, None)
                if result is None:
                    result = compute_covariance_factor(
                        M_arr, tol=eig_tol,)

                gmt_fact.create_dataset(
                    ds_name,
                    data=result.L,
                    compression="gzip",
                    shuffle=True,
                )