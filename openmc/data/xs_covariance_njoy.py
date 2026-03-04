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

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import scipy.linalg as la
from .mf33_njoy import generate_errorr_mf33


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
# Cholesky factor computation (with eigendecomposition fallback)
# -----------------------------------------------------------------------------

def _compute_cholesky_factor(cov_matrix: np.ndarray, tol: float = 1e-10) -> np.ndarray:
    """Compute the lower-triangular Cholesky factor L such that Σ ≈ L L^T.

    If the matrix is symmetric positive-definite, standard Cholesky is used.
    Otherwise, falls back to eigendecomposition with negative-eigenvalue
    clamping (matching the logic in sampling_algorithm.cpp):
        1. Eigendecompose: Σ = V diag(λ) V^T
        2. Zero out eigenvalues with λ/λ_max < tol
        3. Form A = V diag(√λ)
        4. QR-decompose A^T = Q R  →  L = R^T

    Parameters
    ----------
    cov_matrix : np.ndarray
        Square covariance matrix (G x G).
    tol : float
        Eigenvalue tolerance: eigenvalues with λ/λ_max < tol are zeroed.

    Returns
    -------
    np.ndarray
        Lower-triangular factor L (G x G).
    """
    try:
        return la.cholesky(cov_matrix, lower=True)
    except la.LinAlgError:
        pass

    eigenvalues, V = la.eigh(cov_matrix)
    max_eig = eigenvalues.max() if eigenvalues.max() > 0 else 1.0
    eigenvalues[eigenvalues / max_eig < tol] = 0.0
    eigenvalues[eigenvalues < 0.0] = 0.0

    A = V * np.sqrt(eigenvalues)[np.newaxis, :]
    _, R = la.qr(A.T, mode='economic')
    return R.T


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
    cholesky_factors: Optional[Dict[int, Dict[int, np.ndarray]]] = None  # MT -> MT1 -> L

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

        # Pre-compute Cholesky factors for all covariance sub-blocks
        chol: Dict[int, Dict[int, np.ndarray]] = {}
        for mt, sec in reactions.items():
            chol[mt] = {}
            for mt1, M in sec.get("COVS", {}).items():
                chol[mt][mt1] = _compute_cholesky_factor(
                    np.asarray(M, dtype=np.float64)
                )

        if name is None:
            name = Path(endf_path).stem

        return cls(
            name=str(name),
            energy_grid_ev=ek,
            reactions=reactions,
            mat=mat_used,
            temperature_k=float(temperature),
            cholesky_factors=chol,
        )

    def to_hdf5(self, filename: str | Path) -> None:
        """Write covariances to a standalone HDF5 file.

        Layout matches the "covariance/mf33" group that
        :meth:'IncidentNeutron.export_to_hdf5' produces, so that both
        paths share one schema and one reader.
        """
        import h5py

        filename = Path(filename)
        filename.parent.mkdir(parents=True, exist_ok=True)

        with h5py.File(filename, "w") as f:
            self.write_mf33_group(f)

    def write_mf33_group(self, h5_group) -> None:
        """Write the "mf33" sub-tree into an already-open HDF5 group.

        Parameters
        ----------
        h5_group : h5py.Group
            Parent group.  A child group called "mf33" will be created
            (or replaced) directly under it.
        """
        if "mf33" in h5_group:
            del h5_group["mf33"]
        mf33 = h5_group.create_group("mf33")

        mf33.attrs["format"] = np.bytes_("openmc.mf33.v1")
        mf33.attrs["source"] = np.bytes_("njoy errorr")
        mf33.attrs["relative"] = 1  # int flag – portable across languages
        if self.mat is not None:
            mf33.attrs["mat"] = int(self.mat)
        if self.temperature_k is not None:
            mf33.attrs["temperature_k"] = float(self.temperature_k)

        mf33.create_dataset(
            "energy_grid_ev",
            data=np.asarray(self.energy_grid_ev, dtype=np.float64),
        )

        greact = mf33.create_group("reactions")
        gchol = mf33.create_group("cholesky")

        # Compute Cholesky factors on-the-fly if not already cached
        chol = self.cholesky_factors or {}

        for mt, sec in self.reactions.items():
            gmt = greact.create_group(str(int(mt)))
            gmt.attrs["ZA"] = float(sec.get("ZA", 0.0))
            gmt.attrs["AWR"] = float(sec.get("AWR", 0.0))

            gmt_chol = gchol.create_group(str(int(mt)))
            gmt_chol.attrs["ZA"] = float(sec.get("ZA", 0.0))
            gmt_chol.attrs["AWR"] = float(sec.get("AWR", 0.0))

            covs: Dict[int, np.ndarray] = sec.get("COVS", {})
            for mt1, M in covs.items():
                M_arr = np.asarray(M, dtype=np.float64)
                gmt.create_dataset(
                    str(int(mt1)),
                    data=M_arr,
                    compression="gzip",
                    shuffle=True,
                )

                # Use cached L if available, otherwise compute
                L = chol.get(mt, {}).get(mt1, None)
                if L is None:
                    L = _compute_cholesky_factor(M_arr)
                gmt_chol.create_dataset(
                    str(int(mt1)),
                    data=np.asarray(L, dtype=np.float64),
                    compression="gzip",
                    shuffle=True,
                )

    @classmethod
    def from_hdf5(cls, filename: str | Path, name: Optional[str] = None) -> "NeutronXSCovariances":
        """Read covariances back from an HDF5 file.

        Supports two layouts:
        - **standalone** file written by :meth:`to_hdf5` (``mf33`` group
          at root)
        - **embedded** inside an OpenMC incident-neutron file
          (``<nuclide>/covariance/mf33``)
        """
        import h5py

        filename = Path(filename)
        with h5py.File(filename, "r") as f:
            mf33 = cls._find_mf33_group(f)
            return cls._read_mf33_group(mf33, name=name or filename.stem)

    @staticmethod
    def _find_mf33_group(f) -> "h5py.Group":
        """Locate the ``mf33`` group regardless of nesting depth."""
        # Standalone file: /mf33
        if "mf33" in f:
            return f["mf33"]
        # Embedded: /<nuclide>/covariance/mf33
        for key in f:
            g = f[key]
            if isinstance(g, type(f)) and "covariance" in g:  # h5py.Group
                cov = g["covariance"]
                if "mf33" in cov:
                    return cov["mf33"]
        raise KeyError("No mf33 group found in HDF5 file.")

    @classmethod
    def _read_mf33_group(cls, mf33, *, name: str = "") -> "NeutronXSCovariances":
        """Build an instance from an ``mf33`` h5py.Group."""
        mat_val = mf33.attrs.get("mat", None)
        temp_val = mf33.attrs.get("temperature_k", None)

        ek = np.asarray(mf33["energy_grid_ev"][...], dtype=np.float64)

        reactions: Dict[int, Dict[str, Any]] = {}
        greact = mf33["reactions"]
        for mt_str, gmt in greact.items():
            mt = int(mt_str)
            covs: Dict[int, np.ndarray] = {}
            for mt1_str, ds in gmt.items():
                covs[int(mt1_str)] = np.asarray(ds[...], dtype=np.float64)

            reactions[mt] = {
                "MAT": int(mat_val) if mat_val is not None else 0,
                "MF": 33,
                "MT": mt,
                "ZA": float(gmt.attrs.get("ZA", 0.0)),
                "AWR": float(gmt.attrs.get("AWR", 0.0)),
                "COVS": covs,
            }

        # Read pre-computed Cholesky factors if present
        chol: Optional[Dict[int, Dict[int, np.ndarray]]] = None
        if "cholesky" in mf33:
            chol = {}
            for mt_str, gmt_chol in mf33["cholesky"].items():
                mt = int(mt_str)
                chol[mt] = {}
                for mt1_str, ds in gmt_chol.items():
                    chol[mt][int(mt1_str)] = np.asarray(ds[...], dtype=np.float64)

        return cls(
            name=str(name),
            energy_grid_ev=ek,
            reactions=reactions,
            mat=int(mat_val) if mat_val is not None else None,
            temperature_k=float(temp_val) if temp_val is not None else None,
            cholesky_factors=chol,
        )