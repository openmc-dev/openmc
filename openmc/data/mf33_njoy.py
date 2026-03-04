#!/usr/bin/env python3
"""
Minimal NJOY/ERRORR driver to produce MF=33 (multigroup covariance) as text.
Run a short NJOY chain and return dict with keys: mat, ek, tape33.

Assumptions
-----------
- `energy_grid_ev` is a strictly increasing sequence of group boundaries in **eV** (length G+1).
- Output is relative covariances (IRELCO=1) on an explicit group structure (IGN=1).
"""

from __future__ import annotations

import os
import shutil
import subprocess as sp
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Dict, Optional, Sequence, Any, Tuple, List


# ------------------------- small helpers -------------------------

def _read_mat_from_endf(endf_path: Path) -> int:
    """Infer MAT from an ENDF-6 evaluation (text)."""
    first_nonzero = None
    with endf_path.open("r", errors="ignore") as f:
        for line in f:
            if len(line) < 75:
                continue
            try:
                mat = int(line[66:70])
                mf = int(line[70:72])
                mt = int(line[72:75])
            except ValueError:
                continue
            if mat != 0 and first_nonzero is None:
                first_nonzero = mat
            if mat != 0 and mf == 1 and mt == 451:
                return mat
    if first_nonzero is None:
        raise ValueError(f"Could not infer MAT from ENDF file: {endf_path}")
    return first_nonzero


def _validate_energy_grid_ev(ek: Sequence[float]) -> List[float]:
    ek = [float(x) for x in ek]
    if len(ek) < 2:
        raise ValueError("Energy grid must have at least 2 boundaries (G+1).")

    for i in range(1, len(ek)):
        if not (ek[i] > ek[i-1]):
            raise ValueError("Energy grid boundaries must be strictly increasing (in eV).")

    return ek


# ------------------------- NJOY deck builder -------------------------

def _moder_input(nin: int, nout: int) -> str:
    return f"moder\n{nin:d} {nout:d} /\n"


def _reconr_input(endfin: int, pendfout: int, mat: int, tol: float, header: str) -> str:
    return (
        "reconr\n"
        f"{endfin:d} {pendfout:d} /\n"
        f"'{header}'/\n"
        f"{mat:d} 0 0 /\n"
        f"{tol:g} 0. /\n"
        "0/\n"
    )


def _broadr_input(endfin: int, pendfin: int, pendfout: int, mat: int, tol: float, temperature: float) -> str:
    return (
        "broadr\n"
        f"{endfin:d} {pendfin:d} {pendfout:d} /\n"
        f"{mat:d} 1 0 0 0. /\n"
        f"{tol:g} /\n"
        f"{temperature:.1f} /\n"
        "0 /\n"
    )


def _errorr_mf33_input(
    *,
    endfin: int,
    pendfin: int,
    errorrout: int,
    mat: int,
    temperature: float,
    ek: Sequence[float],
) -> str:
    # Fixed / stable defaults
    ign = 1       # explicit group structure
    iwt = 2       # constant weight
    iprint = 0    # quiet
    irelco = 1    # relative covariances
    irespr = 1    # recommended resonance parameter processing
    iread = 0     # don't restrict MT list (process all available in MF=33)

    nk = len(ek) - 1
    return (
        "errorr\n"
        f"{endfin:d} {pendfin:d} 0 {errorrout:d} 0 /\n"
        f"{mat:d} {ign:d} {iwt:d} {iprint:d} {irelco:d} /\n"
        f"{iprint:d} {temperature:.1f} /\n"
        f"{iread:d} 33 {irespr:d}/\n"
        f"{nk:d} /\n"
        + " ".join(f"{x:.5e}" for x in ek) + " /\n"
    )


def _build_deck(mat: int, ek: Sequence[float], temperature: float) -> str:
    # Common convention: negative units are formatted.
    # Use a minimal chain: MODER -> RECONR -> BROADR -> MODER -> ERRORR
    tol = 0.001
    e = 21  # internal endf copy
    p = 22  # pendf after reconr
    b = 23  # pendf after broadr
    out_errorr = 33

    deck = ""
    deck += _moder_input(20, -e)
    deck += _reconr_input(-e, -p, mat=mat, tol=tol, header="openmc_mf33")
    deck += _broadr_input(-e, -p, -b, mat=mat, tol=tol, temperature=temperature)
    deck += _errorr_mf33_input(endfin=-e, pendfin=-b, errorrout=out_errorr, mat=mat, temperature=temperature, ek=ek)
    deck += "stop\n"
    return deck


def _run_njoy(deck: str, endf: Path, exe: Optional[str]) -> Dict[int, str]:
    if exe is None:
        exe = os.environ.get("NJOY")
        if not exe:
            raise ValueError("NJOY executable not provided and NJOY env var is not set.")

    with TemporaryDirectory() as td:
        tmpdir = Path(td)
        shutil.copy(endf, tmpdir / "tape20")

        proc = sp.Popen(
            exe,
            shell=True,
            cwd=str(tmpdir),
            stdin=sp.PIPE,
            stdout=True, # sp.DEVNULL
            stderr=True, #sp.DEVNULL
        )
        proc.communicate(input=deck.encode())
        if proc.returncode != 0:
            raise RuntimeError(f"NJOY failed with return code {proc.returncode} (run with a debug wrapper to capture tapes).")

        tp33 = tmpdir / "tape33"
        if not tp33.exists():
            raise RuntimeError("NJOY ran but did not produce tape33 (ERRORR output).")

        return {33: tp33.read_text()}

def generate_errorr_mf33(
    endf_path: str | Path,
    energy_grid_ev: Sequence[float],
    *,
    njoy_exec: Optional[str] = None,
    mat: Optional[int] = None,
    temperature: float = 293.6,
) -> Dict[str, Any]:
    """Run NJOY/ERRORR and return MF=33 tape33 as text.

    Parameters
    ----------
    endf_path
        Path to ENDF-6 evaluation (text).
    energy_grid_ev
        Explicit group boundaries (G+1 values) in **eV**.
    njoy_exec
        NJOY executable/command. If omitted, uses $NJOY.
    mat
        MAT number. If omitted, inferred from the evaluation.
    temperature
        Processing temperature in K (single temperature).

    Returns
    -------
    dict
        Keys: "mat", "ek" (validated eV grid), "tape33" (text).
    """
    endf_path = Path(endf_path).resolve()
    if not endf_path.exists():
        raise FileNotFoundError(endf_path)

    ek = _validate_energy_grid_ev(energy_grid_ev)
    if mat is None:
        mat = _read_mat_from_endf(endf_path)

    deck = _build_deck(mat=mat, ek=ek, temperature=float(temperature))
    tapes = _run_njoy(deck=deck, endf=endf_path, exe=njoy_exec)

    return {"mat": int(mat), "ek": ek, "tape33": tapes[33]}
