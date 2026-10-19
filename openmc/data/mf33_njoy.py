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
from typing import Any, Dict, List, Optional, Sequence, Union


NJOY_tolerance = 0.01
NJOY_sig0 = [1e10]
NJOY_thermr_emax = 10.0

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
            if mat > 0 and first_nonzero is None:
                first_nonzero = mat
            if mat > 0 and mf == 1 and mt == 451:
                return mat
    if first_nonzero is None:
        raise ValueError(f"Could not infer MAT from ENDF file: {endf_path}")
    return first_nonzero

def _endf_has_mf33(endf_path: Path, mat: int) -> bool:
    """Return True if the ENDF tape contains any MF=33 section for this MAT."""
    with endf_path.open("r", errors="ignore") as f:
        for line in f:
            if len(line) < 75:
                continue
            try:
                m  = int(line[66:70])
                mf = int(line[70:72])
            except ValueError:
                continue
            if m == mat and mf == 33:
                return True
    return False

def _validate_energy_grid_ev(ek: Sequence[float]) -> List[float]:
    ek = [float(x) for x in ek]
    if len(ek) < 2:
        raise ValueError("Energy grid must have at least 2 boundaries (G+1).")
    if ek[0] <= 0.0:
        raise ValueError(
            f"Energy grid lower boundary must be positive (got {ek[0]:g} eV). "
            f"ERRORR cannot integrate from zero energy. Use a small positive "
            f"value like 1e-5 eV instead."
        )
    for i in range(1, len(ek)):
        if not (ek[i] > ek[i-1]):
            raise ValueError("Energy grid boundaries must be strictly increasing (in eV).")
    return ek

# ------------------------- NJOY deck builder -------------------------

def _moder_input(nin: int, nout: int) -> str:
    return f"moder\n{nin:d} {nout:d} /\n"

def _reconr_input(endfin, pendfout, mat, err=NJOY_tolerance):
    return (
        "reconr\n"
        f"{endfin:d} {pendfout:d} /\n"
        f"'mf33'/\n"
        f"{mat:d} 0 0 /\n"
        f"{err:g} 0. /\n"
        "0/\n"
    )

def _thermr_input(endfin, pendfin, pendfout, mat, temperature,
                   err=NJOY_tolerance, emax=NJOY_thermr_emax):
    return (
        "thermr\n"
        f"{endfin:d} {pendfin:d} {pendfout:d} /\n"
        f"0 {mat:d} 20 1 1 0 0 1 221 0 /\n"
        f"{temperature:.1f} /\n"
        f"{err:g} {emax:g} /\n"
    )

def _broadr_input(endfin, pendfin, pendfout, mat, temperature, err=NJOY_tolerance):
    return (
        "broadr\n"
        f"{endfin:d} {pendfin:d} {pendfout:d} /\n"
        f"{mat:d} 1 0 0 0. /\n"
        f"{err:g} /\n"
        f"{temperature:.1f} /\n"
        "0 /\n"
    )

def _unresr_input(endfin, pendfin, pendfout, mat, temperature,
                   sig0=NJOY_sig0):
    nsig0 = len(sig0)
    return (
        "unresr\n"
        f"{endfin:d} {pendfin:d} {pendfout:d} /\n"
        f"{mat:d} 1 {nsig0:d} 0 /\n"
        f"{temperature:.1f} /\n"
        + " ".join(f"{s:.2E}" for s in sig0) + " /\n"
        "0 /\n"
    )

def _purr_input(endfin, pendfin, pendfout, mat, temperature,
                 sig0=NJOY_sig0, bins=20, ladders=32):
    nsig0 = len(sig0)
    return (
        "purr\n"
        f"{endfin:d} {pendfin:d} {pendfout:d} /\n"
        f"{mat:d} 1 {nsig0:d} {bins:d} {ladders:d} 0 /\n"
        f"{temperature:.1f} /\n"
        + " ".join(f"{s:.2E}" for s in sig0) + " /\n"
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
    iwt: int = 6,
    spectrum: Optional[Sequence[float]] = None,
    relative: bool = True,
    irespr: int = 1,
    mt: Optional[Union[int, Sequence[int]]] = None,
    iprint: bool = False,
) -> str:
    irelco = 1 if relative else 0
    iread  = 1 if mt is not None else 0
    iwt_   = 1 if spectrum is not None else iwt
    pf     = int(iprint)
    nk     = len(ek) - 1

    lines = [
        "errorr",
        f"{endfin:d} {pendfin:d} 0 {errorrout:d} 0 /",
        f"{mat:d} 1 {iwt_:d} {pf:d} {irelco} /",   # ign=1 (explicit grid)
        f"{pf:d} {temperature:.1f} /",
        f"{iread:d} 33 {irespr:d} /",
    ]

    # MT list
    if iread == 1:
        mtlist = [mt] if isinstance(mt, int) else list(mt)
        lines.append(f"{len(mtlist):d} 0 /")
        lines.append(" ".join(str(m) for m in mtlist) + " /")

    # Explicit energy grid
    lines.append(f"{nk} /")
    lines.append(" ".join(f"{x:.5e}" for x in ek) + " /")

    # User weight spectrum (TAB1, iwt=1)
    if iwt_ == 1 and spectrum is not None:
        n_pairs = len(spectrum) // 2
        lines.append(f" 0.0 0.0 0 0 1 {n_pairs:d} /")
        lines.append(f"{n_pairs:d} 1 /")
        for k in range(n_pairs):
            lines.append(f"{spectrum[2*k]:.6e} {spectrum[2*k+1]:.6e} /")
        lines.append("/")

    return "\n".join(lines) + "\n"


# ---- deck assembly ---------------------------------------------------------

def _build_deck(
    mat: int,
    ek: Sequence[float],
    temperature: float,
    *,
    err: float = NJOY_tolerance,
    thermr: bool = False,
    unresr: bool = False,
    purr: bool = False,
    iwt: int = 2,
    spectrum: Optional[Sequence[float]] = None,
    relative: bool = True,
    irespr: int = 1,
    mt: Optional[Union[int, Sequence[int]]] = None,
    iprint: bool = False,
) -> str:
    e = 21      # internal formatted ENDF
    p = e + 1   # running PENDF tape number
    deck = _moder_input(20, -e)

    # RECONR
    deck += _reconr_input(-e, -p, mat=mat, err=err)

    # BROADR
    if temperature > 0:
        o = p + 1
        deck += _broadr_input(-e, -p, -o, mat=mat,
                               temperature=temperature, err=err)
        p = o

    # THERMR (free-gas)
    if thermr:
        o = p + 1
        deck += _thermr_input(0, -p, -o, mat=mat,
                               temperature=temperature, err=err)
        p = o

    # UNRESR
    if unresr:
        o = p + 1
        deck += _unresr_input(-e, -p, -o, mat=mat,
                               temperature=temperature)
        p = o

    # PURR
    if purr:
        o = p + 1
        deck += _purr_input(-e, -p, -o, mat=mat,
                             temperature=temperature)
        p = o

    # ERRORR (MF=33)
    deck += _errorr_mf33_input(
        endfin=-e, pendfin=-p, errorrout=33,
        mat=mat, temperature=temperature, ek=ek,
        iwt=iwt, spectrum=spectrum, relative=relative,
        irespr=irespr, mt=mt, iprint=iprint,
    )
    deck += "stop\n"
    return deck


# ---- NJOY runner -----------------------------------------------------------

def _run_njoy(deck: str, endf: Path, exe: Optional[str]) -> Dict[int, str]:
    if exe is None:
        exe = os.environ.get("NJOY")
        if not exe:
            raise ValueError(
                "NJOY executable not provided and $NJOY env var is unset."
            )

    with TemporaryDirectory() as td:
        tmpdir = Path(td)
        shutil.copy(endf, tmpdir / "tape20")
        (tmpdir / "input").write_text(deck)

        proc = sp.Popen(
            exe, shell=True, cwd=str(tmpdir),
            stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE,
        )
        stdout, stderr = proc.communicate(input=deck.encode())

        if proc.returncode != 0:
            # Save work dir for debugging
            dest = Path("njoy_failed_outputs")
            if dest.exists():
                shutil.rmtree(dest)
            shutil.copytree(tmpdir, dest)
            raise RuntimeError(
                f"NJOY failed (rc={proc.returncode}). "
                f"Work dir saved to '{dest}'.\n"
                f"stderr: {stderr.decode(errors='replace')[:2000]}"
            )

        tp33 = tmpdir / "tape33"
        if not tp33.exists():
            raise RuntimeError("NJOY ran but did not produce tape33.")

        return {33: tp33.read_text()}


# ---- public API ------------------------------------------------------------

def generate_errorr_mf33(
    endf_path: str | Path,
    energy_grid_ev: Sequence[float],
    *,
    njoy_exec: Optional[str] = None,
    mat: Optional[int] = None,
    temperature: float = 293.6,
    # processing chain
    err: float = NJOY_tolerance,
    minimal_processing: bool = False,
    thermr: bool = False,
    unresr: bool = True,
    purr: bool = True,
    # ERRORR options
    iwt: int = 6,
    spectrum: Optional[Sequence[float]] = None,
    relative: bool = True,
    irespr: int = 1,
    mt: Optional[Union[int, Sequence[int]]] = None,
    iprint: bool = False,
) -> Dict[str, Any]:
    """Run NJOY/ERRORR and return MF=33 tape33 as text.

    Parameters
    ----------
    endf_path
        Path to ENDF-6 evaluation (text).
    energy_grid_ev
        Explicit group boundaries (G+1 values) in eV.
    njoy_exec
        NJOY executable/command.  Falls back to $NJOY.
    mat
        MAT number.  Inferred from the file when omitted.
    temperature
        Processing temperature in K.
    err
        Reconstruction tolerance for RECONR / BROADR.
    minimal_processing
        If True (default), force thermr/unresr/purr off.
    thermr, unresr, purr
        Individual module toggles.  All ignored when
        ``minimal_processing=True``.
    iwt
        Weight function option (default 2 = constant).
        Ignored when *spectrum* is given.
    spectrum
        User weight function as flat [E0, W0, E1, W1, …].
        Forces iwt=1.
    relative
        Produce relative covariances (default True).
    irespr
        Resonance-parameter covariance processing
        (0 = area sensitivity, 1 = 1%% sensitivity).
    mt
        Restrict ERRORR to specific reaction numbers.
    iprint
        NJOY print flag.

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

    if not _endf_has_mf33(endf_path, mat):
        raise ValueError(
            f"ENDF file {endf_path.name} (MAT={mat}) contains no MF=33 "
            f"covariance data; nothing for ERRORR to process."
        )

    if minimal_processing:
        thermr = unresr = purr = False

    deck = _build_deck(
        mat=mat, ek=ek, temperature=float(temperature),
        err=err, thermr=thermr, unresr=unresr, purr=purr,
        iwt=iwt, spectrum=spectrum, relative=relative,
        irespr=irespr, mt=mt, iprint=iprint,
    )
    print("Running NJOY to produce MF33 covariance data")
    tapes = _run_njoy(deck, endf_path, exe=njoy_exec)

    return {"mat": int(mat), "ek": ek, "tape33": tapes[33]}