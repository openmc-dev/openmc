from pathlib import Path

import numpy as np
import h5py

import openmc.checkvalue as cv
from openmc.data import EV_PER_MEV
from ..data import ATOMIC_NUMBER
from ..function import Tabulated1D

# Embedded NIST-126 data
# Air (Dry Near Sea Level) — NIST Standard Reference Database 126 Table 4 (doi: 10.18434/T4D01F)
# Columns: Energy (MeV), μ_en/ρ (cm^2/g)
_NIST126_AIR = np.array([
    [1.00000e-03, 3.599e03],
    [1.50000e-03, 1.188e03],
    [2.00000e-03, 5.262e02],
    [3.00000e-03, 1.614e02],
    [3.20290e-03, 1.330e02],
    [3.20290e-03, 1.460e02],
    [4.00000e-03, 7.636e01],
    [5.00000e-03, 3.931e01],
    [6.00000e-03, 2.270e01],
    [8.00000e-03, 9.446e00],
    [1.00000e-02, 4.742e00],
    [1.50000e-02, 1.334e00],
    [2.00000e-02, 5.389e-01],
    [3.00000e-02, 1.537e-01],
    [4.00000e-02, 6.833e-02],
    [5.00000e-02, 4.098e-02],
    [6.00000e-02, 3.041e-02],
    [8.00000e-02, 2.407e-02],
    [1.00000e-01, 2.325e-02],
    [1.50000e-01, 2.496e-02],
    [2.00000e-01, 2.672e-02],
    [3.00000e-01, 2.872e-02],
    [4.00000e-01, 2.949e-02],
    [5.00000e-01, 2.966e-02],
    [6.00000e-01, 2.953e-02],
    [8.00000e-01, 2.882e-02],
    [1.00000e00, 2.789e-02],
    [1.25000e00, 2.666e-02],
    [1.50000e00, 2.547e-02],
    [2.00000e00, 2.345e-02],
    [3.00000e00, 2.057e-02],
    [4.00000e00, 1.870e-02],
    [5.00000e00, 1.740e-02],
    [6.00000e00, 1.647e-02],
    [8.00000e00, 1.525e-02],
    [1.00000e01, 1.450e-02],
    [1.50000e01, 1.353e-02],
    [2.00000e01, 1.311e-02],
])

# Registry of embedded tables: (data_source, material) -> ndarray
# Table shape: (N, 2) with columns [Energy (MeV), μen/ρ (cm^2/g)]
_MUEN_TABLES = {
    ("nist126", "air"): _NIST126_AIR,
}


def mass_energy_absorption_coefficient(
    material: str, data_source: str = "nist126"
) -> Tabulated1D:
    r"""Return the mass energy-absorption coefficient as a function of energy.

    The mass energy-absorption coefficient, :math:`\mu_\text{en}/\rho`, is
    defined as the fraction of incident photon energy absorbed in a material per
    unit mass less the energy carried away by scattered photons. It is obtained
    from `NIST Standard Reference Database 126
    <https://doi.org/10.18434/T4D01F>`_: X-Ray Mass Attenuation Coefficients.

    Parameters
    ----------
    material : {'air'}
        Material compound for which to load coefficients.
    data_source : {'nist126'}
        Source library.

    Returns
    -------
    Tabulated1D
        Mass energy-absorption coefficient [cm^2/g] as a function of photon
        energy [eV], using log-log interpolation.

    """
    cv.check_value("material", material, {"air"})
    cv.check_value("data_source", data_source, {"nist126"})

    key = (data_source, material)
    if key not in _MUEN_TABLES:
        available = sorted({m for (ds, m) in _MUEN_TABLES.keys() if ds == data_source})
        raise ValueError(
            f"No mass energy-absorption data for '{material}' in data source "
            f"'{data_source}'. Available materials: {available}"
        )

    data = _MUEN_TABLES[key]
    energy = data[:, 0].copy() * EV_PER_MEV  # MeV -> eV
    mu_en_coeffs = data[:, 1].copy()
    return Tabulated1D(energy, mu_en_coeffs,
                       breakpoints=[len(energy)], interpolation=[5])


# Used in mass_attenuation_coefficient function as a cache.
# Maps atomic number Z (int) -> Tabulated1D of (mu/rho) [cm^2/g] vs E [eV]
_MASS_ATTENUATION: dict[int, object] = {}


def mass_attenuation_coefficient(element):
    r"""Return the photon mass attenuation coefficient as a function of energy.

    The mass energy-absorption coefficient, :math:`\mu_\text{en}/\rho`, is
    defined as the fraction of incident photon energy absorbed in a material per
    unit mass. Values for each element are obtained from `NIST Standard
    Reference Database 8 <https://doi.org/10.18434/T48G6X>`_: XCOM Photon Cross
    Sections Database.

    Parameters
    ----------
    element : str or int
        Element symbol (e.g., 'Fe') or atomic number (e.g., 26).

    Returns
    -------
    Tabulated1D
        Mass attenuation coefficient [cm^2/g] as a function of photon energy
        [eV], using log-log interpolation.

    """
    if not _MASS_ATTENUATION:
        data_file = Path(__file__).with_name('mass_attenuation.h5')
        with h5py.File(data_file, 'r') as f:
            for key, dataset in f.items():
                energies, mu_rho = dataset[()]  # shape (2, N)
                _MASS_ATTENUATION[int(key)] = Tabulated1D(
                    energies, mu_rho,
                    breakpoints=[len(energies)],
                    interpolation=[5]  # log-log
                )

    # Resolve element argument to atomic number
    if isinstance(element, str):
        if element not in ATOMIC_NUMBER:
            raise ValueError(f"'{element}' is not a recognized element symbol")
        Z = ATOMIC_NUMBER[element]
    else:
        Z = int(element)

    if Z not in _MASS_ATTENUATION:
        raise ValueError(f"No mass attenuation data available for Z={Z}")

    return _MASS_ATTENUATION[Z]
