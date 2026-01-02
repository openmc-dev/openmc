import numpy as np

import openmc.checkvalue as cv
from openmc.data import EV_PER_MEV

# Embedded NIST-126 data
# Air (Dry Near Sea Level) — NIST Standard Reference Database 126 Table 4 (doi: 10.18434/T4D01F)
# Columns: Energy (MeV), μen/ρ (cm^2/g)
_NIST126_AIR = np.array(
    [
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
    ],
    dtype=float,
)

# Registry of embedded tables: (data_source, material) -> ndarray
# Table shape: (N, 2) with columns [Energy (MeV), μen/ρ (cm^2/g)]
_MUEN_TABLES = {
    ("nist126", "air"): _NIST126_AIR,
}


def mu_en_coefficients(
    material: str, data_source: str = "nist126"
) -> tuple[np.ndarray, np.ndarray]:
    """Return tabulated mass energy-absorption coefficients.

    Parameters
    ----------
    material : {'air'}
        Material compound for which to load coefficients.
    data_source : {'nist126'}
        Source library.

    Returns
    -------
    energy : numpy.ndarray
        Energies [eV]
    mu_en_coeffs : numpy.ndarray
        Mass energy-absorption coefficients [cm^2/g]
    """
    cv.check_value("material", material, {"air"})
    cv.check_value("data_source", data_source, {"nist126"})

    key = (data_source, material)
    if key not in _MUEN_TABLES:
        available = sorted({m for (ds, m) in _MUEN_TABLES.keys() if ds == data_source})
        raise ValueError(
            f"'{material}' has no embedded μen/ρ table for data source {data_source}. "
            f"Available materials for {data_source}: {available}"
        )

    data = _MUEN_TABLES[key]
    energy = data[:, 0].copy() * EV_PER_MEV  # MeV -> eV
    mu_en_coeffs = data[:, 1].copy()
    return energy, mu_en_coeffs
