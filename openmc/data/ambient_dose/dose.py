from pathlib import Path

import numpy as np

import openmc.checkvalue as cv

_FILES = {
    ('icrp74', 'neutron'): Path('icrp74') / 'neutrons.txt',
    ('icrp74', 'photon'): Path('icrp74') / 'photons.txt',
}

_DOSE_TABLES = {}


def _load_dose_icrp(data_source: str, particle: str):
    """Load effective dose tables from text files.

    Parameters
    ----------
    data_source : {'icrp74'}
        The dose conversion data source to use
    particle : {'neutron', 'photon'}
        Incident particle

    """
    path = Path(__file__).parent / _FILES[data_source, particle]
    data = np.loadtxt(path, skiprows=3, encoding='utf-8')
    data[:, 0] *= 1e6   # Change energies to eV
    _DOSE_TABLES[data_source, particle] = data


def ambient_dose_coefficients(particle, data_source='icrp74'):
    """Return effective dose conversion coefficients.

    This function provides fluence to ambient dose
    (H*(10)) conversion coefficients for various types of external exposures
    based on values in ICRP publications. Corrected values found in a
    corrigendum are used rather than the values in the original report.
    Available libraries include `ICRP Publication 74
    <https://doi.org/10.1016/S0146-6453(96)90010-X>` and `ICRP Publication 116
    <https://doi.org/10.1016/j.icrp.2011.10.001>`.

    Parameters
    ----------
    particle : {'neutron', 'photon'}
        Incident particle
    data_source : {'icrp74'}
        The data source for the effective dose conversion coefficients.

    Returns
    -------
    energy : numpy.ndarray
        Energies at which dose conversion coefficients are given
    dose_coeffs : numpy.ndarray
        Effective dose coefficients in [pSv cm^2] at provided energies. For
        'photon kerma', the coefficients are given in [Sv/Gy].

    """

    cv.check_value('data_source', data_source, {'icrp74'})

    if (data_source, particle) not in _FILES:
        raise ValueError(f"{particle} has no dose data in data source {data_source}.")
    elif (data_source, particle) not in _DOSE_TABLES:
        _load_dose_icrp(data_source, particle)

    # Get all data for selected particle
    data = _DOSE_TABLES[data_source, particle]

    # Pull out energy and dose from table
    energy = data[:, 0].copy()
    dose_coeffs = data[:, 1].copy()
    return energy, dose_coeffs
