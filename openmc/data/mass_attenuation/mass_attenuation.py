from pathlib import Path

import numpy as np

import openmc.checkvalue as cv

_FILES = {
    ('nist126', 'air'): Path('nist126') / 'air.txt',
    ('nist126', 'water'): Path('nist126') / 'water.txt',
}

_MU_TABLES = {}


def _load_mass_attenuation(data_source: str, material: str):
    """Load mass energy attenuation and absorption coefficients from 
    the NIST database stored in the text files.

    Parameters
    ----------
    data_source : {'nist126'}
        The data source to use for the mass attenuation coefficients.
    material : {'air', 'water'}
        Material compound for which to load mass attenuation.

    """
    path = Path(__file__).parent / _FILES[data_source, material]
    data = np.loadtxt(path, skiprows=5, encoding='utf-8')
    data[:, 0] *= 1e6   # Change energies to eV
    _MU_TABLES[data_source, material] = data


def mu_en_coefficients(material, data_source='nist126'):
    """Return mass energy-absorption coefficients.

    This function returns the phtono mass energy-absorption coefficients for 
    various tabulated material compounds.
    Available libraries include `NIST Standard Reference Database 126
    <https://dx.doi.org/10.18434/T4D01F>`.


    Parameters
    ----------
    material : {'air', 'water'}
        Material compound for which to load mass attenuation.
    data_source : {'nist126'}
        The data source to use for the mass attenuation coefficients.

    Returns
    -------
    energy : numpy.ndarray
        Energies at which mass energy-absorption coefficients are given.
    mu_en_coeffs : numpy.ndarray
        mass energy absoroption coefficients [cm^2/g] at provided energies.

    """

    cv.check_value('material', material, {'air','water'})
    cv.check_value('data_source', data_source, {'nist126'})

    if (data_source, material) not in _FILES:
        available_materials = sorted({m for (ds, m) in _FILES if ds == data_source})
        msg = (
            f"'{material}' has no mass energy-absorption coefficients in data source {data_source}. "
            f"Available materials for {data_source} are: {available_materials}"
        )
        raise ValueError(msg)
    elif (data_source, material) not in _MU_TABLES:
        _load_mass_attenuation(data_source, material)

    # Get all data for selected material
    data = _MU_TABLES[data_source, material]

    # mass energy-absorption coefficients are in the third column
    mu_en_index = 2

    # Pull out energy and dose from table
    energy = data[:, 0].copy()
    mu_en_coeffs = data[:, mu_en_index].copy()
    return energy, mu_en_coeffs
