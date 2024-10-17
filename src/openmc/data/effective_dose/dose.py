from pathlib import Path

import numpy as np

import openmc.checkvalue as cv

_FILES = {
    ('icrp74', 'neutron'): Path('icrp74') / 'neutrons.txt',
    ('icrp74', 'photon'): Path('icrp74') / 'photons.txt',
    ('icrp116', 'electron'): Path('icrp116') / 'electrons.txt',
    ('icrp116', 'helium'): Path('icrp116') / 'helium_ions.txt',
    ('icrp116', 'mu-'): Path('icrp116') / 'negative_muons.txt',
    ('icrp116', 'pi-'): Path('icrp116') / 'negative_pions.txt',
    ('icrp116', 'neutron'): Path('icrp116') / 'neutrons.txt',
    ('icrp116', 'photon'): Path('icrp116') / 'photons.txt',
    ('icrp116', 'photon kerma'): Path('icrp116') / 'photons_kerma.txt',
    ('icrp116', 'mu+'): Path('icrp116') / 'positive_muons.txt',
    ('icrp116', 'pi+'): Path('icrp116') / 'positive_pions.txt',
    ('icrp116', 'positron'): Path('icrp116') / 'positrons.txt',
    ('icrp116', 'proton'): Path('icrp116') / 'protons.txt',
}

_DOSE_TABLES = {}


def _load_dose_icrp(data_source: str, particle: str):
    """Load effective dose tables from text files.

    Parameters
    ----------
    data_source : {'icrp74', 'icrp116'}
        The dose conversion data source to use
    particle : {'neutron', 'photon', 'photon kerma', 'electron', 'positron'}
        Incident particle

    """
    path = Path(__file__).parent / _FILES[data_source, particle]
    data = np.loadtxt(path, skiprows=3, encoding='utf-8')
    data[:, 0] *= 1e6   # Change energies to eV
    _DOSE_TABLES[data_source, particle] = data


def dose_coefficients(particle, geometry='AP', data_source='icrp116'):
    """Return effective dose conversion coefficients.

    This function provides fluence (and air kerma) to effective or ambient dose
    (H*(10)) conversion coefficients for various types of external exposures
    based on values in ICRP publications. Corrected values found in a
    corrigendum are used rather than the values in the original report.
    Available libraries include `ICRP Publication 74
    <https://doi.org/10.1016/S0146-6453(96)90010-X>` and `ICRP Publication 116
    <https://doi.org/10.1016/j.icrp.2011.10.001>`.

    For ICRP 74 data, the photon effective dose per fluence is determined by
    multiplying the air kerma per fluence values (Table A.1) by the effective
    dose per air kerma (Table A.17). The neutron effective dose per fluence is
    found in Table A.41. For ICRP 116 data, the photon effective dose per
    fluence is found in Table A.1 and the neutron effective dose per fluence is
    found in Table A.5.

    Parameters
    ----------
    particle : {'neutron', 'photon', 'photon kerma', 'electron', 'positron'}
        Incident particle
    geometry : {'AP', 'PA', 'LLAT', 'RLAT', 'ROT', 'ISO'}
        Irradiation geometry assumed. Refer to ICRP-116 (Section 3.2) for the
        meaning of the options here.
    data_source : {'icrp74', 'icrp116'}
        The data source for the effective dose conversion coefficients.

    Returns
    -------
    energy : numpy.ndarray
        Energies at which dose conversion coefficients are given
    dose_coeffs : numpy.ndarray
        Effective dose coefficients in [pSv cm^2] at provided energies. For
        'photon kerma', the coefficients are given in [Sv/Gy].

    """

    cv.check_value('geometry', geometry, {'AP', 'PA', 'LLAT', 'RLAT', 'ROT', 'ISO'})
    cv.check_value('data_source', data_source, {'icrp74', 'icrp116'})

    if (data_source, particle) not in _FILES:
        raise ValueError(f"{particle} has no dose data in data source {data_source}.")
    elif (data_source, particle) not in _DOSE_TABLES:
        _load_dose_icrp(data_source, particle)

    # Get all data for selected particle
    data = _DOSE_TABLES[data_source, particle]

    # Determine index for selected geometry
    if particle in ('neutron', 'photon', 'proton', 'photon kerma'):
        columns = ('AP', 'PA', 'LLAT', 'RLAT', 'ROT', 'ISO')
    else:
        columns = ('AP', 'PA', 'ISO')
    index = columns.index(geometry)

    # Pull out energy and dose from table
    energy = data[:, 0].copy()
    dose_coeffs = data[:, index + 1].copy()
    return energy, dose_coeffs
