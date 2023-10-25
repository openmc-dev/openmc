from pathlib import Path

import numpy as np

_FILES = (
    ('electron', 'electrons.txt'),
    ('helium', 'helium_ions.txt'),
    ('mu-', 'negative_muons.txt'),
    ('pi-', 'negative_pions.txt'),
    ('neutron', 'neutrons.txt'),
    ('photon', 'photons.txt'),
    ('photon kerma', 'photons_kerma.txt'),
    ('mu+', 'positive_muons.txt'),
    ('pi+', 'positive_pions.txt'),
    ('positron', 'positrons.txt'),
    ('proton', 'protons.txt')
)

_DOSE_TABLES = {'icrp116':{}}


def _load_dose(library='icrp116'):
    """Load effective dose tables from text files

    Parameters
    ----------
    library : {'icrp116'}
        The dose conversion library to use see `ICRP Publication 116 
        <https://doi.org/10.1016/j.icrp.2011.10.001>`_.

    """
    for particle, filename in _FILES:
        path = Path(__file__).parent / library / filename
        data = np.loadtxt(path, skiprows=3, encoding='utf-8')
        data[:, 0] *= 1e6   # Change energies to eV
        _DOSE_TABLES[library][particle] = data


def dose_coefficients(particle, geometry='AP', library='icrp116'):
    """Return effective dose conversion coefficients from ICRP-116

    This function provides fluence (and air kerma) to effective dose conversion
    coefficients for various types of external exposures based on values in ICRP
    publications. Corrected values found in a correigendum are used rather than
    the values in theoriginal report. Avaiable libraries include `ICRP
    Publication 116 <https://doi.org/10.1016/j.icrp.2011.10.001>`_

    Parameters
    ----------
    particle : {'neutron', 'photon', 'photon kerma', 'electron', 'positron'}
        Incident particle
    geometry : {'AP', 'PA', 'LLAT', 'RLAT', 'ROT', 'ISO'}
        Irradiation geometry assumed. Refer to ICRP-116 (Section 3.2) for the
        meaning of the options here.
    library : {'icrp116'}
        The dose conversion library to use.

    Returns
    -------
    energy : numpy.ndarray
        Energies at which dose conversion coefficients are given
    dose_coeffs : numpy.ndarray
        Effective dose coefficients in [pSv cm^2] at provided energies. For
        'photon kerma', the coefficients are given in [Sv/Gy].

    """
    if not _DOSE_TABLES[library]:
        _load_dose(library=library)

    # Get all data for selected particle
    data = _DOSE_TABLES[library].get(particle)
    if data is None:
        raise ValueError(f"{particle} has no effective dose data")

    # Determine index for selected geometry
    if particle in ('neutron', 'photon', 'proton', 'photon kerma'):
        index = ('AP', 'PA', 'LLAT', 'RLAT', 'ROT', 'ISO').index(geometry)
    else:
        index = ('AP', 'PA', 'ISO').index(geometry)

    # Pull out energy and dose from table
    energy = data[:, 0].copy()
    dose_coeffs = data[:, index + 1].copy()
    return energy, dose_coeffs
