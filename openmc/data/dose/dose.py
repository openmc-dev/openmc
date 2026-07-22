from pathlib import Path

import numpy as np

import openmc.checkvalue as cv

_FULL_GEOMETRIES = ('AP', 'PA', 'LLAT', 'RLAT', 'ROT', 'ISO')
_LIMITED_GEOMETRIES = ('AP', 'PA', 'ISO')

_TABLES = {
    ('icrp74', 'effective', 'neutron'): (
        Path('icrp74') / 'neutrons.txt', _FULL_GEOMETRIES),
    ('icrp74', 'effective', 'photon'): (
        Path('icrp74') / 'photons.txt', _FULL_GEOMETRIES),
    ('icrp116', 'effective', 'electron'): (
        Path('icrp116') / 'electrons.txt', _LIMITED_GEOMETRIES),
    ('icrp116', 'effective', 'helium'): (
        Path('icrp116') / 'helium_ions.txt', _LIMITED_GEOMETRIES),
    ('icrp116', 'effective', 'mu-'): (
        Path('icrp116') / 'negative_muons.txt', _LIMITED_GEOMETRIES),
    ('icrp116', 'effective', 'pi-'): (
        Path('icrp116') / 'negative_pions.txt', _LIMITED_GEOMETRIES),
    ('icrp116', 'effective', 'neutron'): (
        Path('icrp116') / 'neutrons.txt', _FULL_GEOMETRIES),
    ('icrp116', 'effective', 'photon'): (
        Path('icrp116') / 'photons.txt', _FULL_GEOMETRIES),
    ('icrp116', 'effective', 'photon kerma'): (
        Path('icrp116') / 'photons_kerma.txt', _FULL_GEOMETRIES),
    ('icrp116', 'effective', 'mu+'): (
        Path('icrp116') / 'positive_muons.txt', _LIMITED_GEOMETRIES),
    ('icrp116', 'effective', 'pi+'): (
        Path('icrp116') / 'positive_pions.txt', _LIMITED_GEOMETRIES),
    ('icrp116', 'effective', 'positron'): (
        Path('icrp116') / 'positrons.txt', _LIMITED_GEOMETRIES),
    ('icrp116', 'effective', 'proton'): (
        Path('icrp116') / 'protons.txt', _FULL_GEOMETRIES),
    ('icrp74', 'ambient', 'neutron'): (
        Path('icrp74') / 'neutrons_H10.txt', None),
    ('icrp74', 'ambient', 'photon'): (
        Path('icrp74') / 'photons_H10.txt', None),
}

_DOSE_TABLES = {}


def _load_dose_table(data_source: str, dose_quantity: str, particle: str):
    """Load dose tables from text files.

    Parameters
    ----------
    data_source : {'icrp74', 'icrp116'}
        The dose conversion data source to use
    dose_quantity : {'effective', 'ambient'}
        Dose quantity to load. 'ambient' corresponds to ambient dose
        equivalent H*(10).
    particle : {'neutron', 'photon', 'photon kerma', 'electron', 'positron'}
        Incident particle

    """
    key = (data_source, dose_quantity, particle)
    path = Path(__file__).parent / _TABLES[key][0]
    data = np.loadtxt(path, skiprows=3, encoding='utf-8')
    data[:, 0] *= 1e6   # Change energies to eV
    _DOSE_TABLES[key] = data


def dose_coefficients(
    particle, geometry='AP', data_source='icrp116', dose_quantity='effective'
):
    """Return dose conversion coefficients.

    This function provides fluence (and air kerma) to effective dose or ambient
    dose equivalent (H*(10)) conversion coefficients for various types of
    external exposures based on values in ICRP publications. Corrected values
    found in a corrigendum are used rather than the values in the original
    report. Available libraries include `ICRP Publication 74
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
        Irradiation geometry assumed for effective dose coefficients. Refer to
        ICRP-116 (Section 3.2) for the meaning of the options here. This
        argument does not apply when ``dose_quantity`` is 'ambient'.
    data_source : {'icrp74', 'icrp116'}
        The data source for the dose conversion coefficients.
    dose_quantity : {'effective', 'ambient'}
        Dose quantity to return. 'effective' returns effective dose
        coefficients; 'ambient' returns ambient dose equivalent (H*(10))
        coefficients.

    Returns
    -------
    energy : numpy.ndarray
        Energies at which dose conversion coefficients are given
    dose_coeffs : numpy.ndarray
        Dose coefficients in [pSv cm^2] at provided energies. For 'photon
        kerma', the coefficients are given in [Sv/Gy].

    """

    cv.check_value('geometry', geometry, _FULL_GEOMETRIES)
    cv.check_value('data_source', data_source, {'icrp74', 'icrp116'})
    cv.check_value('dose_quantity', dose_quantity, {'effective', 'ambient'})

    key = (data_source, dose_quantity, particle)
    if key not in _TABLES:
        available_particles = sorted(
            p for ds, dq, p in _TABLES
            if ds == data_source and dq == dose_quantity
        )
        msg = (
            f"'{particle}' has no {dose_quantity} dose data in data source "
            f"{data_source}. Available particles for {data_source} "
            f"with dose quantity {dose_quantity} are: {available_particles}"
        )
        raise ValueError(msg)
    elif key not in _DOSE_TABLES:
        _load_dose_table(data_source, dose_quantity, particle)

    data = _DOSE_TABLES[key]
    columns = _TABLES[key][1]

    if columns is None:
        if geometry != 'AP':
            raise ValueError(
                "Irradiation geometry is not defined for ambient dose "
                "equivalent coefficients. Use the default geometry='AP'."
            )
        index = 0
    else:
        index = columns.index(geometry)

    energy = data[:, 0].copy()
    dose_coeffs = data[:, index + 1].copy()
    return energy, dose_coeffs
