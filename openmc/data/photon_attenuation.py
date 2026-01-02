import numpy as np

from openmc.exceptions import DataError

from .data import ATOMIC_SYMBOL, ELEMENT_SYMBOL, zam
from .function import Sum, Tabulated1D
from .library import DataLibrary
from .photon import IncidentPhoton


_PHOTON_LIB: DataLibrary | None = None
_PHOTON_DATA: dict[str, IncidentPhoton] = {}


def _get_photon_data(nuclide: str) -> IncidentPhoton | None:
    global _PHOTON_LIB

    if _PHOTON_LIB is None:
        try:
            _PHOTON_LIB = DataLibrary.from_xml()
        except Exception as err:
            raise DataError(
                "A cross section library must be specified with "
                "openmc.config['cross_sections'] in order to load photon data."
            ) from err

    lib = _PHOTON_LIB.get_by_material(nuclide, data_type="photon")
    if lib is None:
        return None

    if nuclide not in _PHOTON_DATA:
        _PHOTON_DATA[nuclide] = IncidentPhoton.from_hdf5(lib["path"])

    return _PHOTON_DATA[nuclide]


def linear_attenuation_xs(element_input: str) -> Sum | None:
    """Return total photon interaction cross section for a nuclide.

    Parameters
    ----------
    element_input : str
        Name of nuclide or element

    Returns
    -------
    openmc.data.Sum or None
        Sum of the relevant photon reaction cross sections as a function of
        photon energy, or None if no photon data exist for *nuclide*.
    """

    try:
        z = zam(element_input)[0]
        element = ATOMIC_SYMBOL[z]
    except (ValueError, KeyError, TypeError):
        if element_input not in ELEMENT_SYMBOL.values():
            raise ValueError(f"Element '{element_input}' not found in ELEMENT_SYMBOL.")
        element = element_input

    photon_data = _get_photon_data(element)
    if photon_data is None:
        return None

    photon_mts = (502, 504, 515, 517, 522)

    xs_list = []
    for reaction in photon_data.reactions.values():
        mt = getattr(reaction, "mt", None)
        if mt not in photon_mts:
            continue

        xs_list.append(reaction.xs)

    if not xs_list:
        return None

    return Sum(xs_list)



