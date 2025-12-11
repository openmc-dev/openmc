import numpy as np

from .function import Sum
from .library import DataLibrary   
from .photon import IncidentPhoton
from openmc.exceptions import DataError

_PHOTON_LIB: DataLibrary | None = None
_PHOTON_DATA: dict[str, IncidentPhoton] = {}


def _get_photon_data(nuclide: str) ->IncidentPhoton | None:
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


def linear_attenuation_xs(nuclide: str, temperature: float) -> Sum | None:
    """Return total photon interaction cross section for a nuclide.

    Parameters
    ----------
    nuclide : str
        Name of nuclide.
    temperature : float
        Temperature in Kelvin.

    Returns
    -------
    openmc.data.Sum or None
        Sum of the relevant photon reaction cross sections as a function of
        photon energy, or None if no photon data exist for *nuclide*.
    """
    photon_data = _get_photon_data(nuclide)
    if photon_data is None:
        return None

    temp_key = f"{int(round(temperature))}K"
    photon_mts = (502, 504, 515, 517, 522)

    xs_list = []
    for reaction in photon_data.reactions.values():
        mt = getattr(reaction, "mt", None)
        if mt not in photon_mts:
            continue

        xs_obj = reaction.xs
        if isinstance(xs_obj, dict):
            if temp_key in xs_obj:
                xs_T = xs_obj[temp_key]
            else:
                # Fall back to closest available temperature
                temps = np.array(
                    [float(t.rstrip("K")) for t in xs_obj.keys()]
                )
                idx = int(np.argmin(np.abs(temps - temperature)))
                sel_key = f"{int(round(temps[idx]))}K"
                xs_T = xs_obj[sel_key]
            xs_list.append(xs_T)
        else:
            xs_list.append(xs_obj)

    if not xs_list:
        return None

    return Sum(xs_list)
