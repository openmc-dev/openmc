#
# import numpy as np
#
# from .function import Sum
# from .library import DataLibrary   # if you need it explicitly
# from . import K_BOLTZMANN
# from .photon import IncidentPhoton
#
#
# def linear_attenuation_xs(nuclide:str, temperature:float) -> Sum | None:
#     """Return a summed photon interaction cross section for a nuclide.
#
#     Parameters
#     ----------
#     nuclide : str
#         Name of nuclide.
#     temperature : float
#         Temperature in Kelvin.
#
#     Returns
#     -------
#     openmc.data.Sum or None
#         Sum of the relevant photon reaction cross sections as a function of
#         photon energy, or None if no photon data exist for nuclide.
#     """
#     strT = f"{int(round(temperature))}K"
#     photon_mts = {502, 504, 515, 516, 522}
#
#     # Load cross section library (uses OPENMC_CROSS_SECTIONS / config)
#     library = DataLibrary.from_xml()
#
#     lib = library.get_by_material(nuclide, data_type="photon")
#     if lib is None:
#         # No photon data for this nuclide; skip it
#         return None
#
#     # Load incident photon data
#     photon_data = IncidentPhoton.from_hdf5(lib["path"])
#
#     xs_list = []
#     # Sum the desired reaction channels to obtain a "total" photon xs
#     for reaction in photon_data.reactions.values():
#         mt = getattr(reaction, "mt", None)
#         if mt not in photon_mts:
#             continue
#
#         xs_obj = reaction.xs
#
#         # resolve xs for the temperature
#         if isinstance(xs_obj, dict):
#             # Try exact temperature match first
#             if strT in xs_obj:
#                 xs_T = xs_obj[strT]
#             else:
#                 # Fall back to nearest temperature if kTs/temperatures exist
#                 xs_T = None
#                 kTs = getattr(photon_data, "kTs", None)
#                 temps = getattr(photon_data, "temperatures", None)
#                 if kTs is not None and temps is not None and len(kTs) == len(temps):
#                     delta_T = np.array(kTs) - temperature * K_BOLTZMANN
#                     idx = int(np.argmin(np.abs(delta_T)))
#                     xs_T = xs_obj[temps[idx]]
#                 # If we still don't have a match, just take the first
#                 # available dataset as a last resort.
#                 if xs_T is None:
#                     xs_T = next(iter(xs_obj.values()))
#
#             xs = xs_T
#         else:
#             xs = xs_obj
#
#         xs_list.append(xs)
#
#     if len(xs_list) == 0:
#         return None
#     else:
#         return Sum(xs_list)
#
#
#
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
