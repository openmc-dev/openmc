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



def material_photon_mass_attenuation_dist(material) -> Sum | None:
    """Return material photon mass attenuation coefficient μ/ρ(E) [cm^2/g].

    the linear attenuation coefficient of the material is given by:
        μ(E) = Σ_el N_el * σ_el(E)
    with N_el in [atom/b-cm] and σ_el(E) in [barn/atom] => μ in [1/cm].

    The mass attenuation coefficients are given by:
        μ/ρ(E) = μ(E) / ρ
    => [1/cm] / [g/cm^3] = [cm^2/g]

    Parameters
    ----------
    material : openmc.Material

    Returns
    -------
    openmc.data.Sum or None
        Sum of Tabulated1D terms giving μ/ρ(E) in [cm^2/g], or None if no photon
        data exist for any constituents.
    """
    el_dens = material.get_element_atom_densities()
    if not el_dens:
        raise ValueError(
            f'For Material ID="{material.id}" no element densities are defined.'
        )

    # Mass density of the material [g/cm^3]
    rho = material.get_mass_density()  # g/cm^3

    if rho is None or rho <= 0.0:
        raise ValueError(
            f'Material ID="{material.id}" has non-positive mass density; '
            "cannot compute mass attenuation coefficient."
        )


    inv_rho = 1.0 / rho
    terms = []

    for el, n_el in el_dens.items():
        xs_sum = linear_attenuation_xs(el)  # barns/atom functions vs E
        if xs_sum is None or n_el == 0.0:
            continue

        scale = float(n_el) * inv_rho  # (atom/b-cm) / (g/cm^3) = (atom*cm^2)/(barn*g)

        for f in xs_sum.functions:
            if not isinstance(f, Tabulated1D):
                raise TypeError(
                    f"Expected Tabulated1D photon XS for element {el}, got {type(f)!r}."
                )
            # keep x, breakpoints, interpolation; scale y. 
            terms.append(
                Tabulated1D(
                    f.x,
                    np.asarray(f.y, dtype=float) * scale,
                    breakpoints=f.breakpoints,
                    interpolation=f.interpolation,
                )
            )

    return Sum(terms) if terms else None

