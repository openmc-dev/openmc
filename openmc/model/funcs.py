from collections.abc import Iterable
from math import sqrt
from operator import attrgetter
from warnings import warn

from openmc import Cylinder, Universe, Cell
from .surface_composite import RectangularPrism, HexagonalPrism
from ..checkvalue import (check_type, check_value, check_length,
                          check_less_than, check_iterable_type)
import openmc.data


ZERO_CELSIUS_TO_KELVIN = 273.15
ZERO_FAHRENHEIT_TO_KELVIN = 459.67
PSI_TO_MPA = 0.006895


def borated_water(boron_ppm, temperature=293., pressure=0.1013, temp_unit='K',
                  press_unit='MPa', density=None, **kwargs):
    """Return a Material with the composition of boron dissolved in water.

    The water density can be determined from a temperature and pressure, or it
    can be set directly.

    The concentration of boron has no effect on the stoichiometric ratio of H
    and O---they are fixed at 2-1.

    Parameters
    ----------
    boron_ppm : float
        The weight fraction in parts-per-million of elemental boron in the
        water.
    temperature : float
        Temperature in [K] used to compute water density.
    pressure : float
        Pressure in [MPa] used to compute water density.
    temp_unit : {'K', 'C', 'F'}
        The units used for the `temperature` argument.
    press_unit : {'MPa', 'psi'}
        The units used for the `pressure` argument.
    density : float
        Water density in [g / cm^3].  If specified, this value overrides 
        the value that is computed from the temperature and pressure arguments.
    **kwargs
        All keyword arguments are passed to the created Material object.

    Returns
    -------
    openmc.Material

    """
    # Perform any necessary unit conversions.
    check_value('temperature unit', temp_unit, ('K', 'C', 'F'))
    if temp_unit == 'K':
        T = temperature
    elif temp_unit == 'C':
        T = temperature + ZERO_CELSIUS_TO_KELVIN
    elif temp_unit == 'F':
        T = (temperature + ZERO_FAHRENHEIT_TO_KELVIN) * (5/9)
    check_value('pressure unit', press_unit, ('MPa', 'psi'))
    if press_unit == 'MPa':
        P = pressure
    elif press_unit == 'psi':
        P = pressure * PSI_TO_MPA

    # Set the density of water, either from an explicitly given density or from
    # temperature and pressure.
    if density is not None:
        water_density = density
    else:
        water_density = openmc.data.water_density(T, P)

    # Compute the density of the solution.
    solution_density = water_density / (1 - boron_ppm * 1e-6)

    # Compute the molar mass of pure water.
    hydrogen = openmc.Element('H')
    oxygen = openmc.Element('O')
    M_H2O = 0.0
    for iso_name, frac, junk in hydrogen.expand(2.0, 'ao'):
        M_H2O += frac * openmc.data.atomic_mass(iso_name)
    for iso_name, frac, junk in oxygen.expand(1.0, 'ao'):
        M_H2O += frac * openmc.data.atomic_mass(iso_name)

    # Compute the molar mass of boron.
    boron = openmc.Element('B')
    M_B = 0.0
    for iso_name, frac, junk in boron.expand(1.0, 'ao'):
        M_B += frac * openmc.data.atomic_mass(iso_name)

    # Compute the number fractions of each element.
    frac_H2O = (1 - boron_ppm * 1e-6) / M_H2O
    frac_H = 2 * frac_H2O
    frac_O = frac_H2O
    frac_B = boron_ppm * 1e-6 / M_B

    # Build the material.
    out = openmc.Material(temperature=T, **kwargs)
    out.add_element('H', frac_H, 'ao')
    out.add_element('O', frac_O, 'ao')
    out.add_element('B', frac_B, 'ao')
    out.set_density('g/cc', solution_density)
    out.add_s_alpha_beta('c_H_in_H2O')
    return out




def rectangular_prism(width, height, axis='z', origin=(0., 0.),
                      boundary_type='transmission', corner_radius=0.):
    warn("The rectangular_prism(...) function has been replaced by the "
         "RectangularPrism(...) class. Future versions of OpenMC will not "
         "accept rectangular_prism.", FutureWarning)
    return -RectangularPrism(
        width=width, height=height, axis=axis, origin=origin,
        boundary_type=boundary_type, corner_radius=corner_radius)


def hexagonal_prism(edge_length=1., orientation='y', origin=(0., 0.),
                    boundary_type='transmission', corner_radius=0.):
    warn("The hexagonal_prism(...) function has been replaced by the "
         "HexagonalPrism(...) class. Future versions of OpenMC will not "
         "accept hexagonal_prism.", FutureWarning)
    return -HexagonalPrism(
        edge_length=edge_length, orientation=orientation, origin=origin,
        boundary_type=boundary_type, corner_radius=corner_radius)


def get_hexagonal_prism(*args, **kwargs):
    warn("get_hexagonal_prism(...) has been renamed hexagonal_prism(...). "
         "Future versions of OpenMC will not accept get_hexagonal_prism.",
         FutureWarning)
    return hexagonal_prism(*args, **kwargs)


cylinder_from_points = Cylinder.from_points


def subdivide(surfaces):
    """Create regions separated by a series of surfaces.

    This function allows regions to be constructed from a set of a surfaces that
    are "in order". For example, if you had four instances of
    :class:`openmc.ZPlane` at z=-10, z=-5, z=5, and z=10, this function would
    return a list of regions corresponding to z < -10, -10 < z < -5, -5 < z < 5,
    5 < z < 10, and 10 < z. That is, for n surfaces, n+1 regions are returned.

    Parameters
    ----------
    surfaces : sequence of openmc.Surface
        Surfaces separating regions

    Returns
    -------
    list of openmc.Region
        Regions formed by the given surfaces

    """
    regions = [-surfaces[0]]
    for s0, s1 in zip(surfaces[:-1], surfaces[1:]):
        regions.append(+s0 & -s1)
    regions.append(+surfaces[-1])
    return regions


def pin(surfaces, items, subdivisions=None, divide_vols=True,
        **kwargs):
    """Convenience function for building a fuel pin

    Parameters
    ----------
    surfaces : iterable of :class:`openmc.Cylinder`
        Cylinders used to define boundaries
        between items. All cylinders must be
        concentric and of the same orientation, e.g.
        all :class:`openmc.ZCylinder`
    items : iterable
        Objects to go between ``surfaces``. These can be anything
        that can fill a :class:`openmc.Cell`, including
        :class:`openmc.Material`, or other :class:`openmc.Universe`
        objects. There must be one more item than surfaces,
        which will span all space outside the final ring.
    subdivisions : None or dict of int to int
        Dictionary describing which rings to subdivide and how
        many times. Keys are indexes of the annular rings
        to be divided. Will construct equal area rings
    divide_vols : bool
        If this evaluates to ``True``, then volumes of subdivided
        :class:`openmc.Material` instances will also be divided by the
        number of divisions.  Otherwise the volume of the
        original material will not be modified before subdivision
    kwargs:
        Additional key-word arguments to be passed to
        :class:`openmc.Universe`, like ``name="Fuel pin"``

    Returns
    -------
    :class:`openmc.Universe`
        Universe of concentric cylinders filled with the desired
        items
    """
    if "cells" in kwargs:
        raise ValueError(
            "Cells will be set by this function, not from input arguments.")
    check_type("items",  items, Iterable)
    check_length("surfaces", surfaces, len(items) - 1, len(items) - 1)
    # Check that all surfaces are of similar orientation
    check_type("surface", surfaces[0], Cylinder)
    surf_type = type(surfaces[0])
    check_iterable_type("surfaces", surfaces[1:], surf_type)

    # Check for increasing radii and equal centers
    if surf_type is openmc.ZCylinder:
        center_getter = attrgetter("x0", "y0")
    elif surf_type is openmc.YCylinder:
        center_getter = attrgetter("x0", "z0")
    elif surf_type is openmc.XCylinder:
        center_getter = attrgetter("z0", "y0")
    else:
        raise TypeError(
            f"Not configured to interpret {surf_type.__name__} surfaces")

    centers = set()
    prev_rad = 0
    for ix, surf in enumerate(surfaces):
        cur_rad = surf.r
        if cur_rad <= prev_rad:
            raise ValueError(
                "Surfaces do not appear to be increasing in radius. "
                "Surface {} at index {} has radius {:7.3e} compared to "
                "previous radius of {:7.5e}".format(
                    surf.id, ix, cur_rad, prev_rad))
        prev_rad = cur_rad
        centers.add(center_getter(surf))

    if len(centers) > 1:
        raise ValueError(
            "Surfaces do not appear to be concentric. The following "
            "centers were found: {}".format(centers))

    if subdivisions is not None:
        check_length("subdivisions", subdivisions, 1, len(surfaces))
        orig_indexes = list(subdivisions.keys())
        check_iterable_type("ring indexes", orig_indexes, int)
        check_iterable_type(
            "number of divisions", list(subdivisions.values()), int)
        for ix in orig_indexes:
            if ix < 0:
                subdivisions[len(surfaces) + ix] = subdivisions.pop(ix)
        # Dissallow subdivision on outer most, infinite region
        check_less_than(
            "outer ring", max(subdivisions), len(surfaces), equality=True)

        # ensure ability to concatenate
        if not isinstance(items, list):
            items = list(items)
        if not isinstance(surfaces, list):
            surfaces = list(surfaces)

        # generate equal area divisions
        # Adding N - 1 new regions
        # N - 2 surfaces are made
        # Original cell is not removed, but not occupies last ring
        for ring_index in reversed(sorted(subdivisions.keys())):
            nr = subdivisions[ring_index]
            new_surfs = []

            lower_rad = 0.0 if ring_index == 0 else surfaces[ring_index - 1].r

            upper_rad = surfaces[ring_index].r

            area_term = (upper_rad ** 2 - lower_rad ** 2) / nr

            for new_index in range(nr - 1):
                lower_rad = sqrt(area_term + lower_rad ** 2)
                new_surfs.append(surf_type(r=lower_rad))

            surfaces = (
                    surfaces[:ring_index] + new_surfs + surfaces[ring_index:])

            filler = items[ring_index]
            if (divide_vols and hasattr(filler, "volume")
                    and filler.volume is not None):
                filler.volume /= nr

            items[ring_index:ring_index] = [
                filler.clone() for _i in range(nr - 1)]

    # Build the universe
    regions = subdivide(surfaces)
    cells = [Cell(fill=f, region=r) for r, f in zip(regions, items)]
    return Universe(cells=cells, **kwargs)
