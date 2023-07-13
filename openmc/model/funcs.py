from collections.abc import Iterable
from functools import partial
from math import sqrt, pi, sin, cos
from numbers import Real
from operator import attrgetter
from warnings import warn

from openmc import (
    XPlane, YPlane, Plane, ZCylinder, Cylinder, XCylinder,
    YCylinder, Universe, Cell)
from ..checkvalue import (
    check_type, check_value, check_length, check_less_than,
    check_iterable_type)
import openmc.data
from openmc.model.surface_composite import (CompositeSurface,CylinderSector)

import numpy as np


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
        Water density in [g / cm^3].  If specified, this value overrides the
        temperature and pressure arguments.
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
    if density is None:
        out = openmc.Material(temperature=T, **kwargs)
    else:
        out = openmc.Material(**kwargs)
    out.add_element('H', frac_H, 'ao')
    out.add_element('O', frac_O, 'ao')
    out.add_element('B', frac_B, 'ao')
    out.set_density('g/cc', solution_density)
    out.add_s_alpha_beta('c_H_in_H2O')
    return out


def rectangular_prism(width, height, axis='z', origin=(0., 0.),
                      boundary_type='transmission', corner_radius=0.):
    """Get an infinite rectangular prism from four planar surfaces.

    .. versionchanged:: 0.11
        This function was renamed from `get_rectangular_prism` to
        `rectangular_prism`.

    Parameters
    ----------
    width: float
        Prism width in units of cm. The width is aligned with the y, x,
        or x axes for prisms parallel to the x, y, or z axis, respectively.
    height: float
        Prism height in units of cm. The height is aligned with the z, z,
        or y axes for prisms parallel to the x, y, or z axis, respectively.
    axis : {'x', 'y', 'z'}
        Axis with which the infinite length of the prism should be aligned.
        Defaults to 'z'.
    origin: Iterable of two floats
        Origin of the prism. The two floats correspond to (y,z), (x,z) or
        (x,y) for prisms parallel to the x, y or z axis, respectively.
        Defaults to (0., 0.).
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic'}
        Boundary condition that defines the behavior for particles hitting the
        surfaces comprising the rectangular prism (default is 'transmission').
    corner_radius: float
        Prism corner radius in units of cm. Defaults to 0.

    Returns
    -------
    openmc.Region
        The inside of a rectangular prism

    """

    check_type('width', width, Real)
    check_type('height', height, Real)
    check_type('corner_radius', corner_radius, Real)
    check_value('axis', axis, ['x', 'y', 'z'])
    check_type('origin', origin, Iterable, Real)

    # Define function to create a plane on given axis
    def plane(axis, name, value):
        cls = globals()['{}Plane'.format(axis.upper())]
        return cls(name='{} {}'.format(name, axis),
                   boundary_type=boundary_type,
                   **{axis + '0': value})

    if axis == 'x':
        x1, x2 = 'y', 'z'
    elif axis == 'y':
        x1, x2 = 'x', 'z'
    else:
        x1, x2 = 'x', 'y'

    # Get cylinder class corresponding to given axis
    cyl = globals()['{}Cylinder'.format(axis.upper())]

    # Create rectangular region
    min_x1 = plane(x1, 'minimum', -width/2 + origin[0])
    max_x1 = plane(x1, 'maximum', width/2 + origin[0])
    min_x2 = plane(x2, 'minimum', -height/2 + origin[1])
    max_x2 = plane(x2, 'maximum', height/2 + origin[1])
    if boundary_type == 'periodic':
        min_x1.periodic_surface = max_x1
        min_x2.periodic_surface = max_x2
    prism = +min_x1 & -max_x1 & +min_x2 & -max_x2

    # Handle rounded corners if given
    if corner_radius > 0.:
        if boundary_type == 'periodic':
            raise ValueError('Periodic boundary conditions not permitted when '
                             'rounded corners are used.')

        args = {'r': corner_radius, 'boundary_type': boundary_type}

        args[x1 + '0'] = origin[0] - width/2 + corner_radius
        args[x2 + '0'] = origin[1] - height/2 + corner_radius
        x1_min_x2_min = cyl(name='{} min {} min'.format(x1, x2), **args)

        args[x1 + '0'] = origin[0] - width/2 + corner_radius
        args[x2 + '0'] = origin[1] - height/2 + corner_radius
        x1_min_x2_min = cyl(name='{} min {} min'.format(x1, x2), **args)

        args[x1 + '0'] = origin[0] - width/2 + corner_radius
        args[x2 + '0'] = origin[1] + height/2 - corner_radius
        x1_min_x2_max = cyl(name='{} min {} max'.format(x1, x2), **args)

        args[x1 + '0'] = origin[0] + width/2 - corner_radius
        args[x2 + '0'] = origin[1] - height/2 + corner_radius
        x1_max_x2_min = cyl(name='{} max {} min'.format(x1, x2), **args)

        args[x1 + '0'] = origin[0] + width/2 - corner_radius
        args[x2 + '0'] = origin[1] + height/2 - corner_radius
        x1_max_x2_max = cyl(name='{} max {} max'.format(x1, x2), **args)

        x1_min = plane(x1, 'min', -width/2 + origin[0] + corner_radius)
        x1_max = plane(x1, 'max', width/2 + origin[0] - corner_radius)
        x2_min = plane(x2, 'min', -height/2 + origin[1] + corner_radius)
        x2_max = plane(x2, 'max', height/2 + origin[1] - corner_radius)

        corners = (+x1_min_x2_min & -x1_min & -x2_min) | \
                  (+x1_min_x2_max & -x1_min & +x2_max) | \
                  (+x1_max_x2_min & +x1_max & -x2_min) | \
                  (+x1_max_x2_max & +x1_max & +x2_max)

        prism = prism & ~corners

    return prism


def get_rectangular_prism(*args, **kwargs):
    warn("get_rectangular_prism(...) has been renamed rectangular_prism(...). "
         "Future versions of OpenMC will not accept get_rectangular_prism.",
         FutureWarning)
    return rectangular_prism(*args, **kwargs)


def hexagonal_prism(edge_length=1., orientation='y', origin=(0., 0.),
                    boundary_type='transmission', corner_radius=0.):
    """Create a hexagon region from six surface planes.

    .. versionchanged:: 0.11
        This function was renamed from `get_hexagonal_prism` to
        `hexagonal_prism`.

    Parameters
    ----------
    edge_length : float
        Length of a side of the hexagon in cm
    orientation : {'x', 'y'}
        An 'x' orientation means that two sides of the hexagon are parallel to
        the x-axis and a 'y' orientation means that two sides of the hexagon are
        parallel to the y-axis.
    origin: Iterable of two floats
        Origin of the prism. Defaults to (0., 0.).
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic'}
        Boundary condition that defines the behavior for particles hitting the
        surfaces comprising the hexagonal prism (default is 'transmission').
    corner_radius: float
        Prism corner radius in units of cm. Defaults to 0.

    Returns
    -------
    openmc.Region
        The inside of a hexagonal prism

    """

    l = edge_length
    x, y = origin

    if orientation == 'y':
        right = XPlane(x + sqrt(3.)/2*l, boundary_type=boundary_type)
        left = XPlane(x - sqrt(3.)/2*l, boundary_type=boundary_type)
        c = sqrt(3.)/3.

        # y = -x/sqrt(3) + a
        upper_right = Plane(a=c, b=1., d=l+x*c+y, boundary_type=boundary_type)

        # y = x/sqrt(3) + a
        upper_left = Plane(a=-c, b=1., d=l-x*c+y, boundary_type=boundary_type)

        # y = x/sqrt(3) - a
        lower_right = Plane(a=-c, b=1., d=-l-x*c+y, boundary_type=boundary_type)

        # y = -x/sqrt(3) - a
        lower_left = Plane(a=c, b=1., d=-l+x*c+y, boundary_type=boundary_type)

        prism = -right & +left & -upper_right & -upper_left & \
                +lower_right & +lower_left

        if boundary_type == 'periodic':
            right.periodic_surface = left
            upper_right.periodic_surface = lower_left
            lower_right.periodic_surface = upper_left

    elif orientation == 'x':
        top = YPlane(y0=y + sqrt(3.)/2*l, boundary_type=boundary_type)
        bottom = YPlane(y0=y - sqrt(3.)/2*l, boundary_type=boundary_type)
        c = sqrt(3.)

        # y = -sqrt(3)*(x - a)
        upper_right = Plane(a=c, b=1., d=c*l+x*c+y, boundary_type=boundary_type)

        # y = sqrt(3)*(x + a)
        lower_right = Plane(a=-c, b=1., d=-c*l-x*c+y,
                            boundary_type=boundary_type)

        # y = -sqrt(3)*(x + a)
        lower_left = Plane(a=c, b=1., d=-c*l+x*c+y, boundary_type=boundary_type)

        # y = sqrt(3)*(x + a)
        upper_left = Plane(a=-c, b=1., d=c*l-x*c+y, boundary_type=boundary_type)

        prism = -top & +bottom & -upper_right & +lower_right & \
                            +lower_left & -upper_left

        if boundary_type == 'periodic':
            top.periodic_surface = bottom
            upper_right.periodic_surface = lower_left
            lower_right.periodic_surface = upper_left

    # Handle rounded corners if given
    if corner_radius > 0.:
        if boundary_type == 'periodic':
            raise ValueError('Periodic boundary conditions not permitted when '
                             'rounded corners are used.')

        c = sqrt(3.)/2
        t = l - corner_radius/c

        # Cylinder with corner radius and boundary type pre-applied
        cyl1 = partial(ZCylinder, r=corner_radius, boundary_type=boundary_type)
        cyl2 = partial(ZCylinder, r=corner_radius/(2*c),
                       boundary_type=boundary_type)

        if orientation == 'x':
            x_min_y_min_in = cyl1(name='x min y min in', x0=x-t/2, y0=y-c*t)
            x_min_y_max_in = cyl1(name='x min y max in', x0=x+t/2, y0=y-c*t)
            x_max_y_min_in = cyl1(name='x max y min in', x0=x-t/2, y0=y+c*t)
            x_max_y_max_in = cyl1(name='x max y max in', x0=x+t/2, y0=y+c*t)
            x_min_in = cyl1(name='x min in', x0=x-t, y0=y)
            x_max_in = cyl1(name='x max in', x0=x+t, y0=y)

            x_min_y_min_out = cyl2(name='x min y min out', x0=x-l/2, y0=y-c*l)
            x_min_y_max_out = cyl2(name='x min y max out', x0=x+l/2, y0=y-c*l)
            x_max_y_min_out = cyl2(name='x max y min out', x0=x-l/2, y0=y+c*l)
            x_max_y_max_out = cyl2(name='x max y max out', x0=x+l/2, y0=y+c*l)
            x_min_out = cyl2(name='x min out', x0=x-l, y0=y)
            x_max_out = cyl2(name='x max out', x0=x+l, y0=y)

            corners = (+x_min_y_min_in & -x_min_y_min_out |
                       +x_min_y_max_in & -x_min_y_max_out |
                       +x_max_y_min_in & -x_max_y_min_out |
                       +x_max_y_max_in & -x_max_y_max_out |
                       +x_min_in & -x_min_out |
                       +x_max_in & -x_max_out)

        elif orientation == 'y':
            x_min_y_min_in = cyl1(name='x min y min in', x0=x-c*t, y0=y-t/2)
            x_min_y_max_in = cyl1(name='x min y max in', x0=x-c*t, y0=y+t/2)
            x_max_y_min_in = cyl1(name='x max y min in', x0=x+c*t, y0=y-t/2)
            x_max_y_max_in = cyl1(name='x max y max in', x0=x+c*t, y0=y+t/2)
            y_min_in = cyl1(name='y min in', x0=x, y0=y-t)
            y_max_in = cyl1(name='y max in', x0=x, y0=y+t)

            x_min_y_min_out = cyl2(name='x min y min out', x0=x-c*l, y0=y-l/2)
            x_min_y_max_out = cyl2(name='x min y max out', x0=x-c*l, y0=y+l/2)
            x_max_y_min_out = cyl2(name='x max y min out', x0=x+c*l, y0=y-l/2)
            x_max_y_max_out = cyl2(name='x max y max out', x0=x+c*l, y0=y+l/2)
            y_min_out = cyl2(name='y min out', x0=x, y0=y-l)
            y_max_out = cyl2(name='y max out', x0=x, y0=y+l)

            corners = (+x_min_y_min_in & -x_min_y_min_out |
                       +x_min_y_max_in & -x_min_y_max_out |
                       +x_max_y_min_in & -x_max_y_min_out |
                       +x_max_y_max_in & -x_max_y_max_out |
                       +y_min_in & -y_min_out |
                       +y_max_in & -y_max_out)

        prism = prism & ~corners

    return prism


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
    
    # Every surface in the list "surfaces" is either a Cylinder or a CylinderSector
    # The region inside each surface is added to the list "regions" depending on what the current and next surface is
    
    regions = []
    
    # If the first surface is a Cylinder, the region inside the cylinder is added to "regions".
    if isinstance(surfaces[0], Cylinder):
        regions.append(-surfaces[0])
       
    # This for loop iterates through "surfaces", adding to regions
    # If the surface is a CylinderSector, the region inside the surface is added to "regions".
    # If the current surface and the next surface is a Cylinder, then the region between the cylinders is added to "regions".
    for s0, s1 in zip(surfaces[:-1], surfaces[1:]):
        if isinstance(s0, CylinderSector):
            regions.append(-s0)
        
        elif isinstance(s0, Cylinder) and isinstance(s1, Cylinder):
            regions.append(+s0 & -s1)
    
    # Add the region outside the last surface.
    if isinstance(surfaces[-1], Cylinder):
        regions.append(+surfaces[-1])
    else:
        regions.append(-surfaces[-1])
    
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
    if surf_type is ZCylinder:
        center_getter = attrgetter("x0", "y0")
    elif surf_type is YCylinder:
        center_getter = attrgetter("x0", "z0")
    elif surf_type is XCylinder:
        center_getter = attrgetter("z0", "y0")
    else:
        raise TypeError(
            "Not configured to interpret {} surfaces".format(
                surf_type.__name__))

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



def pin_radial_azimuthal(surfaces, items, subdivisions_r=None, subdivisions_a=None, rad_div_types=None, implicit_azi_div=None,
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
    subdivisions_r : None or dict of int to int
        Dictionary describing which rings to subdivide radially and how
        many times. Keys are indexes of the annular rings
        to be divided. Will construct equal area or equal radius rings
        depending on value of rad_div_type
    subdivisions_a : None or dict of int to int
        Dictionary describing which rings to subdivide azimuthally and how
        many times. Keys are indexes of the annular rings
        to be divided. Will construct equal area sectors
    rad_div_types : None or dict of int to string
        Dictionary descibing how to subdivide rings radially, either with 
        equal area or equal radius. Keys are indexes of the annular rings
        to be divided. Values should be either "area" or "radius".
        A value of "area" will create equal area rings while 
        a value of "radius" will create equal radius rings.
        The division type will default to equal area.
    implicit_azi_div : None or int
        Value describes how to azimuthally divide the implicit outer region.
        Defaults to None, meaning there are no divisions.
    kwargs:
        Additional key-word arguments to be passed to
        :class:`openmc.Universe`, like ``name="Fuel pin"``

    Returns
    -------
    :class:`openmc.Universe`
        Universe of concentric cylinders and CylinderSectors filled with the desired
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
    if surf_type is ZCylinder:
        center_getter = attrgetter("x0", "y0")
    elif surf_type is YCylinder:
        center_getter = attrgetter("x0", "z0")
    elif surf_type is XCylinder:
        center_getter = attrgetter("z0", "y0")
    else:
        raise TypeError(
            "Not configured to interpret {} surfaces".format(
                surf_type.__name__))

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
    
    items_new = items.copy()
    
    # Divides Cylinders into more rings
    if subdivisions_r is not None:
        check_length("subdivisions_r", subdivisions_r, 1, len(surfaces))
        orig_indexes = list(subdivisions_r.keys())
        check_iterable_type("ring indexes", orig_indexes, int)
        check_iterable_type(
            "number of divisions", list(subdivisions_r.values()), int)
        for ix in orig_indexes:
            if ix < 0:
                subdivisions_r[len(surfaces) + ix] = subdivisions_r.pop(ix)
        # Dissallow subdivision on outer most, infinite region
        check_less_than(
            "outer ring", max(subdivisions_r), len(surfaces), equality=True)

        # ensure ability to concatenate
        if not isinstance(items, list):
            items = list(items)
        if not isinstance(surfaces, list):
            surfaces = list(surfaces)

        # generate equal area divisions
        # Adding N - 1 new regions
        # N - 2 surfaces are made
        # Original cell is not removed, but not occupies last ring
        for ring_index in reversed(sorted(subdivisions_r.keys())):
            nr = subdivisions_r[ring_index]
            if nr < 1:
                raise ValueError("The number of subdivisions must be a positive integer.")
            new_surfs = []

            lower_rad = 0.0 if ring_index == 0 else surfaces[ring_index - 1].r

            upper_rad = surfaces[ring_index].r

            area_term = (upper_rad ** 2 - lower_rad ** 2) / nr

            equal_radius_term = (upper_rad - lower_rad) / nr

            if rad_div_types is None or ring_index not in rad_div_types.keys():
                div_type = "area"
            else:
                div_type = rad_div_types[ring_index]

            for new_index in range(nr - 1):
                # generate equal area divisions
                if div_type == "area":
                    lower_rad = sqrt(area_term + lower_rad ** 2)
                
                # generate equal radius divisions
                else:
                    lower_rad = lower_rad + equal_radius_term
                
                new_surfs.append(surf_type(r=lower_rad))

            surfaces = (
                    surfaces[:ring_index] + new_surfs + surfaces[ring_index:])

            filler = items_new[ring_index]
            items_new[ring_index:ring_index] = [
                filler for _i in range(nr - 1)]
   

    # Loop to correct "subdivisions_a" dictionary after creating more rings
    if subdivisions_r is not None and subdivisions_a is not None:
        subdivisions_a_fixed = {}
        counter = 0
        
        for i in range(len(items)):
            if i in subdivisions_r.keys():
                for j in range(subdivisions_r[i]):
                    if i in subdivisions_a.keys():
                        subdivisions_a_fixed[counter] = subdivisions_a[i]
                    counter += 1
            
            elif i in subdivisions_a.keys():
                subdivisions_a_fixed[counter] = subdivisions_a[i]
                counter += 1
            
            else:
                counter += 1

        subdivisions_a = subdivisions_a_fixed


    # Divides Cylinders into CylinderSectors
    if subdivisions_a is not None:
        check_length("subdivisions_a", subdivisions_a, 1, len(surfaces))
        orig_indexes = list(subdivisions_a.keys())
        check_iterable_type("ring indexes", orig_indexes, int)
        check_iterable_type(
            "number of divisions", list(subdivisions_a.values()), int)
        for ix in orig_indexes:
            if ix < 0:
                subdivisions_a[len(surfaces) + ix] = subdivisions_a.pop(ix)
        # Dissallow subdivision on outer most, infinite region
        check_less_than(
            "outer ring", max(subdivisions_a), len(surfaces), equality=True)

        # ensure ability to concatenate
        if not isinstance(items, list):
            items = list(items)
        if not isinstance(surfaces, list):
            surfaces = list(surfaces)

        # Generate azimuthal divisions
        for ring_index in reversed(sorted(subdivisions_a.keys())):
            ns = subdivisions_a[ring_index]
            if ns < 1:
                raise ValueError("The number of subdivisions must be a positive integer.")
            if ns == 1:
                continue
            new_surfs = []
            center = list(centers)[0]

            lower_rad = 0.0 if ring_index == 0 else surfaces[ring_index - 1].r

            upper_rad = surfaces[ring_index].r
            
            # Creates a CylinderSector (composite surface of two cocentric cylinders and two planes) for each azimuthal region
            spacing = 360.0/ns
            
            angles = [(i * spacing - 45.0) for i in range(ns)]
            angles.append(315.0)
            

            for ix in range(ns):
                lower_angle = angles[ix]
                upper_angle = angles[ix+1]

                if surf_type is ZCylinder:
                    cyl_sec = CylinderSector(r1=lower_rad, r2=upper_rad, theta1=lower_angle, theta2=upper_angle, 
                            center=center, axis='z')
                elif surf_type is YCylinder:
                    cyl_sec = CylinderSector(r1=lower_rad, r2=upper_rad, theta1=lower_angle, theta2=upper_angle, 
                            center=center, axis='y')
                elif surf_type is XCylinder:
                    cyl_sec = CylinderSector(r1=lower_rad, r2=upper_rad, theta1=lower_angle, theta2=upper_angle, 
                            center=center, axis='x')
                
                new_surfs.append(cyl_sec)


            surfaces = (
                    surfaces[:ring_index] + new_surfs + surfaces[ring_index:])

            filler = items_new[ring_index]
            items_new[ring_index:ring_index] = [
                filler for _i in range(ns - 1)]

    if implicit_azi_div is not None and implicit_azi_div != 1:
        ns = implicit_azi_div
        if ns < 1:
            raise ValueError("The number of subdivisions must be a positive integer.")
        new_surfs = []
        center = list(centers)[0]
            
        lower_rad = surfaces[-1].r
        upper_rad = None
            
        # Creates a CylinderSector (composite surface of two cocentric cylinders and two planes) for each azimuthal region
        spacing = 360.0/ns
            
        angles = [(i * spacing - 45.0) for i in range(ns)]
        angles.append(315.0)

        for ix in range(ns):
            lower_angle = angles[ix]
            upper_angle = angles[ix+1]

            if surf_type is ZCylinder:
                cyl_sec = CylinderSector(r1=lower_rad, r2=upper_rad, theta1=lower_angle, theta2=upper_angle, 
                        center=center, axis='z')
            elif surf_type is YCylinder:
                cyl_sec = CylinderSector(r1=lower_rad, r2=upper_rad, theta1=lower_angle, theta2=upper_angle, 
                        center=center, axis='y')
            elif surf_type is XCylinder:
                cyl_sec = CylinderSector(r1=lower_rad, r2=upper_rad, theta1=lower_angle, theta2=upper_angle, 
                        center=center, axis='x')
                
            new_surfs.append(cyl_sec)


        for _s in new_surfs:
            surfaces.append(_s)
            
        filler = items_new[-1]
        for _i in range(ns - 1):
            items_new.append(filler)



    # Build the universe
    regions = subdivide(surfaces)
    cells = [Cell(fill=f, region=r) for r, f in zip(regions, items_new)]
    return Universe(cells=cells, **kwargs)
