from collections import OrderedDict
from collections.abc import Iterable
from math import sqrt
from numbers import Real
from functools import partial

from openmc import XPlane, YPlane, Plane, ZCylinder, Quadric
from openmc.checkvalue import check_type, check_value
import openmc.data


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
        T = temperature + 273.15
    elif temp_unit == 'F':
        T = (temperature + 459.67) * 5.0 / 9.0
    check_value('pressure unit', press_unit, ('MPa', 'psi'))
    if press_unit == 'MPa':
        P = pressure
    elif press_unit == 'psi':
        P = pressure * 0.006895

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


def get_rectangular_prism(width, height, axis='z', origin=(0., 0.),
                          boundary_type='transmission', corner_radius=0.):
    """Get an infinite rectangular prism from four planar surfaces.

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
    prism = +min_x1 & -max_x1 & +min_x2 & -max_x2

    # Handle rounded corners if given
    if corner_radius > 0.:
        args = {'R': corner_radius, 'boundary_type': boundary_type}

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


def get_hexagonal_prism(edge_length=1., orientation='y', origin=(0., 0.),
                        boundary_type='transmission', corner_radius=0.):
    """Create a hexagon region from six surface planes.

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
        right = XPlane(x0=x + sqrt(3.)/2*l, boundary_type=boundary_type)
        left = XPlane(x0=x - sqrt(3.)/2*l, boundary_type=boundary_type)
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


def cylinder_from_points(p1, p2, r, **kwargs):
    """Return cylinder defined by two points passing through its center.

    Parameters
    ----------
    p1, p2 : 3-tuples
        Coordinates of two points that pass through the center of the cylinder
    r : float
        Radius of the cylinder
    kwargs : dict
        Keyword arguments passed to the :class:`openmc.Quadric` constructor

    Returns
    -------
    openmc.Quadric
        Quadric surface representing the cylinder.

    """
    # Get x, y, z coordinates of two points
    x1, y1, z1 = p1
    x2, y2, z2 = p2

    # Define intermediate terms
    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1
    cx = y1*z2 + y2*z1
    cy = -(x1*z2 + x2*z1)
    cz = x1*y2 + x2*y1

    # Given p=(x,y,z), p1=(x1, y1, z1), p2=(x2, y2, z2), the equation for the
    # cylinder can be derived as r = |(p - p1) ⨯ (p - p2)| / |p2 - p1|.
    # Expanding out all terms and grouping according to what Quadric expects
    # gives the following coefficients.
    kwargs['a'] = dy*dy + dz*dz
    kwargs['b'] = dx*dx + dz*dz
    kwargs['c'] = dx*dx + dy*dy
    kwargs['d'] = -2*dx*dy
    kwargs['e'] = -2*dy*dz
    kwargs['f'] = -2*dx*dz
    kwargs['g'] = cy*dz - cz*dy
    kwargs['h'] = cz*dx - cx*dz
    kwargs['j'] = cx*dy - cy*dx
    kwargs['k'] = -(dx*dx + dy*dy + dz*dz)*r*r

    return openmc.Quadric(**kwargs)


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
