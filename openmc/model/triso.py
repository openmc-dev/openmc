import copy
import warnings
import itertools
import random
from abc import ABCMeta, abstractproperty, abstractmethod
from collections import Counter, defaultdict
from collections.abc import Iterable
from heapq import heappush, heappop
from math import pi, sin, cos, floor, log10, sqrt
from numbers import Real
from random import uniform, gauss

import numpy as np
import scipy.spatial

import openmc
from openmc.checkvalue import check_type


MAX_PF_RSP = 0.38
MAX_PF_CRP = 0.64


class TRISO(openmc.Cell):
    """Tristructural-isotopic (TRISO) micro fuel particle

    Parameters
    ----------
    outer_radius : float
        Outer radius of TRISO particle
    fill : openmc.Universe
        Universe which contains all layers of the TRISO particle
    center : Iterable of float
        Cartesian coordinates of the center of the TRISO particle in cm

    Attributes
    ----------
    id : int
        Unique identifier for the TRISO cell
    name : str
        Name of the TRISO cell
    center : numpy.ndarray
        Cartesian coordinates of the center of the TRISO particle in cm
    fill : openmc.Universe
        Universe that contains the TRISO layers
    region : openmc.Region
        Region of space within the TRISO particle

    """

    def __init__(self, outer_radius, fill, center=(0., 0., 0.)):
        self._surface = openmc.Sphere(r=outer_radius)
        super().__init__(fill=fill, region=-self._surface)
        self.center = np.asarray(center)

    @property
    def center(self):
        return self._center

    @center.setter
    def center(self, center):
        check_type('TRISO center', center, Iterable, Real)
        self._surface.x0 = center[0]
        self._surface.y0 = center[1]
        self._surface.z0 = center[2]
        self.translation = center
        self._center = center

    def classify(self, lattice):
        """Determine lattice element indices which might contain the TRISO particle.

        Parameters
        ----------
        lattice : openmc.RectLattice
            Lattice to check

        Returns
        -------
        list of tuple
            (z,y,x) lattice element indices which might contain the TRISO
            particle.

        """

        ll, ur = self.region.bounding_box
        if lattice.ndim == 2:
            (i_min, j_min), p = lattice.find_element(ll)
            (i_max, j_max), p = lattice.find_element(ur)
            return list(np.broadcast(*np.ogrid[
                j_min:j_max+1, i_min:i_max+1]))
        else:
            (i_min, j_min, k_min), p = lattice.find_element(ll)
            (i_max, j_max, k_max), p = lattice.find_element(ur)
            return list(np.broadcast(*np.ogrid[
                k_min:k_max+1, j_min:j_max+1, i_min:i_max+1]))


class _Container(metaclass=ABCMeta):
    """Container in which to pack spheres.

    Parameters
    ----------
    sphere_radius : float
        Radius of spheres to be packed in container.
    center : Iterable of float
        Cartesian coordinates of the center of the container. Default is
        (0., 0., 0.)

    Attributes
    ----------
    sphere_radius : float
        Radius of spheres to be packed in container.
    center : list of float
        Cartesian coordinates of the center of the container. Default is
        (0., 0., 0.)
    cell_length : list of float
        Length in x-, y-, and z- directions of each cell in mesh overlaid on
        domain.
    limits : list of float
        Constraint on where sphere center can be placed.
    volume : float
        Volume of the container.

    """
    def __init__(self, sphere_radius, center=(0., 0., 0.)):
        self._cell_length = None
        self._limits = None

        self.sphere_radius = sphere_radius
        self.center = center

    @property
    def sphere_radius(self):
        return self._sphere_radius

    @property
    def center(self):
        return self._center

    @abstractproperty
    def limits(self):
        pass

    @abstractproperty
    def cell_length(self):
        pass

    @abstractproperty
    def volume(self):
        pass

    @sphere_radius.setter
    def sphere_radius(self, sphere_radius):
        self._sphere_radius = float(sphere_radius)
        self._limits = None
        self._cell_length = None

    @center.setter
    def center(self, center):
        self._center = center

    def mesh_cell(self, p):
        """Calculate the index of the cell in a mesh overlaid on the domain in
        which the given sphere center falls.

        Parameters
        ----------
        p : Iterable of float
            Cartesian coordinates of sphere center.

        Returns
        -------
        tuple of int
            Indices of mesh cell.

        """
        return tuple(int(p[i]/self.cell_length[i]) for i in range(3))

    def nearby_mesh_cells(self, p):
        """Calculates the indices of all cells in a mesh overlaid on the domain
        within one diameter of the given sphere.

        Parameters
        ----------
        p : Iterable of float
            Cartesian coordinates of sphere center.

        Returns
        -------
        list of tuple of int
            Indices of mesh cells.

        """
        d = 2*self.sphere_radius
        r = [[a/self.cell_length[i] for a in [p[i]-d, p[i], p[i]+d]]
             for i in range(3)]
        return list(itertools.product(*({int(x) for x in y} for y in r)))

    @abstractmethod
    def from_region(self, region, sphere_radius):
        """Create a container to pack spheres in based on a region.

        Parameters
        ----------
        region : openmc.Region
            Region to create container from.
        sphere_radius : float
            Outer radius of spheres.

        """
        pass

    @abstractmethod
    def random_point(self):
        """Generate Cartesian coordinates of center of a sphere that is
        contained entirely within the domain with uniform probability.

        Returns
        -------
        list of float
            Cartesian coordinates of sphere center.

        """
        pass

    @abstractmethod
    def repel_spheres(self, p, q, d, d_new):
        """Move spheres p and q apart according to the following
        transformation (accounting for boundary conditions on domain):

            r_i^(n+1) = r_i^(n) + 1/2(d_out^(n+1) - d^(n))
            r_j^(n+1) = r_j^(n) - 1/2(d_out^(n+1) - d^(n))

        Parameters
        ----------
        p, q : numpy.ndarray
            Cartesian coordinates of sphere center.
        d : float
            distance between centers of spheres i and j.
        d_new : float
            final distance between centers of spheres i and j.

        """
        pass


class _RectangularPrism(_Container):
    """Rectangular prism container in which to pack spheres.

    Parameters
    ----------
    width : float
        Prism length along the x-axis
    depth : float
        Prism length along the y-axis
    height : float
        Prism length along the z-axis
    sphere_radius : float
        Radius of spheres to be packed in container.
    center : Iterable of float
        Cartesian coordinates of the center of the container. Default is
        (0., 0., 0.)

    Attributes
    ----------
    width : float
        Prism length along the x-axis
    depth : float
        Prism length along the y-axis
    height : float
        Prism length along the z-axis
    sphere_radius : float
        Radius of spheres to be packed in container.
    center : list of float
        Cartesian coordinates of the center of the container. Default is
        (0., 0., 0.)
    cell_length : list of float
        Length in x-, y-, and z- directions of each cell in mesh overlaid on
        domain.
    limits : list of float
        Minimum and maximum distance in x-, y-, and z-direction where sphere
        center can be placed.
    volume : float
        Volume of the container.

    """

    def __init__(self, width, depth, height, sphere_radius, center=(0., 0., 0.)):
        super().__init__(sphere_radius, center)
        self.width = width
        self.depth = depth
        self.height = height

    @property
    def width(self):
        return self._width

    @property
    def depth(self):
        return self._depth

    @property
    def height(self):
        return self._height

    @property
    def limits(self):
        if self._limits is None:
            c = self.center
            r = self.sphere_radius
            x, y, z = self.width/2, self.depth/2, self.height/2
            self._limits = [[c[0] - x + r, c[1] - y + r, c[2] - z + r],
                            [c[0] + x - r, c[1] + y - r, c[2] + z - r]]
        return self._limits

    @property
    def cell_length(self):
        if self._cell_length is None:
            mesh_length = [self.width, self.depth, self.height]
            self._cell_length = [x/int(x/(4*self.sphere_radius))
                                 for x in mesh_length]
        return self._cell_length

    @property
    def volume(self):
        return self.width*self.depth*self.height

    @width.setter
    def width(self, width):
        self._width = float(width)
        self._limits = None
        self._cell_length = None

    @depth.setter
    def depth(self, depth):
        self._depth = float(depth)
        self._limits = None
        self._cell_length = None

    @height.setter
    def height(self, height):
        self._height = float(height)
        self._limits = None
        self._cell_length = None

    @limits.setter
    def limits(self, limits):
        self._limits = limits

    @classmethod
    def from_region(self, region, sphere_radius):
        check_type('region', region, openmc.Region)

        # Assume the simplest case where the prism volume is the intersection
        # of the half-spaces of six planes
        if not isinstance(region, openmc.Intersection):
            raise ValueError

        if any(not isinstance(node, openmc.Halfspace) for node in region):
            raise ValueError

        if len(region) != 6:
            raise ValueError

        # Sort half-spaces by surface type
        px1, px2, py1, py2, pz1, pz2 = sorted(region, key=lambda x: x.surface.type)

        # Make sure the region consists of the correct surfaces
        if (not isinstance(px1.surface, openmc.XPlane) or
            not isinstance(px2.surface, openmc.XPlane) or
            not isinstance(py1.surface, openmc.YPlane) or
            not isinstance(py2.surface, openmc.YPlane) or
            not isinstance(pz1.surface, openmc.ZPlane) or
            not isinstance(pz2.surface, openmc.ZPlane)):
            raise ValueError

        # Make sure the half-spaces are on the correct side of the surfaces
        ll, ur = region.bounding_box
        if any(x in ll or x in ur for x in (-np.inf, np.inf)):
            raise ValueError

        # Calculate the parameters for the container
        width, depth, height = ur - ll
        center = ll + [width/2, depth/2, height/2]

        # The region is the volume of a rectangular prism, so create container
        return _RectangularPrism(width, depth, height, sphere_radius, center)

    def random_point(self):
        ll, ul = self.limits
        return [uniform(ll[0], ul[0]),
                uniform(ll[1], ul[1]),
                uniform(ll[2], ul[2])]

    def repel_spheres(self, p, q, d, d_new):
        # Moving each sphere distance 's' away from the other along the line
        # joining the sphere centers will ensure their final distance is
        # equal to the outer diameter
        s = (d_new - d)/2

        v = (p - q)/d
        p += s*v
        q -= s*v

        # Enforce the rigid boundary by moving each sphere back along the
        # surface normal until it is completely within the container if it
        # overlaps the surface
        p[:] = np.clip(p, self.limits[0], self.limits[1])
        q[:] = np.clip(q, self.limits[0], self.limits[1])


class _Cylinder(_Container):
    """Cylindrical container in which to pack spheres.

    Parameters
    ----------
    length : float
        Length of the cylindrical container.
    radius : float
        Radius of the cylindrical container.
    axis : string
        Axis along which the length of the cylinder is aligned.
    sphere_radius : float
        Radius of spheres to be packed in container.
    center : Iterable of float
        Cartesian coordinates of the center of the container. Default is
        (0., 0., 0.)

    Attributes
    ----------
    length : float
        Length of the cylindrical container.
    radius : float
        Radius of the cylindrical container.
    axis : string
        Axis along which the length of the cylinder is aligned.
    sphere_radius : float
        Radius of spheres to be packed in container.
    center : list of float
        Cartesian coordinates of the center of the container. Default is
        (0., 0., 0.)
    shift : list of int
        Rolled indices of the x-, y-, and z- coordinates of a sphere so the
        configuration is aligned with the correct axis. No shift corresponds to
        a cylinder along the z-axis.
    cell_length : list of float
        Length in x-, y-, and z- directions of each cell in mesh overlaid on
        domain.
    limits : list of float
        Maximum radial distance and minimum and maximum distance in the
        direction parallel to the axis where sphere center can be placed.
    volume : float
        Volume of the container.

    """

    def __init__(self, length, radius, axis, sphere_radius, center=(0., 0., 0.)):
        super().__init__(sphere_radius, center)
        self._shift = None
        self.length = length
        self.radius = radius
        self.axis = axis

    @property
    def length(self):
        return self._length

    @property
    def radius(self):
        return self._radius

    @property
    def axis(self):
        return self._axis

    @property
    def shift(self):
        if self._shift is None:
            if self.axis == 'x':
                self._shift = [1, 2, 0]
            elif self.axis == 'y':
                self._shift = [2, 0, 1]
            else:
                self._shift = [0, 1, 2]
        return self._shift

    @property
    def limits(self):
        if self._limits is None:
            z0 = self.center[self.shift[2]]
            z = self.length/2
            r = self.sphere_radius
            self._limits = [[z0 - z + r], [z0 + z - r, self.radius - r]]
        return self._limits

    @property
    def cell_length(self):
        if self._cell_length is None:
            h = 4*self.sphere_radius
            i, j, k = self.shift
            self._cell_length = [None]*3
            self._cell_length[i] = 2*self.radius/int(2*self.radius/h)
            self._cell_length[j] = 2*self.radius/int(2*self.radius/h)
            self._cell_length[k] = self.length/int(self.length/h)
        return self._cell_length

    @property
    def volume(self):
        return self.length*pi*self.radius**2

    @length.setter
    def length(self, length):
        self._length = float(length)
        self._limits = None
        self._cell_length = None

    @radius.setter
    def radius(self, radius):
        self._radius = float(radius)
        self._limits = None
        self._cell_length = None

    @axis.setter
    def axis(self, axis):
        self._axis = axis
        self._shift = None

    @limits.setter
    def limits(self, limits):
        self._limits = limits

    @classmethod
    def from_region(self, region, sphere_radius):
        check_type('region', region, openmc.Region)

        # Assume the simplest case where the cylinder volume is the
        # intersection of the half-spaces of a cylinder and two planes
        if not isinstance(region, openmc.Intersection):
            raise ValueError

        if any(not isinstance(node, openmc.Halfspace) for node in region):
            raise ValueError

        if len(region) != 3:
            raise ValueError

        # Identify the axis that the cylinder lies along
        axis = region[0].surface.type[0]

        # Make sure the region is composed of a cylinder and two planes on the
        # same axis
        count = Counter(node.surface.type for node in region)
        if count[axis + '-cylinder'] != 1 or count[axis + '-plane'] != 2:
            raise ValueError

        # Sort the half-spaces by surface type
        cyl, p1, p2 = sorted(region, key=lambda x: x.surface.type)

        # Calculate the parameters for a cylinder along the x-axis
        if axis == 'x':
            if p1.surface.x0 > p2.surface.x0:
                p1, p2 = p2, p1
            length = p2.surface.x0 - p1.surface.x0
            center = (p1.surface.x0 + length/2, cyl.surface.y0, cyl.surface.z0)

        # Calculate the parameters for a cylinder along the y-axis
        elif axis == 'y':
            if p1.surface.y0 > p2.surface.y0:
                p1, p2 = p2, p1
            length = p2.surface.y0 - p1.surface.y0
            center = (cyl.surface.x0, p1.surface.y0 + length/2, cyl.surface.z0)

        # Calculate the parameters for a cylinder along the z-axis
        else:
            if p1.surface.z0 > p2.surface.z0:
                p1, p2 = p2, p1
            length = p2.surface.z0 - p1.surface.z0
            center = (cyl.surface.x0, cyl.surface.y0, p1.surface.z0 + length/2)

        # Make sure the half-spaces are on the correct side of the surfaces
        if cyl.side != '-' or p1.side != '+' or p2.side != '-':
            raise ValueError

        radius = cyl.surface.r

        # The region is the volume of a cylinder, so create container
        return _Cylinder(length, radius, axis, sphere_radius, center)

    def random_point(self):
        ll, ul = self.limits
        r = sqrt(uniform(0, ul[1]**2))
        t = uniform(0, 2*pi)
        i, j, k = self.shift
        p = [None]*3
        p[i] = r*cos(t) + self.center[i]
        p[j] = r*sin(t) + self.center[j]
        p[k] = uniform(ll[0], ul[0])
        return p

    def repel_spheres(self, p, q, d, d_new):
        # Moving each sphere distance 's' away from the other along the line
        # joining the sphere centers will ensure their final distance is
        # equal to the outer diameter
        s = (d_new - d)/2

        v = (p - q)/d
        p += s*v
        q -= s*v

        # Enforce the rigid boundary by moving each sphere back along the
        # surface normal until it is completely within the container if it
        # overlaps the surface
        ll, ul = self.limits
        c = self.center
        i, j, k = self.shift

        r = sqrt((p[i] - c[i])**2 + (p[j] - c[j])**2)
        if r > ul[1]:
            p[i] = (p[i] - c[i])*ul[1]/r + c[i]
            p[j] = (p[j] - c[j])*ul[1]/r + c[j]
        p[k] = np.clip(p[k], ll[0], ul[0])

        r = sqrt((q[i] - c[i])**2 + (q[j] - c[j])**2)
        if r > ul[1]:
            q[i] = (q[i] - c[i])*ul[1]/r + c[i]
            q[j] = (q[j] - c[j])*ul[1]/r + c[j]
        q[k] = np.clip(q[k], ll[0], ul[0])


class _SphericalShell(_Container):
    """Spherical shell container in which to pack spheres.

    Parameters
    ----------
    radius : float
        Outer radius of the spherical shell container.
    inner_radius : float
        Inner radius of the spherical shell container.
    center : Iterable of float
        Cartesian coordinates of the center of the container. Default is
        (0., 0., 0.)

    Attributes
    ----------
    radius : float
        Outer radius of the spherical shell container.
    inner_radius : float
        Inner radius of the spherical shell container.
    sphere_radius : float
        Radius of spheres to be packed in container.
    center : list of float
        Cartesian coordinates of the center of the container. Default is
        (0., 0., 0.)
    cell_length : list of float
        Length in x-, y-, and z- directions of each cell in mesh overlaid on
        domain.
    limits : list of float
        Minimum and maximum radial distance where sphere center can be placed.
    volume : float
        Volume of the container.

    """

    def __init__(self, radius, inner_radius, sphere_radius,
                 center=(0., 0., 0.)):
        super().__init__(sphere_radius, center)
        self.radius = radius
        self.inner_radius = inner_radius

    @property
    def radius(self):
        return self._radius

    @property
    def inner_radius(self):
        return self._inner_radius

    @property
    def limits(self):
        if self._limits is None:
            r_max = self.radius - self.sphere_radius
            if self.inner_radius == 0:
                r_min = 0
            else:
                r_min = self.inner_radius + self.sphere_radius
            self._limits = [[r_min], [r_max]]
        return self._limits

    @property
    def cell_length(self):
        if self._cell_length is None:
            mesh_length = 3*[2*self.radius]
            self._cell_length = [x/int(x/(4*self.sphere_radius))
                                 for x in mesh_length]
        return self._cell_length

    @property
    def volume(self):
        return 4/3*pi*(self.radius**3 - self.inner_radius**3)

    @radius.setter
    def radius(self, radius):
        self._radius = float(radius)
        self._limits = None
        self._cell_length = None

    @inner_radius.setter
    def inner_radius(self, inner_radius):
        self._inner_radius = float(inner_radius)
        self._limits = None

    @limits.setter
    def limits(self, limits):
        self._limits = limits

    @classmethod
    def from_region(self, region, sphere_radius):
        check_type('region', region, openmc.Region)

        # First check if the region is the volume inside a sphere. Assume the
        # simplest case where the sphere volume is the negative half-space of a
        # sphere.
        if (isinstance(region, openmc.Halfspace)
            and isinstance(region.surface, openmc.Sphere)
            and region.side == '-'):

            # The region is the volume of a sphere, so create container
            radius = region.surface.r
            center = (region.surface.x0, region.surface.y0, region.surface.z0)

            return _SphericalShell(radius, 0., sphere_radius, center)

        # Next check for a spherical shell volume. Assume the simplest case
        # where the spherical shell volume is the intersection of the
        # half-spaces of two spheres.
        if not isinstance(region, openmc.Intersection):
            raise ValueError

        if any(not isinstance(node, openmc.Halfspace) for node in region):
            raise ValueError

        if len(region) != 2:
            raise ValueError

        if any(not isinstance(node.surface, openmc.Sphere) for node in region):
            raise ValueError

        s1, s2 = sorted(region, key=lambda x: x.surface.r)
        radius = s2.surface.r
        inner_radius = s1.surface.r
        center = (s1.surface.x0, s1.surface.y0, s1.surface.z0)

        if center != (s2.surface.x0, s2.surface.y0, s2.surface.z0):
            raise ValueError

        if s1.side != '+' or s2.side != '-':
            raise ValueError

        # The region is the volume of a spherical shell, so create container
        return _SphericalShell(radius, inner_radius, sphere_radius, center)

    def random_point(self):
        c = self.center
        ll, ul = self.limits
        x, y, z = (gauss(0, 1), gauss(0, 1), gauss(0, 1))
        r = (uniform(ll[0]**3, ul[0]**3)**(1/3)/sqrt(x**2 + y**2 + z**2))
        return [r*x + c[0],  r*y + c[1], r*z + c[2]]

    def repel_spheres(self, p, q, d, d_new):
        # Moving each sphere distance 's' away from the other along the line
        # joining the sphere centers will ensure their final distance is
        # equal to the outer diameter
        s = (d_new - d)/2

        v = (p - q)/d
        p += s*v
        q -= s*v

        # Enforce the rigid boundary by moving each sphere back along the
        # surface normal until it is completely within the container if it
        # overlaps the surface
        c = self.center
        ll, ul = self.limits

        r = sqrt((p[0] - c[0])**2 + (p[1] - c[1])**2 + (p[2] - c[2])**2)
        if r > ul[0]:
            p[:] = (p - c)*ul[0]/r + c
        elif r < ll[0]:
            p[:] = (p - c)*ll[0]/r + c

        r = sqrt((q[0] - c[0])**2 + (q[1] - c[1])**2 + (q[2] - c[2])**2)
        if r > ul[0]:
            q[:] = (q - c)*ul[0]/r + c
        elif r < ll[0]:
            q[:] = (q - c)*ll[0]/r + c


def create_triso_lattice(trisos, lower_left, pitch, shape, background):
    """Create a lattice containing TRISO particles for optimized tracking.

    Parameters
    ----------
    trisos : list of openmc.model.TRISO
        List of TRISO particles to put in lattice
    lower_left : Iterable of float
        Lower-left Cartesian coordinates of the lattice
    pitch : Iterable of float
        Pitch of the lattice elements in the x-, y-, and z-directions
    shape : Iterable of float
        Number of lattice elements in the x-, y-, and z-directions
    background : openmc.Material
        A background material that is used anywhere within the lattice but
        outside a TRISO particle

    Returns
    -------
    lattice : openmc.RectLattice
        A lattice containing the TRISO particles

    """

    lattice = openmc.RectLattice()
    lattice.lower_left = lower_left
    lattice.pitch = pitch

    indices = list(np.broadcast(*np.ogrid[:shape[2], :shape[1], :shape[0]]))
    triso_locations = {idx: [] for idx in indices}
    for t in trisos:
        for idx in t.classify(lattice):
            if idx in triso_locations:
                # Create copy of TRISO particle with materials preserved and
                # different cell/surface IDs
                t_copy = copy.copy(t)
                t_copy.id = None
                t_copy.fill = t.fill
                t_copy._surface = openmc.Sphere(r=t._surface.r,
                                                x0=t._surface.x0,
                                                y0=t._surface.y0,
                                                z0=t._surface.z0)
                t_copy.region = -t_copy._surface
                triso_locations[idx].append(t_copy)
            else:
                warnings.warn('TRISO particle is partially or completely '
                              'outside of the lattice.')

    # Create universes
    universes = np.empty(shape[::-1], dtype=openmc.Universe)
    for idx, triso_list in sorted(triso_locations.items()):
        if len(triso_list) > 0:
            outside_trisos = openmc.Intersection(~t.region for t in triso_list)
            background_cell = openmc.Cell(fill=background, region=outside_trisos)
        else:
            background_cell = openmc.Cell(fill=background)

        u = openmc.Universe()
        u.add_cell(background_cell)
        for t in triso_list:
            u.add_cell(t)
            iz, iy, ix = idx
            t.center = lattice.get_local_coordinates(t.center, (ix, iy, iz))

        if len(shape) == 2:
            universes[-1 - idx[0], idx[1]] = u
        else:
            universes[idx[0], -1 - idx[1], idx[2]] = u
    lattice.universes = universes

    # Set outer universe
    background_cell = openmc.Cell(fill=background)
    lattice.outer = openmc.Universe(cells=[background_cell])

    return lattice


def _random_sequential_pack(domain, num_spheres):
    """Random sequential packing of spheres within a container.

    Parameters
    ----------
    domain : openmc.model._Container
        Container in which to pack spheres.
    num_spheres : int
        Number of spheres to pack.

    Returns
    ------
    numpy.ndarray
        Cartesian coordinates of centers of spheres.

    """

    sqd = (2*domain.sphere_radius)**2
    spheres = []
    mesh = defaultdict(list)

    for i in range(num_spheres):
        # Randomly sample new center coordinates while there are any overlaps
        while True:
            p = domain.random_point()
            idx = domain.mesh_cell(p)
            if any((p[0]-q[0])**2 + (p[1]-q[1])**2 + (p[2]-q[2])**2 < sqd
                   for q in mesh[idx]):
                continue
            else:
                break
        spheres.append(p)

        for idx in domain.nearby_mesh_cells(p):
            mesh[idx].append(p)

    return np.array(spheres)


def _close_random_pack(domain, spheres, contraction_rate):
    """Close random packing of spheres using the Jodrey-Tory algorithm.

    Parameters
    ----------
    domain : openmc.model._Container
        Container in which to pack spheres.
    spheres : numpy.ndarray
        Initial Cartesian coordinates of centers of spheres.
    contraction_rate : float
        Contraction rate of outer diameter.

    """

    def add_rod(d, i, j):
        """Add a new rod to the priority queue.

        Parameters
        ----------
        d : float
            distance between centers of spheres i and j.
        i, j : int
            Index of spheres in spheres array.

        """

        rod = [d, i, j]
        rods_map[i] = (j, rod)
        rods_map[j] = (i, rod)
        heappush(rods, rod)

    def remove_rod(i):
        """Mark the rod containing sphere i as removed.

        Parameters
        ----------
        i : int
            Index of sphere in spheres array.

        """

        if i in rods_map:
            j, rod = rods_map.pop(i)
            del rods_map[j]
            rod[1] = removed
            rod[2] = removed

    def pop_rod():
        """Remove and return the shortest rod.

        Returns
        -------
        d : float
            distance between centers of spheres i and j.
        i, j : int
            Index of spheres in spheres array.

        """

        while rods:
            d, i, j = heappop(rods)
            if i != removed and j != removed:
                del rods_map[i]
                del rods_map[j]
                return d, i, j
        return None, None, None

    def create_rod_list():
        """Generate sorted list of rods (distances between sphere centers).

        Rods are arranged in a heap where each element contains the rod length
        and the sphere indices. A rod between spheres p and q is only
        included if the distance between p and q could not be changed by the
        elimination of a greater overlap, i.e. q has no nearer neighbors than p.

        A mapping of sphere ids to rods is maintained in 'rods_map'. Each key
        in the dict is the id of a sphere that is in the rod list, and the
        value is the id of its nearest neighbor and the rod that contains them.
        The dict is used to find rods in the priority queue and to mark removed
        rods so rods can be "removed" without breaking the heap structure
        invariant.

        """

        # Create KD tree for quick nearest neighbor search
        tree = scipy.spatial.cKDTree(spheres)

        # Find distance to nearest neighbor and index of nearest neighbor for
        # all spheres
        d, n = tree.query(spheres, k=2)
        d = d[:, 1]
        n = n[:, 1]

        # Array of sphere indices, indices of nearest neighbors, and
        # distances to nearest neighbors
        a = np.vstack((list(range(n.size)), n, d)).T

        # Sort along second column and swap first and second columns to create
        # array of nearest neighbor indices, indices of spheres they are
        # nearest neighbors of, and distances between them
        b = a[a[:, 1].argsort()]
        b[:, [0, 1]] = b[:, [1, 0]]

        # Find the intersection between 'a' and 'b': a list of spheres who
        # are each other's nearest neighbors and the distance between them
        r = list({tuple(x) for x in a} & {tuple(x) for x in b})

        # Remove duplicate rods and sort by distance
        r = map(list, set([(x[2], int(min(x[0:2])), int(max(x[0:2])))
                            for x in r]))

        # Clear priority queue and add rods
        del rods[:]
        rods_map.clear()
        for d, i, j in r:
            if d < outer_diameter and not np.isclose(d, outer_diameter, atol=1.0e-14):
                add_rod(d, i, j)

    def update_mesh(i):
        """Update which mesh cells the sphere is in based on new sphere
        center coordinates.

        'mesh'/'mesh_map' is a two way dictionary used to look up which
        spheres are located within one diameter of a given mesh cell and
        which mesh cells a given sphere center is within one diameter of.
        This is used to speed up the nearest neighbor search.

        Parameters
        ----------
        i : int
            Index of sphere in spheres array.

        """

        # Determine which mesh cells the sphere is in and remove the
        # sphere id from those cells
        for idx in mesh_map[i]:
            mesh[idx].remove(i)
        del mesh_map[i]

        # Determine which mesh cells are within one diameter of sphere's
        # center and add this sphere to the list of spheres in those cells
        for idx in domain.nearby_mesh_cells(spheres[i]):
            mesh[idx].add(i)
            mesh_map[i].add(idx)

    def reduce_outer_diameter():
        """Reduce the outer diameter so that at the (i+1)-st iteration it is:

            d_out^(i+1) = d_out^(i) - (1/2)^(j) * d_out0 * k / n,

        where k is the contraction rate, n is the number of spheres, and

            j = floor(-log10(pf_out - pf_in)).

        Returns
        -------
        float
            New outer diameter

        """

        inner_pf = 4/3*pi*(inner_diameter/2)**3*num_spheres/domain.volume
        outer_pf = 4/3*pi*(outer_diameter/2)**3*num_spheres/domain.volume

        j = floor(-log10(outer_pf - inner_pf))
        return (outer_diameter - 0.5**j * contraction_rate *
                initial_outer_diameter / num_spheres)

    def nearest(i):
        """Find index of nearest neighbor of sphere i.

        Parameters
        ----------
        i : int
            Index in spheres array of sphere for which to find nearest
            neighbor.

        Returns
        -------
        int
            Index in spheres array of nearest neighbor of i
        float
            distance between i and nearest neighbor.

        """

        # Need the second nearest neighbor of i since the nearest neighbor
        # will be itself. Using argpartition, the k-th nearest neighbor is
        # placed at index k.
        idx = list(mesh[domain.mesh_cell(spheres[i])])
        dists = scipy.spatial.distance.cdist([spheres[i]], spheres[idx])[0]
        if dists.size > 1:
            j = dists.argpartition(1)[1]
            return idx[j], dists[j]
        else:
            return None, None

    def update_rod_list(i):
        """Update the rod list with the new nearest neighbors of sphere since
        its overlap was eliminated.

        Parameters
        ----------
        i : int
            Index of sphere in spheres array.

        """

        # If the nearest neighbor k of sphere i has no nearer neighbors,
        # remove the rod currently containing k from the rod list and add rod
        # k-i, keeping the rod list sorted
        k, d_ik = nearest(i)
        if (k and nearest(k)[0] == i and d_ik < outer_diameter
            and not np.isclose(d, outer_diameter, atol=1.0e-14)):
            remove_rod(k)
            add_rod(d_ik, i, k)

    num_spheres = len(spheres)
    diameter = 2*domain.sphere_radius

    # Flag for marking rods that have been removed from priority queue
    removed = -1

    # Outer diameter initially set to arbitrary value that yields pf of 1
    initial_outer_diameter = 2*(domain.volume/(num_spheres*4/3*pi))**(1/3)

    # Inner and outer diameter of spheres will change during packing
    outer_diameter = initial_outer_diameter
    inner_diameter = 0.

    # List of rods arranged in a heap and mapping of sphere ids to rods
    rods = []
    rods_map = {}

    # Initialize two-way dictionary that identifies which spheres are near a
    # given mesh cell and which mesh cells a sphere is near
    mesh = defaultdict(set)
    mesh_map = defaultdict(set)
    for i in range(num_spheres):
        for idx in domain.nearby_mesh_cells(spheres[i]):
            mesh[idx].add(i)
            mesh_map[i].add(idx)

    while True:
        # Rebuild the sorted list of rods according to the current sphere
        # configuration
        create_rod_list()

        # Set the inner diameter to the shortest center-to-center distance
        # between any two spheres
        if rods:
            inner_diameter = rods[0][0]

        # Reached the desired sphere radius
        if inner_diameter >= diameter:
            break

        # The algorithm converged before reaching the desired sphere radius.
        # This can happen when the desired packing fraction is close to the
        # packing fraction limit. The packing fraction is a random variable
        # that is determined by the sphere locations and the contraction
        # rate. A higher packing fraction can be achieved with a smaller
        # contraction rate, though at the cost of a longer simulation time --
        # the number of iterations needed to remove all overlaps is inversely
        # proportional to the contraction rate.
        if inner_diameter >= outer_diameter or not rods:
            warnings.warn('Close random pack converged before reaching true '
                          'sphere radius; some spheres may overlap. Try '
                          'reducing contraction rate or packing fraction.')
            break

        while True:
            d, i, j = pop_rod()
            if not d:
                break
            outer_diameter = reduce_outer_diameter()
            domain.repel_spheres(spheres[i], spheres[j], d, outer_diameter)
            update_mesh(i)
            update_mesh(j)
            update_rod_list(i)
            update_rod_list(j)
            if not rods:
                break
            inner_diameter = rods[0][0]
            if inner_diameter >= diameter or inner_diameter >= outer_diameter:
                break


def pack_spheres(radius, region, pf=None, num_spheres=None, initial_pf=0.3,
                 contraction_rate=1.e-3, seed=1):
    """Generate a random, non-overlapping configuration of spheres within a
    container.

    Parameters
    ----------
    radius : float
        Outer radius of spheres.
    region : openmc.Region
        Container in which the spheres are packed. Supported shapes are
        rectangular prism, cylinder, sphere, and spherical shell.
    pf : float
        Packing fraction of the spheres. One of 'pf' and 'num_spheres' must
        be specified; the other will be calculated. If both are specified, 'pf'
        takes precedence over 'num_spheres'.
    num_spheres : int
        Number of spheres to pack in the domain. One of 'num_spheres' and 'pf'
        must be specified; the other will be calculated.
    initial_pf : float, optional
        Packing fraction used to initialize the configuration of spheres in
        the domain. Default value is 0.3. It is not recommended to set the
        initial packing fraction much higher than 0.3 as the random sequential
        packing algorithm becomes prohibitively slow as it approaches its limit
        (~0.38).
    contraction_rate : float, optional
        Contraction rate of the outer diameter. Higher packing fractions can be
        reached using a smaller contraction rate, but the algorithm will take
        longer to converge.
    seed : int, optional
        RNG seed.

    Returns
    ------
    numpy.ndarray
        Cartesian coordinates of sphere centers.

    Notes
    -----
    The sphere configuration is generated using a combination of random
    sequential packing (RSP) and close random packing (CRP). RSP performs
    better than CRP for lower packing fractions (pf), but it becomes
    prohibitively slow as it approaches its packing limit (~0.38). CRP can
    achieve higher pf of up to ~0.64 and scales better with increasing pf.

    If the desired pf is below some threshold for which RSP will be faster than
    CRP ('initial_packing_fraction'), only RSP is used. If a higher pf is
    required, spheres with a radius smaller than the desired final radius
    (and therefore with a smaller pf) are initialized within the domain using
    RSP. This initial configuration of spheres is then used as a starting
    point for CRP using Jodrey and Tory's algorithm [1]_.

    In RSP, sphere centers are placed one by one at random, and placement
    attempts for a sphere are made until the sphere is not overlapping any
    others. This implementation of the algorithm uses a mesh over the domain
    to speed up the nearest neighbor search by only searching for a sphere's
    neighbors within that mesh cell.

    In CRP, each sphere is assigned two diameters, an inner and an outer,
    which approach each other during the simulation. The inner diameter,
    defined as the minimum center-to-center distance, is the true diameter of
    the spheres and defines the pf. At each iteration the worst overlap
    between spheres based on outer diameter is eliminated by moving the
    spheres apart along the line joining their centers. Iterations continue
    until the two diameters converge or until the desired pf is reached.

    References
    ----------
    .. [1] W. S. Jodrey and E. M. Tory, "Computer simulation of close random
       packing of equal spheres", Phys. Rev. A 32 (1985) 2347-2351.

    """
    # Seed RNG
    random.seed(seed)

    # Create container with the correct shape based on the supplied region
    domain = None
    for cls in _Container.__subclasses__():
        try:
            domain = cls.from_region(region, radius)
        except ValueError:
            pass

    if not domain:
        raise ValueError('Could not map region {} to a container: supported '
                         'container shapes are rectangular prism, cylinder, '
                         'sphere, and spherical shell.'.format(region))

    # Determine the packing fraction/number of spheres
    volume = 4/3*pi*radius**3
    if pf is None and num_spheres is None:
        raise ValueError('`pf` or `num_spheres` must be specified.')
    elif pf is None:
        num_spheres = int(num_spheres)
        pf = volume*num_spheres/domain.volume
    else:
        pf = float(pf)
        num_spheres = int(pf*domain.volume//volume)

    # Make sure initial packing fraction is less than packing fraction
    if initial_pf > pf:
        initial_pf = pf

    # Check packing fraction for close random packing
    if pf > MAX_PF_CRP:
        raise ValueError('Packing fraction {0} is greater than the limit for '
                         'close random packing, {1}'.format(pf, MAX_PF_CRP))

    # Check packing fraction for random sequential packing
    if initial_pf > MAX_PF_RSP:
        raise ValueError('Initial packing fraction {0} is greater than the '
                         'limit for random sequential packing, '
                         '{1}'.format(initial_pf, MAX_PF_RSP))

    # Calculate the sphere radius used in the initial random sequential
    # packing from the initial packing fraction
    initial_radius = (3/4*initial_pf*domain.volume/(pi*num_spheres))**(1/3)
    domain.sphere_radius = initial_radius

    # Recalculate the limits for the initial random sequential packing using
    # the desired final sphere radius to ensure spheres are fully contained
    # within the domain during the close random pack
    domain.limits = [[x - initial_radius + radius for x in domain.limits[0]],
                     [x + initial_radius - radius for x in domain.limits[1]]]

    # Generate non-overlapping spheres for an initial inner radius using
    # random sequential packing algorithm
    spheres = _random_sequential_pack(domain, num_spheres)

    # Use the sphere configuration produced in random sequential packing as a
    # starting point for close random pack with the desired final sphere
    # radius
    if initial_pf != pf:
        domain.sphere_radius = radius
        _close_random_pack(domain, spheres, contraction_rate)

    return spheres
