import copy
import warnings
import itertools
import random
from collections import defaultdict
from collections.abc import Iterable
from numbers import Real
from random import uniform, gauss
from heapq import heappush, heappop
from math import pi, sin, cos, floor, log10, sqrt
from abc import ABCMeta, abstractproperty, abstractmethod

import numpy as np
import scipy.spatial

import openmc
import openmc.checkvalue as cv


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
        self._surface = openmc.Sphere(R=outer_radius)
        super().__init__(fill=fill, region=-self._surface)
        self.center = np.asarray(center)

    @property
    def center(self):
        return self._center

    @center.setter
    def center(self, center):
        cv.check_type('TRISO center', center, Iterable, Real)
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


class _Domain(metaclass=ABCMeta):
    """Container in which to pack particles.

    Parameters
    ----------
    particle_radius : float
        Radius of particles to be packed in container.
    center : Iterable of float
        Cartesian coordinates of the center of the container. Default is
        [0., 0., 0.]

    Attributes
    ----------
    particle_radius : float
        Radius of particles to be packed in container.
    center : list of float
        Cartesian coordinates of the center of the container. Default is
        [0., 0., 0.]
    cell_length : list of float
        Length in x-, y-, and z- directions of each cell in mesh overlaid on
        domain.
    limits : list of float
        Minimum and maximum position in x-, y-, and z-directions where particle
        center can be placed.
    volume : float
        Volume of the container.

    """
    def __init__(self, particle_radius, center=[0., 0., 0.]):
        self._cell_length = None
        self._limits = None

        self.particle_radius = particle_radius
        self.center = center

    @property
    def particle_radius(self):
        return self._particle_radius

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

    @particle_radius.setter
    def particle_radius(self, particle_radius):
        self._particle_radius = float(particle_radius)
        self._limits = None
        self._cell_length = None

    @center.setter
    def center(self, center):
        if np.asarray(center).size != 3:
            raise ValueError('Unable to set domain center to {} since it must '
                             'be of length 3'.format(center))
        self._center = [float(x) for x in center]
        self._limits = None
        self._cell_length = None

    def mesh_cell(self, p):
        """Calculate the index of the cell in a mesh overlaid on the domain in
        which the given particle center falls.

        Parameters
        ----------
        p : Iterable of float
            Cartesian coordinates of particle center.

        Returns
        -------
        tuple of int
            Indices of mesh cell.

        """
        return tuple(int(p[i]/self.cell_length[i]) for i in range(3))

    def nearby_mesh_cells(self, p):
        """Calculates the indices of all cells in a mesh overlaid on the domain
        within one diameter of the given particle.

        Parameters
        ----------
        p : Iterable of float
            Cartesian coordinates of particle center.

        Returns
        -------
        list of tuple of int
            Indices of mesh cells.

        """
        d = 2*self.particle_radius
        r = [[a/self.cell_length[i] for a in [p[i]-d, p[i], p[i]+d]]
             for i in range(3)]
        return list(itertools.product(*({int(x) for x in y} for y in r)))

    @abstractmethod
    def random_point(self):
        """Generate Cartesian coordinates of center of a particle that is
        contained entirely within the domain with uniform probability.

        Returns
        -------
        list of float
            Cartesian coordinates of particle center.

        """
        pass


class _CubicDomain(_Domain):
    """Cubic container in which to pack particles.

    Parameters
    ----------
    length : float
        Length of each side of the cubic container.
    particle_radius : float
        Radius of particles to be packed in container.
    center : Iterable of float
        Cartesian coordinates of the center of the container. Default is
        [0., 0., 0.]

    Attributes
    ----------
    length : float
        Length of each side of the cubic container.
    particle_radius : float
        Radius of particles to be packed in container.
    center : list of float
        Cartesian coordinates of the center of the container. Default is
        [0., 0., 0.]
    cell_length : list of float
        Length in x-, y-, and z- directions of each cell in mesh overlaid on
        domain.
    limits : list of float
        Minimum and maximum position in x-, y-, and z-directions where particle
        center can be placed.
    volume : float
        Volume of the container.

    """

    def __init__(self, length, particle_radius, center=[0., 0., 0.]):
        super().__init__(particle_radius, center)
        self.length = length

    @property
    def length(self):
        return self._length

    @property
    def limits(self):
        if self._limits is None:
            xlim = self.length/2 - self.particle_radius
            self._limits = [[x - xlim for x in self.center],
                            [x + xlim for x in self.center]]
        return self._limits

    @property
    def cell_length(self):
        if self._cell_length is None:
            mesh_length = [self.length, self.length, self.length]
            self._cell_length = [x/int(x/(4*self.particle_radius))
                                 for x in mesh_length]
        return self._cell_length

    @property
    def volume(self):
        return self.length**3

    @length.setter
    def length(self, length):
        self._length = float(length)
        self._limits = None
        self._cell_length = None

    @limits.setter
    def limits(self, limits):
        self._limits = limits

    def random_point(self):
        return [uniform(self.limits[0][0], self.limits[1][0]),
                uniform(self.limits[0][1], self.limits[1][1]),
                uniform(self.limits[0][2], self.limits[1][2])]


class _CylindricalDomain(_Domain):
    """Cylindrical container in which to pack particles.

    Parameters
    ----------
    length : float
        Length along z-axis of the cylindrical container.
    radius : float
        Radius of the cylindrical container.
    center : Iterable of float
        Cartesian coordinates of the center of the container. Default is
        [0., 0., 0.]

    Attributes
    ----------
    length : float
        Length along z-axis of the cylindrical container.
    radius : float
        Radius of the cylindrical container.
    particle_radius : float
        Radius of particles to be packed in container.
    center : list of float
        Cartesian coordinates of the center of the container. Default is
        [0., 0., 0.]
    cell_length : list of float
        Length in x-, y-, and z- directions of each cell in mesh overlaid on
        domain.
    limits : list of float
        Minimum and maximum position in x-, y-, and z-directions where particle
        center can be placed.
    volume : float
        Volume of the container.

    """

    def __init__(self, length, radius, particle_radius, center=[0., 0., 0.]):
        super().__init__(particle_radius, center)
        self.length = length
        self.radius = radius

    @property
    def length(self):
        return self._length

    @property
    def radius(self):
        return self._radius

    @property
    def limits(self):
        if self._limits is None:
            xlim = self.length/2 - self.particle_radius
            rlim = self.radius - self.particle_radius
            self._limits = [[self.center[0] - rlim, self.center[1] - rlim,
                             self.center[2] - xlim],
                            [self.center[0] + rlim, self.center[1] + rlim,
                             self.center[2] + xlim]]
        return self._limits

    @property
    def cell_length(self):
        if self._cell_length is None:
            mesh_length = [2*self.radius, 2*self.radius, self.length]
            self._cell_length = [x/int(x/(4*self.particle_radius))
                                 for x in mesh_length]
        return self._cell_length

    @property
    def volume(self):
        return self.length * pi * self.radius**2

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

    @limits.setter
    def limits(self, limits):
        self._limits = limits

    def random_point(self):
        r = sqrt(uniform(0, (self.radius - self.particle_radius)**2))
        t = uniform(0, 2*pi)
        return [r*cos(t) + self.center[0], r*sin(t) + self.center[1],
                uniform(self.limits[0][2], self.limits[1][2])]


class _SphericalDomain(_Domain):
    """Spherical container in which to pack particles.

    Parameters
    ----------
    radius : float
        Radius of the spherical container.
    center : Iterable of float
        Cartesian coordinates of the center of the container. Default is
        [0., 0., 0.]

    Attributes
    ----------
    radius : float
        Radius of the spherical container.
    particle_radius : float
        Radius of particles to be packed in container.
    center : list of float
        Cartesian coordinates of the center of the container. Default is
        [0., 0., 0.]
    cell_length : list of float
        Length in x-, y-, and z- directions of each cell in mesh overlaid on
        domain.
    limits : list of float
        Minimum and maximum position in x-, y-, and z-directions where particle
        center can be placed.
    volume : float
        Volume of the container.

    """

    def __init__(self, radius, particle_radius, center=[0., 0., 0.]):
        super().__init__(particle_radius, center)
        self.radius = radius

    @property
    def radius(self):
        return self._radius

    @property
    def limits(self):
        if self._limits is None:
            rlim = self.radius - self.particle_radius
            self._limits = [[x - rlim for x in self.center],
                            [x + rlim for x in self.center]]
        return self._limits

    @property
    def cell_length(self):
        if self._cell_length is None:
            mesh_length = [2*self.radius, 2*self.radius, 2*self.radius]
            self._cell_length = [x/int(x/(4*self.particle_radius))
                                 for x in mesh_length]
        return self._cell_length

    @property
    def volume(self):
        return 4/3 * pi * self.radius**3

    @radius.setter
    def radius(self, radius):
        self._radius = float(radius)
        self._limits = None
        self._cell_length = None

    @limits.setter
    def limits(self, limits):
        self._limits = limits

    def random_point(self):
        x = (gauss(0, 1), gauss(0, 1), gauss(0, 1))
        r = (uniform(0, (self.radius - self.particle_radius)**3)**(1/3) /
             sqrt(x[0]**2 + x[1]**2 + x[2]**2))
        return [r*x[i] + self.center[i] for i in range(3)]


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
            if idx in sorted(triso_locations):
                # Create copy of TRISO particle with materials preserved and
                # different cell/surface IDs
                t_copy = copy.deepcopy(t)
                t_copy.id = None
                t_copy.fill = t.fill
                t_copy._surface.id = None
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


def _random_sequential_pack(domain, n_particles):
    """Random sequential packing of particles within a container.

    Parameters
    ----------
    domain : openmc.model._Domain
        Container in which to pack particles.
    n_particles : int
        Number of particles to pack.

    Returns
    ------
    numpy.ndarray
        Cartesian coordinates of centers of particles.

    """

    sqd = (2*domain.particle_radius)**2
    particles = []
    mesh = defaultdict(list)

    for i in range(n_particles):
        # Randomly sample new center coordinates while there are any overlaps
        while True:
            p = domain.random_point()
            idx = domain.mesh_cell(p)
            if any((p[0]-q[0])**2 + (p[1]-q[1])**2 + (p[2]-q[2])**2 < sqd
                   for q in mesh[idx]):
                continue
            else:
                break
        particles.append(p)

        for idx in domain.nearby_mesh_cells(p):
            mesh[idx].append(p)

    return np.array(particles)


def _close_random_pack(domain, particles, contraction_rate):
    """Close random packing of particles using the Jodrey-Tory algorithm.

    Parameters
    ----------
    domain : openmc.model._Domain
        Container in which to pack particles.
    particles : numpy.ndarray
        Initial Cartesian coordinates of centers of particles.
    contraction_rate : float
        Contraction rate of outer diameter.

    """

    def add_rod(d, i, j):
        """Add a new rod to the priority queue.

        Parameters
        ----------
        d : float
            distance between centers of particles i and j.
        i, j : int
            Index of particles in particles array.

        """

        rod = [d, i, j]
        rods_map[i] = (j, rod)
        rods_map[j] = (i, rod)
        heappush(rods, rod)

    def remove_rod(i):
        """Mark the rod containing particle i as removed.

        Parameters
        ----------
        i : int
            Index of particle in particles array.

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
            distance between centers of particles i and j.
        i, j : int
            Index of particles in particles array.

        """

        while rods:
            d, i, j = heappop(rods)
            if i != removed and j != removed:
                del rods_map[i]
                del rods_map[j]
                return d, i, j

    def create_rod_list():
        """Generate sorted list of rods (distances between particle centers).

        Rods are arranged in a heap where each element contains the rod length
        and the particle indices. A rod between particles p and q is only
        included if the distance between p and q could not be changed by the
        elimination of a greater overlap, i.e. q has no nearer neighbors than p.

        A mapping of particle ids to rods is maintained in 'rods_map'. Each key
        in the dict is the id of a particle that is in the rod list, and the
        value is the id of its nearest neighbor and the rod that contains them.
        The dict is used to find rods in the priority queue and to mark removed
        rods so rods can be "removed" without breaking the heap structure
        invariant.

        """

        # Create KD tree for quick nearest neighbor search
        tree = scipy.spatial.cKDTree(particles)

        # Find distance to nearest neighbor and index of nearest neighbor for
        # all particles
        d, n = tree.query(particles, k=2)
        d = d[:,1]
        n = n[:,1]

        # Array of particle indices, indices of nearest neighbors, and
        # distances to nearest neighbors
        a = np.vstack((list(range(n.size)), n, d)).T

        # Sort along second column and swap first and second columns to create
        # array of nearest neighbor indices, indices of particles they are
        # nearest neighbors of, and distances between them
        b = a[a[:,1].argsort()]
        b[:,[0, 1]] = b[:,[1, 0]]

        # Find the intersection between 'a' and 'b': a list of particles who
        # are each other's nearest neighbors and the distance between them
        r = list({tuple(x) for x in a} & {tuple(x) for x in b})

        # Remove duplicate rods and sort by distance
        r = map(list, set([(x[2], int(min(x[0:2])), int(max(x[0:2])))
                            for x in r]))

        # Clear priority queue and add rods
        del rods[:]
        rods_map.clear()
        for d, i, j in r:
            add_rod(d, i, j)

        # Inner diameter is set initially to the shortest center-to-center
        # distance between any two particles
        if rods:
            inner_diameter[0] = rods[0][0]

    def update_mesh(i):
        """Update which mesh cells the particle is in based on new particle
        center coordinates.

        'mesh'/'mesh_map' is a two way dictionary used to look up which
        particles are located within one diameter of a given mesh cell and
        which mesh cells a given particle center is within one diameter of.
        This is used to speed up the nearest neighbor search.

        Parameters
        ----------
        i : int
            Index of particle in particles array.

        """

        # Determine which mesh cells the particle is in and remove the
        # particle id from those cells
        for idx in mesh_map[i]:
            mesh[idx].remove(i)
        del mesh_map[i]

        # Determine which mesh cells are within one diameter of particle's
        # center and add this particle to the list of particles in those cells
        for idx in domain.nearby_mesh_cells(particles[i]):
            mesh[idx].add(i)
            mesh_map[i].add(idx)

    def reduce_outer_diameter():
        """Reduce the outer diameter so that at the (i+1)-st iteration it is:

            d_out^(i+1) = d_out^(i) - (1/2)^(j) * d_out0 * k / n,

        where k is the contraction rate, n is the number of particles, and

            j = floor(-log10(pf_out - pf_in)).

        """

        inner_pf = (4/3 * pi * (inner_diameter[0]/2)**3 * n_particles /
                    domain.volume)
        outer_pf = (4/3 * pi * (outer_diameter[0]/2)**3 * n_particles /
                    domain.volume)

        j = floor(-log10(outer_pf - inner_pf))
        outer_diameter[0] = (outer_diameter[0] - 0.5**j * contraction_rate *
                             initial_outer_diameter / n_particles)


    def repel_particles(i, j, d):
        """Move particles p and q apart according to the following
        transformation (accounting for reflective boundary conditions on
        domain):

            r_i^(n+1) = r_i^(n) + 1/2(d_out^(n+1) - d^(n))
            r_j^(n+1) = r_j^(n) - 1/2(d_out^(n+1) - d^(n))

        Parameters
        ----------
        i, j : int
            Index of particles in particles array.
        d : float
            distance between centers of particles i and j.

        """

        # Moving each particle distance 'r' away from the other along the line
        # joining the particle centers will ensure their final distance is equal
        # to the outer diameter
        r = (outer_diameter[0] - d)/2

        v = (particles[i] - particles[j])/d
        particles[i] += r*v
        particles[j] -= r*v

        # Apply reflective boundary conditions
        particles[i] = particles[i].clip(domain.limits[0], domain.limits[1])
        particles[j] = particles[j].clip(domain.limits[0], domain.limits[1])

        update_mesh(i)
        update_mesh(j)

    def nearest(i):
        """Find index of nearest neighbor of particle i.

        Parameters
        ----------
        i : int
            Index in particles array of particle for which to find nearest
            neighbor.

        Returns
        -------
        int
            Index in particles array of nearest neighbor of i
        float
            distance between i and nearest neighbor.

        """

        # Need the second nearest neighbor of i since the nearest neighbor
        # will be itself. Using argpartition, the k-th nearest neighbor is
        # placed at index k.
        idx = list(mesh[domain.mesh_cell(particles[i])])
        dists = scipy.spatial.distance.cdist([particles[i]], particles[idx])[0]
        if dists.size > 1:
            j = dists.argpartition(1)[1]
            return idx[j], dists[j]
        else:
            return None, None

    def update_rod_list(i, j):
        """Update the rod list with the new nearest neighbors of particles i
        and j since their overlap was eliminated.

        Parameters
        ----------
        i, j : int
            Index of particles in particles array.

        """

        # If the nearest neighbor k of particle i has no nearer neighbors,
        # remove the rod currently containing k from the rod list and add rod
        # k-i, keeping the rod list sorted
        k, d_ik = nearest(i)
        if k and nearest(k)[0] == i:
            remove_rod(k)
            add_rod(d_ik, i, k)
        l, d_jl = nearest(j)
        if l and nearest(l)[0] == j:
            remove_rod(l)
            add_rod(d_jl, j, l)

        # Set inner diameter to the shortest distance between two particle
        # centers
        if rods:
            inner_diameter[0] = rods[0][0]

    n_particles = len(particles)
    diameter = 2*domain.particle_radius

    # Flag for marking rods that have been removed from priority queue
    removed = -1

    # Outer diameter initially set to arbitrary value that yields pf of 1
    initial_outer_diameter = 2*(domain.volume/(n_particles*4/3*pi))**(1/3)

    # Inner and outer diameter of particles will change during packing
    outer_diameter = [initial_outer_diameter]
    inner_diameter = [0]

    rods = []
    rods_map = {}
    mesh = defaultdict(set)
    mesh_map = defaultdict(set)

    for i in range(n_particles):
        for idx in domain.nearby_mesh_cells(particles[i]):
            mesh[idx].add(i)
            mesh_map[i].add(idx)

    while True:
        create_rod_list()
        if inner_diameter[0] >= diameter:
            break
        while True:
            d, i, j = pop_rod()
            reduce_outer_diameter()
            repel_particles(i, j, d)
            update_rod_list(i, j)
            if inner_diameter[0] >= diameter or not rods:
                break


def pack_trisos(radius, fill, domain_shape='cylinder', domain_length=None,
                domain_radius=None, domain_center=[0., 0., 0.],
                n_particles=None, packing_fraction=None,
                initial_packing_fraction=0.3, contraction_rate=1/400, seed=1):
    """Generate a random, non-overlapping configuration of TRISO particles
    within a container.

    Parameters
    ----------
    radius : float
        Outer radius of TRISO particles.
    fill : openmc.Universe
        Universe which contains all layers of the TRISO particle.
    domain_shape : {'cube', 'cylinder', or 'sphere'}
        Geometry of the container in which the TRISO particles are packed.
    domain_length : float
        Length of the container (if cube or cylinder).
    domain_radius : float
        Radius of the container (if cylinder or sphere).
    domain_center : Iterable of float
        Cartesian coordinates of the center of the container.
    n_particles : int
        Number of TRISO particles to pack in the domain. Exactly one of
        'n_particles' and 'packing_fraction' should be specified -- the other
        will be calculated.
    packing_fraction : float
        Packing fraction of particles. Exactly one of 'n_particles' and
        'packing_fraction' should be specified -- the other will be calculated.
    initial_packing_fraction : float, optional
        Packing fraction used to initialize the configuration of particles in
        the domain. Default value is 0.3. It is not recommended to set the
        initial packing fraction much higher than 0.3 as the random sequential
        packing algorithm becomes prohibitively slow as it approaches its limit
        (~0.38).
    contraction_rate : float, optional
        Contraction rate of outer diameter. This can affect the speed of the
        close random packing algorithm. Default value is 1/400.
    seed : int, optional
        RNG seed.

    Returns
    -------
    trisos : list of openmc.model.TRISO
        List of TRISO particles in the domain.

    Notes
    -----
    The particle configuration is generated using a combination of random
    sequential packing (RSP) and close random packing (CRP). RSP performs
    better than CRP for lower packing fractions (pf), but it becomes
    prohibitively slow as it approaches its packing limit (~0.38). CRP can
    achieve higher pf of up to ~0.64 and scales better with increasing pf.

    If the desired pf is below some threshold for which RSP will be faster than
    CRP ('initial_packing_fraction'), only RSP is used. If a higher pf is
    required, particles with a radius smaller than the desired final radius
    (and therefore with a smaller pf) are initialized within the domain using
    RSP. This initial configuration of particles is then used as a starting
    point for CRP using Jodrey and Tory's algorithm [1]_.

    In RSP, particle centers are placed one by one at random, and placement
    attempts for a particle are made until the particle is not overlapping any
    others. This implementation of the algorithm uses a mesh over the domain
    to speed up the nearest neighbor search by only searching for a particle's
    neighbors within that mesh cell.

    In CRP, each particle is assigned two diameters, and inner and an outer,
    which approach each other during the simulation. The inner diameter,
    defined as the minimum center-to-center distance, is the true diameter of
    the particles and defines the pf. At each iteration the worst overlap
    between particles based on outer diameter is eliminated by moving the
    particles apart along the line joining their centers. Iterations continue
    until the two diameters converge or until the desired pf is reached.

    References
    ----------
    .. [1] W. S. Jodrey and E. M. Tory, "Computer simulation of close random
       packing of equal spheres", Phys. Rev. A 32 (1985) 2347-2351.

    """

    # Check for valid container geometry and dimensions
    if domain_shape not in ['cube', 'cylinder', 'sphere']:
        raise ValueError('Unable to set domain_shape to "{}". Only "cube", '
                         '"cylinder", and "sphere" are '
                         'supported."'.format(domain_shape))
    if not domain_length and domain_shape in ['cube', 'cylinder']:
        raise ValueError('"domain_length" must be specified for {} domain '
                         'geometry '.format(domain_shape))
    if not domain_radius and domain_shape in ['cylinder', 'sphere']:
        raise ValueError('"domain_radius" must be specified for {} domain '
                         'geometry '.format(domain_shape))

    if domain_shape == 'cube':
        domain = _CubicDomain(length=domain_length, particle_radius=radius,
                              center=domain_center)
    elif domain_shape == 'cylinder':
        domain = _CylindricalDomain(length=domain_length, radius=domain_radius,
                                    particle_radius=radius, center=domain_center)
    elif domain_shape == 'sphere':
        domain = _SphericalDomain(radius=domain_radius, particle_radius=radius,
                                  center=domain_center)

    # Calculate the packing fraction if the number of particles is specified;
    # otherwise, calculate the number of particles from the packing fraction.
    if ((n_particles is None and packing_fraction is None) or
        (n_particles is not None and packing_fraction is not None)):
        raise ValueError('Exactly one of "n_particles" and "packing_fraction" '
                         'must be specified.')
    elif packing_fraction is None:
        n_particles = int(n_particles)
        packing_fraction = 4/3*pi*radius**3*n_particles / domain.volume
    elif n_particles is None:
        packing_fraction = float(packing_fraction)
        n_particles = int(packing_fraction*domain.volume // (4/3*pi*radius**3))

    # Check for valid packing fractions for each algorithm
    if packing_fraction >= 0.64:
        raise ValueError('Packing fraction of {} is greater than the '
                         'packing fraction limit for close random '
                         'packing (0.64)'.format(packing_fraction))
    if initial_packing_fraction >= 0.38:
        raise ValueError('Initial packing fraction of {} is greater than the '
                         'packing fraction limit for random sequential'
                         'packing (0.38)'.format(initial_packing_fraction))
    if initial_packing_fraction > packing_fraction:
        initial_packing_fraction = packing_fraction
        if packing_fraction > 0.3:
            initial_packing_fraction = 0.3

    random.seed(seed)

    # Calculate the particle radius used in the initial random sequential
    # packing from the initial packing fraction
    initial_radius = (3/4 * initial_packing_fraction * domain.volume /
                      (pi * n_particles))**(1/3)
    domain.particle_radius = initial_radius

    # Recalculate the limits for the initial random sequential packing using
    # the desired final particle radius to ensure particles are fully contained
    # within the domain during the close random pack
    domain.limits = [[x - initial_radius + radius for x in domain.limits[0]],
                     [x + initial_radius - radius for x in domain.limits[1]]]

    # Generate non-overlapping particles for an initial inner radius using
    # random sequential packing algorithm
    particles = _random_sequential_pack(domain, n_particles)

    # Use the particle configuration produced in random sequential packing as a
    # starting point for close random pack with the desired final particle
    # radius
    if initial_packing_fraction != packing_fraction:
        domain.particle_radius = radius
        _close_random_pack(domain, particles, contraction_rate)

    trisos = []
    for p in particles:
        trisos.append(TRISO(radius, fill, p))
    return trisos
