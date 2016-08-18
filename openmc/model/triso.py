from __future__ import division
import copy
from collections import Iterable, defaultdict
from numbers import Real
import warnings
import itertools
import random
from random import uniform, gauss
from heapq import heappush, heappop
from math import pi, sin, cos, floor, log10, sqrt

import numpy as np
import scipy.spatial
from scipy.spatial.distance import cdist

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
        super(TRISO, self).__init__(fill=fill, region=-self._surface)
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
            outside_trisos = openmc.Intersection(*[~t.region for t in triso_list])
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


def pack_trisos(radius, fill, domain_shape='cylinder', domain_length=None,
                domain_radius=None, domain_center=(0., 0., 0.),
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
    sequential packing (RSP) and close random packing (CRP). RSP is faster than
    CRP for lower packing fractions (pf), but it becomes prohibitively slow as
    it approaches its packing limit (~0.38). CRP can achieve higher pf of up to
    ~0.64 and scales better with increasing pf.

    If the desired pf is below some threshold for which RSP performs better
    than CRP ('initial_packing_fraction'), only RSP is used. If a higher pf is
    required, particles with a radius smaller than the desired final radius
    (and therefore with a smaller pf) are initialized within the domain using
    RSP. This initial configuration of particles is then used as a starting
    point for CRP using Jodrey and Tory's algorithm [1]_.

    In RSP, particle centers are placed one by one at random, and placement
    attempts for a particle are made until the particle is not overlapping any
    others. This implementation of the algorithm uses a lattice over the domain
    to speed up the nearest neighbor search by only searching for a particle's
    neighbors within that lattice cell.

    In CRP, each particle is assigned two diameters, and inner and an outer,
    which approach each other during the simulation. The inner diameter,
    defined as the minimum center-to-center distance, is the true diameter of
    the particles and defines the pf. At each iteration the worst overlap
    between particles based on outer diameter is eliminated by moving the
    particles apart along the line joining their centers and the outer diameter
    is decreased. Iterations continue until the two diameters converge or until
    the desired pf is reached.

    References
    ----------
    .. [1] W. S. Jodrey and E. M. Tory, "Computer simulation of close random
       packing of equal spheres", Phys. Rev. A 32 (1985) 2347-2351.

    """

    def get_domain_volume():
        """Calculates the volume of the container in which the TRISO particles
        are packed.

        Returns
        -------
        float
            Volume of the domain.

        """

        if domain_shape is 'cube':
            return domain_length**3
        elif domain_shape is 'cylinder':
            return domain_length * pi * domain_radius**2
        elif domain_shape is 'sphere':
            return 4/3 * pi * domain_radius**3


    def get_cell_length(radius):
        """Calculates the length of a lattice element in x-, y-, and
        z-directions.

        Parameters
        ----------
        radius : float
            Radius of the particle.

        Returns
        -------
        tuple of float
            Length of lattice cell in x-, y-, and z-directions.

        """

        if domain_length:
            m = domain_length/int(domain_length/(4*radius))
        if domain_radius:
            n = 2*domain_radius/int(domain_radius/(2*radius))

        if domain_shape is 'cube':
            return (m, m, m)
        elif domain_shape is 'cylinder':
            return (n, n, m)
        elif domain_shape is 'sphere':
            return (n, n, n)


    def get_boundary_extremes():
        """Calculates the minimum and maximum positions in x-, y-, and
        z-directions where a particle center can be placed within the domain.

        Returns
        -------
        llim, ulim : tuple of float
            Minimum and maximum position in x-, y-, and z-directions where
            particle center can be placed.

        """

        if domain_length:
            x_min = radius
            x_max = domain_length - radius
        if domain_radius:
            r_min = radius - domain_radius
            r_max = domain_radius - radius

        if domain_shape is 'cube':
            return (x_min, x_min, x_min), (x_max, x_max, x_max)
        elif domain_shape is 'cylinder':
            return (r_min, r_min, x_min), (r_max, r_max, x_max)
        elif domain_shape is 'sphere':
            return (r_min, r_min, r_min), (r_max, r_max, r_max)


    def get_particle_offset():
        """Calculates the offset in x-, y-, and z-directions of the particle
        center based on the domain center

        Returns
        -------
        tuple of float
            Amount to offset particle center in x-, y-, and z-directions

        """

        if domain_shape is 'cube':
            return np.array(domain_center) - domain_length/2
        elif domain_shape is 'cylinder':
            return np.array(domain_center) - (0, 0, domain_length/2)
        elif domain_shape is 'sphere':
            return np.array(domain_center)


    def inner_packing_fraction():
        """Calculates the true packing fraction of the particles based on the
        inner diameter.

        Returns
        -------
        float
            Packing fraction calculated from inner diameter.

        """

        return (4/3 * pi * (inner_diameter[0]/2)**3 * n_particles /
                domain_volume)


    def outer_packing_fraction():
        """Calculates the nominal packing fraction of the particles based on
        the outer diameter.

        Returns
        -------
        float
            Packing fraction calculated from outer diameter.

        """

        return (4/3 * pi * (outer_diameter[0]/2)**3 * n_particles /
                domain_volume)


    def random_point_cube():
        """Generate Cartesian coordinates of center of a particle that is
        contained entirely within cubic domain with uniform probability.

        Returns
        -------
        list of float
            Cartesian coordinates of particle center.

        """

        return [uniform(llim[0], ulim[0]),
                uniform(llim[0], ulim[0]),
                uniform(llim[0], ulim[0])]


    def random_point_cylinder():
        """Generate Cartesian coordinates of center of a particle that is
        contained entirely within cylindrical domain with uniform probability
        (see http://mathworld.wolfram.com/DiskPointPicking.html for generating
        random points on a disk).

        Returns
        -------
        list of float
            Cartesian coordinates of particle center.

        """

        r = sqrt(uniform(0, ulim[0]**2))
        t = uniform(0, 2*pi)
        return [r*cos(t), r*sin(t), uniform(llim[2], ulim[2])]


    def random_point_sphere():
        """Generate Cartesian coordinates of center of a particle that is
        contained entirely within spherical domain with uniform probability.

        Returns
        -------
        list of float
            Cartesian coordinates of particle center.

        """

        x = (gauss(0, 1), gauss(0, 1), gauss(0, 1))
        r = (uniform(0, ulim[0]**3)**(1/3) / sqrt(x[0]**2 + x[1]**2 + x[2]**2))
        return [r*i for i in x]


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
            rod[1] = None
            rod[2] = None


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
            if i is not None and j is not None:
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
        a = np.dstack(([i for i in range(len(n))], n, d))[0]

        # Array of nearest neighbor indices, indices of particles they are
        # nearest neighbors of, and distances between them
        b = a[a[:,1].argsort()]
        b[:,[0, 1]] = b[:,[1, 0]]

        # Find the intersection between 'a' and 'b': a list of particles who
        # are each other's nearest neighbors and the distance between them
        r = [x for x in {tuple(x) for x in a} & {tuple(x) for x in b}]

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


    def reduce_outer_diameter():
        """Reduce the outer diameter so that at the (i+1)-st iteration it is:

            d_out^(i+1) = d_out^(i) - (1/2)^(j) * d_out0 * k / n,

        where k is the contraction rate, n is the number of particles, and

            j = floor(-log10(pf_out - pf_in)).

        """

        j = floor(-log10(outer_packing_fraction() - inner_packing_fraction()))
        outer_diameter[0] = (outer_diameter[0] - 0.5**j *
                             initial_outer_diameter * contraction_rate /
                             n_particles)


    def update_mesh(i):
        """Update which lattice cells the particle is in based on new particle
        center coordinates.

        'mesh'/'mesh_map' is a two way dictionary used to look up which
        particles are located within one diameter of a given lattice cell and
        which lattice cells a given particle center is within one diameter of.
        This is used to speed up the nearest neighbor search.

        Parameters
        ----------
        i : int
            Index of particle in particles array.

        """

        # Determine which lattice cells the particle is in and remove the
        # particle id from those cells
        for idx in mesh_map[i]:
            mesh[idx].remove(i)
        del mesh_map[i]

        # Determine which lattice cells are within one diameter of particle's
        # center and add this particle to the list of particles in those cells
        for idx in cell_list(particles[i], diameter):
            mesh[idx].add(i)
            mesh_map[i].add(idx)


    def apply_boundary_conditions(i, j):
        """Apply reflective boundary conditions to particles i and j.

        Parameters
        ----------
        i, j : int
            Index of particles in particles array.

        """

        for k in range(3):
            if particles[i][k] < llim[k]:
                particles[i][k] = llim[k]
            elif particles[i][k] > ulim[k]:
                particles[i][k] = ulim[k]
            if particles[j][k] < llim[k]:
                particles[j][k] = llim[k]
            elif particles[j][k] > ulim[k]:
                particles[j][k] = ulim[k]


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
        apply_boundary_conditions(i, j)

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
        idx = list(mesh[cell_index(particles[i])])
        dists = cdist([particles[i]], particles[idx])[0]
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


    def cell_index_cube(p, cl=None):
        """Calculate the index of the lattice cell in which the given particle
        center falls.

        Parameters
        ----------
        p : list of float
            Cartesian coordinates of particle center.
        cl : list of float
            Length of the lattice cells in x-, y-, and z-directions.

        Returns
        -------
        tuple of int
            Indices of lattice cell.

        """

        if cl is None:
            cl = cell_length

        return tuple(int(p[i]/cl[i]) for i in range(3))


    def cell_index_cylinder(p, cl=None):
        """Calculate the index of the lattice cell in which the given particle
        center falls.

        Parameters
        ----------
        p : list of float
            Cartesian coordinates of particle center.
        cl : list of float
            Length of the lattice cells in x-, y-, and z-directions.

        Returns
        -------
        tuple of int
            Indices of lattice cell.

        """

        if cl is None:
            cl = cell_length

        return (int((p[0] + domain_radius)/cl[0]),
                int((p[1] + domain_radius)/cl[1]), int(p[2]/cl[2]))


    def cell_index_sphere(p, cl=None):
        """Calculate the index of the lattice cell in which the given particle
        center falls.

        Parameters
        ----------
        p : list of float
            Cartesian coordinates of particle center.
        cl : list of float
            Length of the lattice cells in x-, y-, and z-directions.

        Returns
        -------
        tuple of int
            Indices of lattice cell.

        """

        if cl is None:
            cl = cell_length

        return tuple(int((p[i] + domain_radius)/cl[i]) for i in range(3))


    def cell_list_cube(p, d, cl=None):
        """Return the indices of all cells within the given distance of the
        point.

        Parameters
        ----------
        p : list of float
            Cartesian coordinates of particle center.
        d : float
            Find all lattice cells that are within a radius of length 'd' of
            the particle center.
        cl : list of float
            Length of the lattice cells in x-, y-, and z-directions.

        Returns
        -------
        list of tuple of int
            Indices of lattice cells.

        """

        if cl is None:
            cl = cell_length

        r = [[a/cl[i] for a in [p[i]-d, p[i], p[i]+d] if a > 0 and
              a < domain_length] for i in range(3)]

        return list(itertools.product(*({int(i) for i in j} for j in r)))


    def cell_list_cylinder(p, d, cl=None):
        """Return the indices of all cells within the given distance of the
        point.

        Parameters
        ----------
        p : list of float
            Cartesian coordinates of particle center.
        d : float
            Find all lattice cells that are within a radius of length 'd' of
            the particle center.
        cl : list of float
            Length of the lattice cells in x-, y-, and z-directions.

        Returns
        -------
        list of tuple of int
            Indices of lattice cells.

        """

        if cl is None:
            cl = cell_length

        x, y = [[(a + domain_radius)/cl[i] for a in [p[i]-d, p[i], p[i]+d]
                 if a > -domain_radius and a < domain_radius] for i in range(2)]

        z = [a/cl[2] for a in [p[2]-d, p[2], p[2]+d] if a > 0
             and a < domain_length]

        return list(itertools.product(*({int(i) for i in j} for j in (x, y, z))))


    def cell_list_sphere(p, d, cl=None):
        """Return the indices of all cells within the given distance of the
        point.

        Parameters
        ----------
        p : list of float
            Cartesian coordinates of particle center.
        d : float
            Find all lattice cells that are within a radius of length 'd' of
            the particle center.
        cl : list of float
            Length of the lattice cells in x-, y-, and z-directions.

        Returns
        -------
        list of tuple of int
            Indices of lattice cells.

        """

        if cl is None:
            cl = cell_length

        r = [[(a + domain_radius)/cl[i] for a in [p[i]-d, p[i], p[i]+d]
              if a > -domain_radius and a < domain_radius] for i in range(3)]

        return list(itertools.product(*({int(i) for i in j} for j in r)))


    def random_sequential_pack():
        """Random sequential packing of particles whose radius is determined by
        initial packing fraction.

        Returns
        ------
        numpy.ndarray
            Cartesian coordinates of centers of TRISO particles.

        """

        # Set parameters for initial random sequential packing of particles.
        r = (3/4*initial_packing_fraction*domain_volume/(pi*n_particles))**(1/3)
        d = 2*r
        sqd = d**2
        cl = get_cell_length(r)

        particles = []
        mesh = defaultdict(list)

        for i in range(n_particles):
            # Randomly sample new center coordinates while there are any overlaps
            while True:
                p = random_point()
                idx = cell_index(p, cl)
                if any((p[0]-q[0])**2 + (p[1]-q[1])**2 + (p[2]-q[2])**2 < sqd
                       for q in mesh[idx]):
                    continue
                else:
                    break
            particles.append(p)

            for idx in cell_list(p, d, cl):
                mesh[idx].append(p)

        return np.array(particles)


    def close_random_pack():
        """Close random packing of particles using the Jodrey-Tory algorithm.

        """

        for i in range(n_particles):
            for idx in cell_list(particles[i], diameter):
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

    domain_volume = get_domain_volume()
    llim, ulim = get_boundary_extremes()
    offset = get_particle_offset()

    # Calculate the packing fraction if the number of particles is specified;
    # otherwise, calculate the number of particles from the packing fraction.
    if ((n_particles is None and packing_fraction is None) or
        (n_particles is not None and packing_fraction is not None)):
        raise ValueError('Exactly one of "n_particles" and "packing_fraction" '
                         'must be specified.')
    elif packing_fraction is None:
        packing_fraction = 4/3*pi*radius**3*n_particles / domain_volume
    elif n_particles is None:
        n_particles = int(packing_fraction*domain_volume // (4/3*pi*radius**3))

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

    # Set domain dependent functions
    if domain_shape is 'cube':
        random_point = random_point_cube
        cell_list = cell_list_cube
        cell_index = cell_index_cube
    elif domain_shape is 'cylinder':
        random_point = random_point_cylinder
        cell_list = cell_list_cylinder
        cell_index = cell_index_cylinder
    elif domain_shape is 'sphere':
        random_point = random_point_sphere
        cell_list = cell_list_sphere
        cell_index = cell_index_sphere

    random.seed(seed)

    # Generate non-overlapping particles for an initial inner radius using
    # random sequential packing algorithm
    particles = random_sequential_pack()

    # Use the particle configuration produced in random sequential packing as a
    # starting point for close random pack with the desired final particle radius
    if initial_packing_fraction != packing_fraction:
        diameter = 2*radius
        cell_length = get_cell_length(radius)

        # Outer diameter initially set to arbitrary value that yields pf of 1
        initial_outer_diameter = 2*(domain_volume/(n_particles*4/3*pi))**(1/3)

        # Inner and outer diameter of particles will change during packing
        outer_diameter = [initial_outer_diameter]
        inner_diameter = [0]

        rods = []
        rods_map = {}
        mesh = defaultdict(set)
        mesh_map = defaultdict(set)

        close_random_pack()

    trisos = []
    for p in particles:
        trisos.append(TRISO(radius, fill, p + offset))

    return trisos
