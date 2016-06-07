import copy
from collections import Iterable
from numbers import Real

import numpy as np

import openmc
import openmc.checkvalue as cv

class TRISO(object):
    """Tristructural-isotopic (TRISO) micro fuel particle

    Parameters
    ----------
    outer_radius : int
        Outer radius of TRISO particle
    inner_univ : openmc.Universe
        Universe which contains all layers of the TRISO particle
    center : Iterable of float
        Cartesian coordinates of the center of the TRISO particle in cm

    Attributes
    ----------
    cell : opemc.Cell
        Cell which contains the TRISO universe
    center : numpy.ndarray
        Cartesian coordinates of the center of the TRISO particle in cm
    outside : openmc.Region
        Region of space outside of the TRISO particle
    bounding_box : tuple of numpy.ndarray
        Lower-left and upper-right coordinates of an axis-aligned bounding box
        for the TRISO particle

    """

    def __init__(self, outer_radius, inner_univ, center=(0., 0., 0.)):
        self._surface = openmc.Sphere(R=outer_radius)
        self._cell = openmc.Cell(fill=inner_univ, region=-self._surface)
        self.center = np.asarray(center)

    @property
    def bounding_box(self):
        return self.cell.region.bounding_box

    @property
    def cell(self):
        return self._cell

    @property
    def center(self):
        return self._center

    @property
    def outside(self):
        return +self._surface

    @center.setter
    def center(self, center):
        cv.check_type('TRISO center', center, Iterable, Real)
        self._surface.x0 = center[0]
        self._surface.y0 = center[1]
        self._surface.z0 = center[2]
        self.cell.translation = center
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

        ll, ur = self.bounding_box
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
                t_copy.cell.id = None
                t_copy.cell.fill = t.cell.fill
                t_copy._surface.id = None
                triso_locations[idx].append(t_copy)

    # Create universes
    universes = np.empty(shape[::-1], dtype=openmc.Universe)
    for idx, triso_list in sorted(triso_locations.items()):
        if len(triso_list) > 0:
            outside_trisos = openmc.Intersection(*[t.outside for t in triso_list])
            background_cell = openmc.Cell(fill=background, region=outside_trisos)
        else:
            background_cell = openmc.Cell(fill=background)

        u = openmc.Universe()
        u.add_cell(background_cell)
        for t in triso_list:
            u.add_cell(t.cell)
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
