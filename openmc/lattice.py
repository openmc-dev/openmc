from abc import ABCMeta
from collections import OrderedDict
from collections.abc import Iterable
from copy import deepcopy
from math import sqrt, floor
from numbers import Real
from xml.etree import ElementTree as ET

import numpy as np
import warnings
import types

import openmc.checkvalue as cv
import openmc
from openmc._xml import get_text
from openmc.mixin import IDManagerMixin


class Lattice(IDManagerMixin, metaclass=ABCMeta):
    """A repeating structure wherein each element is a universe.

    Parameters
    ----------
    lattice_id : int, optional
        Unique identifier for the lattice. If not specified, an identifier will
        automatically be assigned.
    name : str, optional
        Name of the lattice. If not specified, the name is the empty string.

    Attributes
    ----------
    id : int
        Unique identifier for the lattice
    name : str
        Name of the lattice
    pitch : Iterable of float
        Pitch of the lattice in each direction in cm
    outer : openmc.Universe
        A universe to fill all space outside the lattice
    universes : Iterable of Iterable of openmc.Universe
        A two-or three-dimensional list/array of universes filling each element
        of the lattice

    """

    next_id = 1
    used_ids = openmc.Universe.used_ids

    def __init__(self, lattice_id=None, name=''):
        # Initialize Lattice class attributes
        self.id = lattice_id
        self.name = name
        self._pitch = None
        self._outer = None
        self._universes = None

    @property
    def name(self):
        return self._name

    @property
    def pitch(self):
        return self._pitch

    @property
    def outer(self):
        return self._outer

    @property
    def universes(self):
        return self._universes

    @name.setter
    def name(self, name):
        if name is not None:
            cv.check_type('lattice name', name, str)
            self._name = name
        else:
            self._name = ''

    @outer.setter
    def outer(self, outer):
        cv.check_type('outer universe', outer, openmc.Universe)
        self._outer = outer

    @staticmethod
    def from_hdf5(group, universes):
        """Create lattice from HDF5 group

        Parameters
        ----------
        group : h5py.Group
            Group in HDF5 file
        universes : dict
            Dictionary mapping universe IDs to instances of
            :class:`openmc.Universe`.

        Returns
        -------
        openmc.Lattice
            Instance of lattice subclass

        """
        lattice_type = group['type'][()].decode()
        if lattice_type == 'rectangular':
            return openmc.RectLattice.from_hdf5(group, universes)
        elif lattice_type == 'hexagonal':
            return openmc.HexLattice.from_hdf5(group, universes)
        else:
            raise ValueError('Unkown lattice type: {}'.format(lattice_type))

    def get_unique_universes(self):
        """Determine all unique universes in the lattice

        Returns
        -------
        universes : collections.OrderedDict
            Dictionary whose keys are universe IDs and values are
            :class:`openmc.Universe` instances

        """

        univs = OrderedDict()
        for k in range(len(self._universes)):
            for j in range(len(self._universes[k])):
                if isinstance(self._universes[k][j], openmc.Universe):
                    u = self._universes[k][j]
                    univs[u._id] = u
                else:
                    for i in range(len(self._universes[k][j])):
                        u = self._universes[k][j][i]
                        assert isinstance(u, openmc.Universe)
                        univs[u._id] = u

        if self.outer is not None:
            univs[self.outer._id] = self.outer

        return univs

    def get_nuclides(self):
        """Returns all nuclides in the lattice

        Returns
        -------
        nuclides : list of str
            List of nuclide names

        """

        nuclides = []

        # Get all unique Universes contained in each of the lattice cells
        unique_universes = self.get_unique_universes()

        # Append all Universes containing each cell to the dictionary
        for universe in unique_universes.values():
            for nuclide in universe.get_nuclides():
                if nuclide not in nuclides:
                    nuclides.append(nuclide)

        return nuclides

    def get_all_cells(self, memo=None):
        """Return all cells that are contained within the lattice

        Returns
        -------
        cells : collections.OrderedDict
            Dictionary whose keys are cell IDs and values are :class:`Cell`
            instances

        """
        cells = OrderedDict()

        if memo and self in memo:
            return cells

        if memo is not None:
            memo.add(self)

        unique_universes = self.get_unique_universes()

        for universe in unique_universes.values():
            cells.update(universe.get_all_cells(memo))

        return cells

    def get_all_materials(self, memo=None):
        """Return all materials that are contained within the lattice

        Returns
        -------
        materials : collections.OrderedDict
            Dictionary whose keys are material IDs and values are
            :class:`Material` instances

        """

        materials = OrderedDict()

        # Append all Cells in each Cell in the Universe to the dictionary
        cells = self.get_all_cells(memo)
        for cell in cells.values():
            materials.update(cell.get_all_materials(memo))

        return materials

    def get_all_universes(self):
        """Return all universes that are contained within the lattice

        Returns
        -------
        universes : collections.OrderedDict
            Dictionary whose keys are universe IDs and values are
            :class:`Universe` instances

        """

        # Initialize a dictionary of all Universes contained by the Lattice
        # in each nested Universe level
        all_universes = OrderedDict()

        # Get all unique Universes contained in each of the lattice cells
        unique_universes = self.get_unique_universes()

        # Add the unique Universes filling each Lattice cell
        all_universes.update(unique_universes)

        # Append all Universes containing each cell to the dictionary
        for universe in unique_universes.values():
            all_universes.update(universe.get_all_universes())

        return all_universes

    def get_universe(self, idx):
        r"""Return universe corresponding to a lattice element index

        Parameters
        ----------
        idx : Iterable of int
            Lattice element indices. For a rectangular lattice, the indices are
            given in the :math:`(x,y)` or :math:`(x,y,z)` coordinate system. For
            hexagonal lattices, they are given in the :math:`x,\alpha` or
            :math:`x,\alpha,z` coordinate systems for "y" orientations and
            :math:`\alpha,y` or :math:`\alpha,y,z` coordinate systems for "x"
            orientations.

        Returns
        -------
        openmc.Universe
            Universe with given indices

        """
        idx_u = self.get_universe_index(idx)
        if self.ndim == 2:
            return self.universes[idx_u[0]][idx_u[1]]
        else:
            return self.universes[idx_u[0]][idx_u[1]][idx_u[2]]

    def find(self, point):
        """Find cells/universes/lattices which contain a given point

        Parameters
        ----------
        point : 3-tuple of float
            Cartesian coordinates of the point

        Returns
        -------
        list
            Sequence of universes, cells, and lattices which are traversed to
            find the given point

        """
        idx, p = self.find_element(point)
        if self.is_valid_index(idx):
            u = self.get_universe(idx)
        else:
            if self.outer is not None:
                u = self.outer
            else:
                return []
        return [(self, idx)] + u.find(p)

    def clone(self, clone_materials=True, clone_regions=True, memo=None):
        """Create a copy of this lattice with a new unique ID, and clones
        all universes within this lattice.

        Parameters
        ----------
        clone_materials : bool
            Whether to create separate copies of the materials filling cells
            contained in this lattice and its outer universe.
        clone_regions : bool
            Whether to create separate copies of the regions bounding cells
            contained in this lattice and its outer universe.
        memo : dict or None
            A nested dictionary of previously cloned objects. This parameter
            is used internally and should not be specified by the user.

        Returns
        -------
        clone : openmc.Lattice
            The clone of this lattice

        """

        if memo is None:
            memo = {}

        # If no memoize'd clone exists, instantiate one
        if self not in memo:
            clone = deepcopy(self)
            clone.id = None

            if self.outer is not None:
                clone.outer = self.outer.clone(clone_materials, clone_regions,
                     memo)

            # Assign universe clones to the lattice clone
            for i in self.indices:
                if isinstance(self, RectLattice):
                    clone.universes[i] = self.universes[i].clone(
                         clone_materials, clone_regions, memo)
                else:
                    if self.ndim == 2:
                        clone.universes[i[0]][i[1]] = \
                            self.universes[i[0]][i[1]].clone(clone_materials,
                                 clone_regions, memo)
                    else:
                        clone.universes[i[0]][i[1]][i[2]] = \
                            self.universes[i[0]][i[1]][i[2]].clone(
                            clone_materials, clone_regions, memo)

            # Memoize the clone
            memo[self] = clone

        return memo[self]


class RectLattice(Lattice):
    """A lattice consisting of rectangular prisms.

    To completely define a rectangular lattice, the
    :attr:`RectLattice.lower_left` :attr:`RectLattice.pitch`,
    :attr:`RectLattice.outer`, and :attr:`RectLattice.universes` properties need
    to be set.

    Most methods for this class use a natural indexing scheme wherein elements
    are assigned an index corresponding to their position relative to the
    (x,y,z) axes in a Cartesian coordinate system, i.e., an index of (0,0,0) in
    the lattice gives the element whose x, y, and z coordinates are the
    smallest. However, note that when universes are assigned to lattice elements
    using the :attr:`RectLattice.universes` property, the array indices do not
    correspond to natural indices.

    Parameters
    ----------
    lattice_id : int, optional
        Unique identifier for the lattice. If not specified, an identifier will
        automatically be assigned.
    name : str, optional
        Name of the lattice. If not specified, the name is the empty string.

    Attributes
    ----------
    id : int
        Unique identifier for the lattice
    name : str
        Name of the lattice
    pitch : Iterable of float
        Pitch of the lattice in the x, y, and (if applicable) z directions in
        cm.
    outer : openmc.Universe
        A universe to fill all space outside the lattice
    universes : Iterable of Iterable of openmc.Universe
        A two- or three-dimensional list/array of universes filling each element
        of the lattice. The first dimension corresponds to the z-direction (if
        applicable), the second dimension corresponds to the y-direction, and
        the third dimension corresponds to the x-direction. Note that for the
        y-direction, a higher index corresponds to a lower physical
        y-value. Each z-slice in the array can be thought of as a top-down view
        of the lattice.
    lower_left : Iterable of float
        The Cartesian coordinates of the lower-left corner of the lattice. If
        the lattice is two-dimensional, only the x- and y-coordinates are
        specified.
    indices : list of tuple
        A list of all possible (z,y,x) or (y,x) lattice element indices. These
        indices correspond to indices in the :attr:`RectLattice.universes`
        property.
    ndim : int
        The number of dimensions of the lattice
    shape : Iterable of int
        An array of two or three integers representing the number of lattice
        cells in the x- and y- (and z-) directions, respectively.

    """

    def __init__(self, lattice_id=None, name=''):
        super().__init__(lattice_id, name)

        # Initialize Lattice class attributes
        self._lower_left = None

    def __repr__(self):
        string = 'RectLattice\n'
        string += '{: <16}=\t{}\n'.format('\tID', self._id)
        string += '{: <16}=\t{}\n'.format('\tName', self._name)
        string += '{: <16}=\t{}\n'.format('\tShape', self.shape)
        string += '{: <16}=\t{}\n'.format('\tLower Left', self._lower_left)
        string += '{: <16}=\t{}\n'.format('\tPitch', self._pitch)
        string += '{: <16}=\t{}\n'.format(
            '\tOuter', self._outer._id if self._outer is not None else None)

        string += '{: <16}\n'.format('\tUniverses')

        # Lattice nested Universe IDs
        for i, universe in enumerate(np.ravel(self._universes)):
            string += '{} '.format(universe._id)

            # Add a newline character every time we reach end of row of cells
            if (i + 1) % self.shape[0] == 0:
                string += '\n'

        string = string.rstrip('\n')

        return string

    @property
    def indices(self):
        if self.ndim == 2:
            return list(np.broadcast(*np.ogrid[
                :self.shape[1], :self.shape[0]]))
        else:
            return list(np.broadcast(*np.ogrid[
                :self.shape[2], :self.shape[1], :self.shape[0]]))

    @property
    def _natural_indices(self):
        """Iterate over all possible (x,y) or (x,y,z) lattice element indices.

        This property is used when constructing distributed cell and material
        paths. Most importantly, the iteration order matches that used on the
        Fortran side.

        """
        if self.ndim == 2:
            nx, ny = self.shape
            for iy in range(ny):
                for ix in range(nx):
                    yield (ix, iy)
        else:
            nx, ny, nz = self.shape
            for iz in range(nz):
                for iy in range(ny):
                    for ix in range(nx):
                        yield (ix, iy, iz)

    @property
    def lower_left(self):
        return self._lower_left

    @property
    def ndim(self):
        if self.pitch is not None:
            return len(self.pitch)
        else:
            raise ValueError('Number of dimensions cannot be determined until '
                             'the lattice pitch has been set.')

    @property
    def shape(self):
        return self._universes.shape[::-1]

    @lower_left.setter
    def lower_left(self, lower_left):
        cv.check_type('lattice lower left corner', lower_left, Iterable, Real)
        cv.check_length('lattice lower left corner', lower_left, 2, 3)
        self._lower_left = lower_left

    @Lattice.pitch.setter
    def pitch(self, pitch):
        cv.check_type('lattice pitch', pitch, Iterable, Real)
        cv.check_length('lattice pitch', pitch, 2, 3)
        for dim in pitch:
            cv.check_greater_than('lattice pitch', dim, 0.0)
        self._pitch = pitch

    @Lattice.universes.setter
    def universes(self, universes):
        cv.check_iterable_type('lattice universes', universes, openmc.Universe,
                               min_depth=2, max_depth=3)
        self._universes = np.asarray(universes)

    def find_element(self, point):
        """Determine index of lattice element and local coordinates for a point

        Parameters
        ----------
        point : Iterable of float
            Cartesian coordinates of point

        Returns
        -------
        2- or 3-tuple of int
            A tuple of the corresponding (x,y,z) lattice element indices
        3-tuple of float
            Carestian coordinates of the point in the corresponding lattice
            element coordinate system

        """
        ix = floor((point[0] - self.lower_left[0])/self.pitch[0])
        iy = floor((point[1] - self.lower_left[1])/self.pitch[1])
        if self.ndim == 2:
            idx = (ix, iy)
        else:
            iz = floor((point[2] - self.lower_left[2])/self.pitch[2])
            idx = (ix, iy, iz)
        return idx, self.get_local_coordinates(point, idx)

    def get_local_coordinates(self, point, idx):
        """Determine local coordinates of a point within a lattice element

        Parameters
        ----------
        point : Iterable of float
            Cartesian coordinates of point
        idx : Iterable of int
            (x,y,z) indices of lattice element. If the lattice is 2D, the z
            index can be omitted.

        Returns
        -------
        3-tuple of float
            Cartesian coordinates of point in the lattice element coordinate
            system

        """
        x = point[0] - (self.lower_left[0] + (idx[0] + 0.5)*self.pitch[0])
        y = point[1] - (self.lower_left[1] + (idx[1] + 0.5)*self.pitch[1])
        if self.ndim == 2:
            z = point[2]
        else:
            z = point[2] - (self.lower_left[2] + (idx[2] + 0.5)*self.pitch[2])
        return (x, y, z)

    def get_universe_index(self, idx):
        """Return index in the universes array corresponding
        to a lattice element index

        Parameters
        ----------
        idx : Iterable of int
            Lattice element indices in the :math:`(x,y,z)` coordinate system

        Returns
        -------
        2- or 3-tuple of int
            Indices used when setting the :attr:`RectLattice.universes` property

        """
        max_y = self.shape[1] - 1
        if self.ndim == 2:
            x, y = idx
            return (max_y - y, x)
        else:
            x, y, z = idx
            return (z, max_y - y, x)

    def is_valid_index(self, idx):
        """Determine whether lattice element index is within defined range

        Parameters
        ----------
        idx : Iterable of int
            Lattice element indices in the :math:`(x,y,z)` coordinate system

        Returns
        -------
        bool
            Whether index is valid

        """
        if self.ndim == 2:
            return (0 <= idx[0] < self.shape[0] and
                    0 <= idx[1] < self.shape[1])
        else:
            return (0 <= idx[0] < self.shape[0] and
                    0 <= idx[1] < self.shape[1] and
                    0 <= idx[2] < self.shape[2])

    def discretize(self, strategy="degenerate",
                   universes_to_ignore=[],
                   materials_to_clone=[],
                   lattice_neighbors=[], key=lambda univ: univ.id):
        """Discretize the lattice with either a degenerate or a local neighbor
        symmetry strategy

        'Degenerate' clones every universe in the lattice, thus making them all
        uniquely defined. This is typically required if depletion or thermal
        hydraulics will make every universe's environment unique.

        'Local neighbor symmetry' groups universes with similar neighborhoods.
        These clusters of cells and materials provide increased convergence
        speed to multi-group cross sections tallies. The user can specify
        the lattice's neighbors to discriminate between two sides of a
        lattice for example.

        Parameters
        ----------
        strategy : {'degenerate', 'lns'}
            Which strategy to adopt when discretizing the lattice
        universes_to_ignore : Iterable of Universe
            Lattice universes that need not be discretized
        materials_to_clone : Iterable of Material
            List of materials that should be cloned when discretizing
        lattice_neighbors : Iterable of Universe
            List of the lattice's neighbors. By default, if present, the
            lattice outer universe will be used. The neighbors should be
            ordered as follows [top left, top, top right, left, right,
            bottom left, bottom, bottom right]
        key : function
            Function of argument a universe that is used to extract a
            comparison key. This function will be called on each universe's
            neighbors in the lattice to form a neighbor pattern. This pattern
            is then used to identify unique neighbor symmetries.
        """

        # Check routine inputs
        if self.ndim != 2:
            raise NotImplementedError("LNS discretization is not implemented "
                                      "for 1D and 3D lattices")

        cv.check_value('strategy', strategy, ('degenerate', 'lns'))
        cv.check_type('universes_to_ignore', universes_to_ignore, Iterable,
                      openmc.Universe)
        cv.check_type('materials_to_clone', materials_to_clone, Iterable,
                      openmc.Material)
        cv.check_type('lattice_neighbors', lattice_neighbors, Iterable,
                      openmc.Universe)
        cv.check_value('number of lattice_neighbors', len(lattice_neighbors),
                       (0, 8))
        cv.check_type('key', key, types.FunctionType)

        # Use outer universe if neighbors are missing and outer is defined
        if self.outer is not None and len(lattice_neighbors) == 0:
            lattice_neighbors = [key(self.outer) for i in range(8)]
        elif len(lattice_neighbors) == 8:
            lattice_neighbors = [key(universe) for universe in
                                 lattice_neighbors]

        # Dictionary that will keep track of where each pattern appears, how
        # it was rotated and/or symmetrized
        patterns = {}

        # Initialize pattern array
        pattern = np.empty(shape=(3, 3), dtype=type(key(self.universes[0][0])))

        # Define an auxiliary function that returns a universe's neighbors
        # that are outside the lattice
        def find_edge_neighbors(pattern, i, j):

            # If no neighbors have been specified, start with an empty array
            if len(lattice_neighbors) == 0:
                return

            # Left edge
            if i == 0:
                pattern[:, 0] = lattice_neighbors[3]
                if j == 0:
                    pattern[0, 0] = lattice_neighbors[0]
                elif j == self.shape[1] - 1:
                    pattern[2, 0] = lattice_neighbors[5]

            # Bottom edge
            if j == 0:
                pattern[0, 1] = lattice_neighbors[1]
                if i != 0:
                    pattern[0, 0] = lattice_neighbors[1]
                if i != self.shape[0] - 1:
                    pattern[0, 2] = lattice_neighbors[1]

            # Right edge
            if i == self.shape[0] - 1:
                pattern[:, 2] = lattice_neighbors[4]
                if j == 0:
                    pattern[0, 2] = lattice_neighbors[2]
                elif j == self.shape[1] - 1:
                    pattern[2, 2] = lattice_neighbors[7]

            # Top edge
            if j == self.shape[1] - 1:
                pattern[2, 1] = lattice_neighbors[6]
                if i != 0:
                    pattern[2, 0] = lattice_neighbors[6]
                if i != self.shape[0] - 1:
                    pattern[2, 2] = lattice_neighbors[6]

        # Define an auxiliary function that returns a universe's neighbors
        # among the universes inside the lattice
        def find_lattice_neighbors(pattern, i, j):

            # Away from left edge
            if i != 0:
                if j > 0:
                    pattern[0, 0] = key(self.universes[j-1][i-1])
                pattern[1, 0] = key(self.universes[j][i-1])
                if j < self.shape[1] - 1:
                    pattern[2, 0] = key(self.universes[j+1][i-1])

            # Away from bottom edge
            if j != 0:
                if i > 0:
                    pattern[0, 0] = key(self.universes[j-1][i-1])
                pattern[0, 1] = key(self.universes[j-1][i])
                if i < self.shape[0] - 1:
                    pattern[0, 2] = key(self.universes[j-1][i+1])

            # Away from right edge
            if i != self.shape[0] - 1:
                if j > 0:
                    pattern[0, 2] = key(self.universes[j-1][i+1])
                pattern[1, 2] = key(self.universes[j][i+1])
                if j < self.shape[1] - 1:
                    pattern[2, 2] = key(self.universes[j+1][i+1])

            # Away from top edge
            if j != self.shape[1] - 1:
                if i > 0:
                    pattern[2, 0] = key(self.universes[j+1][i-1])
                pattern[2, 1] = key(self.universes[j+1][i])
                if i < self.shape[0] - 1:
                    pattern[2, 2] = key(self.universes[j+1][i+1])

        # Analyze lattice, find unique patterns in groups of universes
        for j in range(self.shape[1]):
            for i in range(self.shape[0]):

                # Skip universes to ignore
                if self.universes[j][i] in universes_to_ignore:
                    continue

                # Create a neighborhood pattern based on the universe's
                # neighbors in the grid, and lattice's neighbors at the edges

                # Degenerate discretization has all universes be different
                if strategy == "degenerate":
                    patterns[(i, j)] = {'locations': [(i, j)]}
                    continue

                # Find neighbors among lattice's neighbors at the edges
                find_edge_neighbors(pattern, i, j)

                # Find neighbors among the lattice's universes
                find_lattice_neighbors(pattern, i, j)

                pattern[1, 1] = key(self.universes[j][i])

                # Look for pattern in dictionary of patterns found
                found = False
                for known_pattern, pattern_data in patterns.items():

                    # Look at all rotations of pattern
                    for rot in range(4):
                        if not found and tuple(map(tuple, pattern)) ==\
                             known_pattern:
                            found = True

                            # Save location of the pattern in the lattice
                            pattern_data['locations'].append((i, j))

                        # Rotate pattern
                        pattern = np.rot90(pattern)

                    # Look at transpose of pattern and its rotations
                    pattern = np.transpose(pattern)
                    for rot in range(4):
                        if not found and tuple(map(tuple, pattern)) ==\
                             known_pattern:
                            found = True

                            # Save location of the pattern in the lattice
                            pattern_data['locations'].append((i, j))

                        # Rotate pattern
                        pattern = np.rot90(pattern)

                    # Transpose pattern back for the next search
                    pattern = np.transpose(pattern)

                # Create new pattern and add to the patterns dictionary
                if not found:
                    patterns[tuple(map(tuple, pattern))] =\
                         {'locations': [(i, j)]}

        # Discretize lattice
        for pattern, pattern_data in patterns.items():

            first_pos = pattern_data['locations'][0]

            # Create a clone of the universe, without cloning materials
            new_universe = self.universes[first_pos[1]][first_pos[0]].clone(
                 clone_materials=False, clone_regions=False)

            # Replace only the materials in materials_to_clone
            for material in materials_to_clone:
                material_cloned = False

                for cell in new_universe.get_all_cells().values():

                    if cell.fill_type == 'material':
                        if cell.fill.id == material.id:

                            # Only a single clone of each material is necessary
                            if not material_cloned:
                                material_clone = material.clone()
                                material_cloned = True

                            cell.fill = material_clone
                    elif cell.fill_type == 'distribmat':
                        raise(ValueError, "Lattice discretization should not "
                              "be used with distributed materials")
                    elif len(cell.temperature) > 1 or len(cell.fill) > 1:
                        raise(ValueError, "Lattice discretization should not "
                              "be used with distributed cells")

            # Rebuild lattice from list of locations with this pattern
            for index, location in enumerate(pattern_data['locations']):
                self.universes[location[1]][location[0]] = new_universe

    def create_xml_subelement(self, xml_element, memo=None):
        """Add the lattice xml representation to an incoming xml element

        Parameters
        ----------
        xml_element : xml.etree.ElementTree.Element
            XML element to be added to

        memo : set or None
            A set of object id's representing geometry entities already
            written to the xml_element. This parameter is used internally
            and should not be specified by users.

        Returns
        -------
        None

        """
        # If the element already contains the Lattice subelement, then return
        if memo and self in memo:
            return
        if memo is not None:
            memo.add(self)

        lattice_subelement = ET.Element("lattice")
        lattice_subelement.set("id", str(self._id))

        if len(self._name) > 0:
            lattice_subelement.set("name", str(self._name))

        # Export the Lattice cell pitch
        pitch = ET.SubElement(lattice_subelement, "pitch")
        pitch.text = ' '.join(map(str, self._pitch))

        # Export the Lattice outer Universe (if specified)
        if self._outer is not None:
            outer = ET.SubElement(lattice_subelement, "outer")
            outer.text = '{0}'.format(self._outer._id)
            self._outer.create_xml_subelement(xml_element, memo)

        # Export Lattice cell dimensions
        dimension = ET.SubElement(lattice_subelement, "dimension")
        dimension.text = ' '.join(map(str, self.shape))

        # Export Lattice lower left
        lower_left = ET.SubElement(lattice_subelement, "lower_left")
        lower_left.text = ' '.join(map(str, self._lower_left))

        # Export the Lattice nested Universe IDs - column major for Fortran
        universe_ids = '\n'

        # 3D Lattices
        if self.ndim == 3:
            for z in range(self.shape[2]):
                for y in range(self.shape[1]):
                    for x in range(self.shape[0]):
                        universe = self._universes[z][y][x]

                        # Append Universe ID to the Lattice XML subelement
                        universe_ids += '{0} '.format(universe._id)

                        # Create XML subelement for this Universe
                        universe.create_xml_subelement(xml_element, memo)

                    # Add newline character when we reach end of row of cells
                    universe_ids += '\n'

                # Add newline character when we reach end of row of cells
                universe_ids += '\n'

        # 2D Lattices
        else:
            for y in range(self.shape[1]):
                for x in range(self.shape[0]):
                    universe = self._universes[y][x]

                    # Append Universe ID to Lattice XML subelement
                    universe_ids += '{0} '.format(universe._id)

                    # Create XML subelement for this Universe
                    universe.create_xml_subelement(xml_element, memo)

                # Add newline character when we reach end of row of cells
                universe_ids += '\n'

        # Remove trailing newline character from Universe IDs string
        universe_ids = universe_ids.rstrip('\n')

        universes = ET.SubElement(lattice_subelement, "universes")
        universes.text = universe_ids

        # Append the XML subelement for this Lattice to the XML element
        xml_element.append(lattice_subelement)

    @classmethod
    def from_xml_element(cls, elem, get_universe):
        """Generate rectangular lattice from XML element

        Parameters
        ----------
        elem : xml.etree.ElementTree.Element
            `<lattice>` element
        get_universe : function
            Function returning universe (defined in
            :meth:`openmc.Geometry.from_xml`)

        Returns
        -------
        RectLattice
            Rectangular lattice

        """
        lat_id = int(get_text(elem, 'id'))
        name = get_text(elem, 'name')
        lat = cls(lat_id, name)
        lat.lower_left = [float(i)
                          for i in get_text(elem, 'lower_left').split()]
        lat.pitch = [float(i) for i in get_text(elem, 'pitch').split()]
        outer = get_text(elem, 'outer')
        if outer is not None:
            lat.outer = get_universe(int(outer))

        # Get array of universes
        dimension = get_text(elem, 'dimension').split()
        shape = np.array(dimension, dtype=int)[::-1]
        uarray = np.array([get_universe(int(i)) for i in
                           get_text(elem, 'universes').split()])
        uarray.shape = shape
        lat.universes = uarray
        return lat

    @classmethod
    def from_hdf5(cls, group, universes):
        """Create rectangular lattice from HDF5 group

        Parameters
        ----------
        group : h5py.Group
            Group in HDF5 file
        universes : dict
            Dictionary mapping universe IDs to instances of
            :class:`openmc.Universe`.

        Returns
        -------
        openmc.RectLattice
            Rectangular lattice

        """
        dimension = group['dimension'][...]
        lower_left = group['lower_left'][...]
        pitch = group['pitch'][...]
        outer = group['outer'][()]
        universe_ids = group['universes'][...]

        # Create the Lattice
        lattice_id = int(group.name.split('/')[-1].lstrip('lattice '))
        name = group['name'][()].decode() if 'name' in group else ''
        lattice = cls(lattice_id, name)
        lattice.lower_left = lower_left
        lattice.pitch = pitch

        # If the Universe specified outer the Lattice is not void
        if outer >= 0:
            lattice.outer = universes[outer]

        # Build array of Universe pointers for the Lattice
        uarray = np.empty(universe_ids.shape, dtype=openmc.Universe)

        for z in range(universe_ids.shape[0]):
            for y in range(universe_ids.shape[1]):
                for x in range(universe_ids.shape[2]):
                    uarray[z, y, x] = universes[universe_ids[z, y, x]]

        # Use 2D NumPy array to store lattice universes for 2D lattices
        if len(dimension) == 2:
            uarray = np.squeeze(uarray)
            uarray = np.atleast_2d(uarray)

        # Set the universes for the lattice
        lattice.universes = uarray

        return lattice


class HexLattice(Lattice):
    r"""A lattice consisting of hexagonal prisms.

    To completely define a hexagonal lattice, the :attr:`HexLattice.center`,
    :attr:`HexLattice.pitch`, :attr:`HexLattice.universes`, and
    :attr:`HexLattice.outer` properties need to be set.

    Most methods for this class use a natural indexing scheme wherein elements
    are assigned an index corresponding to their position relative to skewed
    :math:`(x,\alpha,z)` or :math:`(\alpha,y,z)` bases, depending on the lattice
    orientation, as described fully in :ref:`hexagonal_indexing`. However, note
    that when universes are assigned to lattice elements using the
    :attr:`HexLattice.universes` property, the array indices do not correspond
    to natural indices.

    .. versionchanged:: 0.11
        The orientation of the lattice can now be changed with the
        :attr:`orientation` attribute.

    Parameters
    ----------
    lattice_id : int, optional
        Unique identifier for the lattice. If not specified, an identifier will
        automatically be assigned.
    name : str, optional
        Name of the lattice. If not specified, the name is the empty string.

    Attributes
    ----------
    id : int
        Unique identifier for the lattice
    name : str
        Name of the lattice
    pitch : Iterable of float
        Pitch of the lattice in cm. The first item in the iterable specifies the
        pitch in the radial direction and, if the lattice is 3D, the second item
        in the iterable specifies the pitch in the axial direction.
    outer : openmc.Universe
        A universe to fill all space outside the lattice
    universes : Nested Iterable of openmc.Universe
        A two- or three-dimensional list/array of universes filling each element
        of the lattice. Each sub-list corresponds to one ring of universes and
        should be ordered from outermost ring to innermost ring. The universes
        within each sub-list are ordered from the "top" and proceed in a
        clockwise fashion. The :meth:`HexLattice.show_indices` method can be
        used to help figure out indices for this property.
    center : Iterable of float
        Coordinates of the center of the lattice. If the lattice does not have
        axial sections then only the x- and y-coordinates are specified
    indices : list of tuple
        A list of all possible (z,r,i) or (r,i) lattice element indices that are
        possible, where z is the axial index, r is in the ring index (starting
        from the outermost ring), and i is the index with a ring starting from
        the top and proceeding clockwise.
    orientation : {'x', 'y'}
        str by default 'y' orientation of main lattice diagonal another option
        - 'x'
    num_rings : int
        Number of radial ring positions in the xy-plane
    num_axial : int
        Number of positions along the z-axis.

    """

    def __init__(self, lattice_id=None, name=''):
        super().__init__(lattice_id, name)

        # Initialize Lattice class attributes
        self._num_rings = None
        self._num_axial = None
        self._center = None
        self._orientation = 'y'

    def __repr__(self):
        string = 'HexLattice\n'
        string += '{0: <16}{1}{2}\n'.format('\tID', '=\t', self._id)
        string += '{0: <16}{1}{2}\n'.format('\tName', '=\t', self._name)
        string += '{0: <16}{1}{2}\n'.format('\tOrientation', '=\t',
                                            self._orientation)
        string += '{0: <16}{1}{2}\n'.format('\t# Rings', '=\t', self._num_rings)
        string += '{0: <16}{1}{2}\n'.format('\t# Axial', '=\t', self._num_axial)
        string += '{0: <16}{1}{2}\n'.format('\tCenter', '=\t',
                                            self._center)
        string += '{0: <16}{1}{2}\n'.format('\tPitch', '=\t', self._pitch)

        if self._outer is not None:
            string += '{0: <16}{1}{2}\n'.format('\tOuter', '=\t',
                                                self._outer._id)
        else:
            string += '{0: <16}{1}{2}\n'.format('\tOuter', '=\t',
                                                self._outer)

        string += '{0: <16}\n'.format('\tUniverses')

        if self._num_axial is not None:
            slices = [self._repr_axial_slice(x) for x in self._universes]
            string += '\n'.join(slices)

        else:
            string += self._repr_axial_slice(self._universes)

        return string

    @property
    def num_rings(self):
        return self._num_rings

    @property
    def orientation(self):
        return self._orientation

    @property
    def num_axial(self):
        return self._num_axial

    @property
    def center(self):
        return self._center

    @property
    def indices(self):
        if self.num_axial is None:
            return [(r, i) for r in range(self.num_rings)
                    for i in range(max(6*(self.num_rings - 1 - r), 1))]
        else:
            return [(z, r, i) for z in range(self.num_axial)
                    for r in range(self.num_rings)
                    for i in range(max(6*(self.num_rings - 1 - r), 1))]

    @property
    def _natural_indices(self):
        """Iterate over all possible (x,alpha) or (x,alpha,z) lattice element
        indices.

        This property is used when constructing distributed cell and material
        paths. Most importantly, the iteration order matches that used on the
        Fortran side.

        """
        r = self.num_rings
        if self.num_axial is None:
            for a in range(-r + 1, r):
                for x in range(-r + 1, r):
                    idx = (x, a)
                    if self.is_valid_index(idx):
                        yield idx
        else:
            for z in range(self.num_axial):
                for a in range(-r + 1, r):
                    for x in range(-r + 1, r):
                        idx = (x, a, z)
                        if self.is_valid_index(idx):
                            yield idx

    @property
    def ndim(self):
        return 2 if isinstance(self.universes[0][0], openmc.Universe) else 3

    @center.setter
    def center(self, center):
        cv.check_type('lattice center', center, Iterable, Real)
        cv.check_length('lattice center', center, 2, 3)
        self._center = center

    @orientation.setter
    def orientation(self, orientation):
        cv.check_value('orientation', orientation.lower(), ('x', 'y'))
        self._orientation = orientation.lower()

    @Lattice.pitch.setter
    def pitch(self, pitch):
        cv.check_type('lattice pitch', pitch, Iterable, Real)
        cv.check_length('lattice pitch', pitch, 1, 2)
        for dim in pitch:
            cv.check_greater_than('lattice pitch', dim, 0)
        self._pitch = pitch

    @Lattice.universes.setter
    def universes(self, universes):
        cv.check_iterable_type('lattice universes', universes, openmc.Universe,
                               min_depth=2, max_depth=3)
        self._universes = universes

        # NOTE: This routine assumes that the user creates a "ragged" list of
        # lists, where each sub-list corresponds to one ring of Universes.
        # The sub-lists are ordered from outermost ring to innermost ring.
        # The Universes within each sub-list are ordered from the "top" in a
        # clockwise fashion.

        # Set the number of axial positions.
        if self.ndim == 3:
            self._num_axial = len(self._universes)
        else:
            self._num_axial = None

        # Set the number of rings and make sure this number is consistent for
        # all axial positions.
        if self.ndim == 3:
            self._num_rings = len(self._universes[0])
            for rings in self._universes:
                if len(rings) != self._num_rings:
                    msg = 'HexLattice ID={0:d} has an inconsistent number of ' \
                          'rings per axial positon'.format(self._id)
                    raise ValueError(msg)

        else:
            self._num_rings = len(self._universes)

        # Make sure there are the correct number of elements in each ring.
        if self.ndim == 3:
            for axial_slice in self._universes:
                # Check the center ring.
                if len(axial_slice[-1]) != 1:
                    msg = 'HexLattice ID={0:d} has the wrong number of ' \
                          'elements in the innermost ring. Only 1 element is ' \
                          'allowed in the innermost ring.'.format(self._id)
                    raise ValueError(msg)

                # Check the outer rings.
                for r in range(self._num_rings-1):
                    if len(axial_slice[r]) != 6*(self._num_rings - 1 - r):
                        msg = 'HexLattice ID={0:d} has the wrong number of ' \
                              'elements in ring number {1:d} (counting from the '\
                              'outermost ring). This ring should have {2:d} ' \
                              'elements.'.format(self._id, r,
                                                 6*(self._num_rings - 1 - r))
                        raise ValueError(msg)

        else:
            axial_slice = self._universes
            # Check the center ring.
            if len(axial_slice[-1]) != 1:
                msg = 'HexLattice ID={0:d} has the wrong number of ' \
                      'elements in the innermost ring.  Only 1 element is ' \
                      'allowed in the innermost ring.'.format(self._id)
                raise ValueError(msg)

            # Check the outer rings.
            for r in range(self._num_rings-1):
                if len(axial_slice[r]) != 6*(self._num_rings - 1 - r):
                    msg = 'HexLattice ID={0:d} has the wrong number of ' \
                          'elements in ring number {1:d} (counting from the '\
                          'outermost ring).  This ring should have {2:d} ' \
                          'elements.'.format(self._id, r,
                                             6*(self._num_rings - 1 - r))
                    raise ValueError(msg)

    def find_element(self, point):
        r"""Determine index of lattice element and local coordinates for a point

        Parameters
        ----------
        point : Iterable of float
            Cartesian coordinates of point

        Returns
        -------
        3-tuple of int
            Indices of corresponding lattice element in :math:`(x,\alpha,z)`
            or :math:`(\alpha,y,z)` bases
        numpy.ndarray
            Carestian coordinates of the point in the corresponding lattice
            element coordinate system

        """
        # Convert coordinates to skewed bases
        x = point[0] - self.center[0]
        y = point[1] - self.center[1]
        if self._num_axial is None:
            iz = 1
        else:
            z = point[2] - self.center[2]
            iz = floor(z/self.pitch[1] + 0.5*self.num_axial)
        if self._orientation == 'x':
            alpha = y - x*sqrt(3.)
            i1 = floor(-alpha/(sqrt(3.0) * self.pitch[0]))
            i2 = floor(y/(sqrt(0.75) * self.pitch[0]))
        else:
            alpha = y - x/sqrt(3.)
            i1 = floor(x/(sqrt(0.75) * self.pitch[0]))
            i2 = floor(alpha/self.pitch[0])
        # Check four lattice elements to see which one is closest based on local
        # coordinates
        indices = [(i1, i2, iz), (i1 + 1, i2, iz), (i1, i2 + 1, iz),
                   (i1 + 1, i2 + 1, iz)]
        d_min = np.inf

        for idx in indices:
            p = self.get_local_coordinates(point, idx)
            d = p[0]**2 + p[1]**2
            if d < d_min:
                d_min = d
                idx_min = idx
                p_min = p

        return idx_min, p_min

    def get_local_coordinates(self, point, idx):
        r"""Determine local coordinates of a point within a lattice element

        Parameters
        ----------
        point : Iterable of float
            Cartesian coordinates of point
        idx : Iterable of int
            Indices of lattice element in :math:`(x,\alpha,z)`
            or :math:`(\alpha,y,z)` bases

        Returns
        -------
        3-tuple of float
            Cartesian coordinates of point in the lattice element coordinate
            system

        """
        if self._orientation == 'x':
            x = point[0] - (self.center[0] + (idx[0] + 0.5*idx[1])*self.pitch[0])
            y = point[1] - (self.center[1] + sqrt(0.75)*self.pitch[0]*idx[1])
        else:
            x = point[0] - (self.center[0] + sqrt(0.75)*self.pitch[0]*idx[0])
            y = point[1] - (self.center[1] + (0.5*idx[0] + idx[1])*self.pitch[0])

        if self._num_axial is None:
            z = point[2]
        else:
            z = point[2] - (self.center[2] + (idx[2] + 0.5 - 0.5*self.num_axial) *
                            self.pitch[1])
        return (x, y, z)

    def get_universe_index(self, idx):
        r"""Return index in the universes array corresponding
        to a lattice element index

        Parameters
        ----------
        idx : Iterable of int
            Lattice element indices in the :math:`(x,\alpha,z)` coordinate
            system in 'y' orientation case, or indices in the
            :math:`(\alpha,y,z)` coordinate system in 'x' one

        Returns
        -------
            2- or 3-tuple of int
            Indices used when setting the :attr:`HexLattice.universes`
            property

        """

        # First we determine which ring the index corresponds to.
        x = idx[0]
        a = idx[1]
        z = -a - x
        g = max(abs(x), abs(a), abs(z))

        # Next we use a clever method to figure out where along the ring we are.
        i_ring = self._num_rings - 1 - g
        if x >= 0:
            if a >= 0:
                i_within = x
            else:
                i_within = 2*g + z
        else:
            if a <= 0:
                i_within = 3*g - x
            else:
                i_within = 5*g - z

        if self._orientation == 'x' and g > 0:
            i_within = (i_within + 5*g) % (6*g)

        if self.num_axial is None:
            return (i_ring, i_within)
        else:
            return (idx[2], i_ring, i_within)

    def is_valid_index(self, idx):
        r"""Determine whether lattice element index is within defined range

        Parameters
        ----------
        idx : Iterable of int
            Lattice element indices in the both :math:`(x,\alpha,z)`
            and :math:`(\alpha,y,z)` coordinate system

        Returns
        -------
        bool
            Whether index is valid

        """
        x = idx[0]
        y = idx[1]
        z = 0 - y - x
        g = max(abs(x), abs(y), abs(z))
        if self.num_axial is None:
            return g < self.num_rings
        else:
            return g < self.num_rings and 0 <= idx[2] < self.num_axial

    def create_xml_subelement(self, xml_element, memo=None):
        # If this subelement has already been written, return
        if memo and self in memo:
            return
        if memo is not None:
            memo.add(self)

        lattice_subelement = ET.Element("hex_lattice")
        lattice_subelement.set("id", str(self._id))

        if len(self._name) > 0:
            lattice_subelement.set("name", str(self._name))

        # Export the Lattice cell pitch
        pitch = ET.SubElement(lattice_subelement, "pitch")
        pitch.text = ' '.join(map(str, self._pitch))

        # Export the Lattice outer Universe (if specified)
        if self._outer is not None:
            outer = ET.SubElement(lattice_subelement, "outer")
            outer.text = '{0}'.format(self._outer._id)
            self._outer.create_xml_subelement(xml_element, memo)

        lattice_subelement.set("n_rings", str(self._num_rings))
        # If orientation is "x" export it to XML
        if self._orientation == 'x':
            lattice_subelement.set("orientation", "x")

        if self._num_axial is not None:
            lattice_subelement.set("n_axial", str(self._num_axial))

        # Export Lattice cell center
        center = ET.SubElement(lattice_subelement, "center")
        center.text = ' '.join(map(str, self._center))

        # Export the Lattice nested Universe IDs.

        # 3D Lattices
        if self._num_axial is not None:
            slices = []
            for z in range(self._num_axial):
                # Initialize the center universe.
                universe = self._universes[z][-1][0]
                universe.create_xml_subelement(xml_element, memo)

                # Initialize the remaining universes.
                for r in range(self._num_rings-1):
                    for theta in range(6*(self._num_rings - 1 - r)):
                        universe = self._universes[z][r][theta]
                        universe.create_xml_subelement(xml_element, memo)

                # Get a string representation of the universe IDs.
                slices.append(self._repr_axial_slice(self._universes[z]))

            # Collapse the list of axial slices into a single string.
            universe_ids = '\n'.join(slices)

        # 2D Lattices
        else:
            # Initialize the center universe.
            universe = self._universes[-1][0]
            universe.create_xml_subelement(xml_element, memo)

            # Initialize the remaining universes.
            for r in range(self._num_rings - 1):
                for theta in range(6*(self._num_rings - 1 - r)):
                    universe = self._universes[r][theta]
                    universe.create_xml_subelement(xml_element, memo)

            # Get a string representation of the universe IDs.
            universe_ids = self._repr_axial_slice(self._universes)

        universes = ET.SubElement(lattice_subelement, "universes")
        universes.text = '\n' + universe_ids

        # Append the XML subelement for this Lattice to the XML element
        xml_element.append(lattice_subelement)

    @classmethod
    def from_xml_element(cls, elem, get_universe):
        """Generate hexagonal lattice from XML element

        Parameters
        ----------
        elem : xml.etree.ElementTree.Element
            `<hex_lattice>` element
        get_universe : function
            Function returning universe (defined in
            :meth:`openmc.Geometry.from_xml`)

        Returns
        -------
        HexLattice
            Hexagonal lattice

        """
        lat_id = int(get_text(elem, 'id'))
        name = get_text(elem, 'name')
        lat = cls(lat_id, name)
        lat.center = [float(i) for i in get_text(elem, 'center').split()]
        lat.pitch = [float(i) for i in get_text(elem, 'pitch').split()]
        lat.orientation = get_text(elem, 'orientation', 'y')
        outer = get_text(elem, 'outer')
        if outer is not None:
            lat.outer = get_universe(int(outer))

        # Get nested lists of universes
        lat._num_rings = n_rings = int(get_text(elem, 'n_rings'))
        lat._num_axial = n_axial = int(get_text(elem, 'n_axial', 1))

        # Create empty nested lists for one axial level
        univs = [[None for _ in range(max(6*(n_rings - 1 - r), 1))]
                 for r in range(n_rings)]
        if n_axial > 1:
            univs = [deepcopy(univs) for i in range(n_axial)]

        # Get flat array of universes
        uarray = np.array([get_universe(int(i)) for i in
                           get_text(elem, 'universes').split()])

        # Fill nested lists
        j = 0
        for z in range(n_axial):
            # Get list for a single axial level
            axial_level = univs[z] if n_axial > 1 else univs

            if lat.orientation == 'y':
                # Start iterating from top
                x, alpha = 0, n_rings - 1
                while True:
                    # Set entry in list based on (x,alpha,z) coordinates
                    _, i_ring, i_within = lat.get_universe_index((x, alpha, z))
                    axial_level[i_ring][i_within] = uarray[j]

                    # Move to the right
                    x += 2
                    alpha -= 1
                    if not lat.is_valid_index((x, alpha, z)):
                        # Move down in y direction
                        alpha += x - 1
                        x = 1 - x
                        if not lat.is_valid_index((x, alpha, z)):
                            # Move to the right
                            x += 2
                            alpha -= 1
                            if not lat.is_valid_index((x, alpha, z)):
                                # Reached the bottom
                                break
                    j += 1
            else:
                # Start iterating from top
                alpha, y = 1 - n_rings, n_rings - 1
                while True:
                    # Set entry in list based on (alpha,y,z) coordinates
                    _, i_ring, i_within = lat.get_universe_index((alpha, y, z))
                    axial_level[i_ring][i_within] = uarray[j]

                    # Move to the right
                    alpha += 1
                    if not lat.is_valid_index((alpha, y, z)):
                        # Move down to next row
                        alpha = 1 - n_rings
                        y -= 1

                        # Check if we've reached the bottom
                        if y == -n_rings:
                            break

                        while not lat.is_valid_index((alpha, y, z)):
                            # Move to the right
                            alpha += 1
                    j += 1

        lat.universes = univs
        return lat

    def _repr_axial_slice(self, universes):
        """Return string representation for the given 2D group of universes.

        The 'universes' argument should be a list of lists of universes where
        each sub-list represents a single ring.  The first list should be the
        outer ring.
        """
        if self._orientation == 'x':
            return self._repr_axial_slice_x(universes)
        else:
            return self._repr_axial_slice_y(universes)

    def _repr_axial_slice_x(self, universes):
        """Return string representation for the given 2D group of universes
        in 'x' orientation case.

        The 'universes' argument should be a list of lists of universes where
        each sub-list represents a single ring.  The first list should be the
        outer ring.
        """

        # Find the largest universe ID and count the number of digits so we can
        # properly pad the output string later.
        largest_id = max([max([univ._id for univ in ring])
                          for ring in universes])
        n_digits = len(str(largest_id))
        pad = ' '*n_digits
        id_form = '{: ^' + str(n_digits) + 'd}'

        # Initialize the list for each row.
        rows = [[] for i in range(2*self._num_rings - 1)]
        middle = self._num_rings - 1

        # Start with the degenerate first ring.
        universe = universes[-1][0]
        rows[middle] = [id_form.format(universe._id)]

        # Add universes one ring at a time.
        for r in range(1, self._num_rings):
            # r_prime increments down while r increments up.
            r_prime = self._num_rings - 1 - r
            theta = 0
            y = middle

            # Climb down the bottom-right
            for i in range(r):
                # Add the universe.
                universe = universes[r_prime][theta]
                rows[y].append(id_form.format(universe._id))

                # Translate the indices.
                y += 1
                theta += 1

            # Climb left across the bottom
            for i in range(r):
                # Add the universe.
                universe = universes[r_prime][theta]
                rows[y].insert(0, id_form.format(universe._id))

                # Translate the indices.
                theta += 1

            # Climb up the bottom-left
            for i in range(r):
                # Add the universe.
                universe = universes[r_prime][theta]
                rows[y].insert(0, id_form.format(universe._id))

                # Translate the indices.
                y -= 1
                theta += 1

            # Climb up the top-left
            for i in range(r):
                # Add the universe.
                universe = universes[r_prime][theta]
                rows[y].insert(0, id_form.format(universe._id))

                # Translate the indices.
                y -= 1
                theta += 1

            # Climb right across the top
            for i in range(r):
                # Add the universe.
                universe = universes[r_prime][theta]
                rows[y].append(id_form.format(universe._id))

                # Translate the indices.
                theta += 1

            # Climb down the top-right
            for i in range(r):
                # Add the universe.
                universe = universes[r_prime][theta]
                rows[y].append(id_form.format(universe._id))

                # Translate the indices.
                y += 1
                theta += 1

        # Flip the rows and join each row into a single string.
        rows = [pad.join(x) for x in rows]

        # Pad the beginning of the rows so they line up properly.
        for y in range(self._num_rings - 1):
            rows[y] = (self._num_rings - 1 - y)*pad + rows[y]
            rows[-1 - y] = (self._num_rings - 1 - y)*pad + rows[-1 - y]

        # Join the rows together and return the string.
        universe_ids = '\n'.join(rows)
        return universe_ids

    def _repr_axial_slice_y(self, universes):
        """Return string representation for the given 2D group of universes in
        'y' orientation case..

        The 'universes' argument should be a list of lists of universes where
        each sub-list represents a single ring.  The first list should be the
        outer ring.
        """

        # Find the largest universe ID and count the number of digits so we can
        # properly pad the output string later.
        largest_id = max([max([univ._id for univ in ring])
                          for ring in universes])
        n_digits = len(str(largest_id))
        pad = ' '*n_digits
        id_form = '{: ^' + str(n_digits) + 'd}'

        # Initialize the list for each row.
        rows = [[] for i in range(1 + 4 * (self._num_rings-1))]
        middle = 2 * (self._num_rings - 1)

        # Start with the degenerate first ring.
        universe = universes[-1][0]
        rows[middle] = [id_form.format(universe._id)]

        # Add universes one ring at a time.
        for r in range(1, self._num_rings):
            # r_prime increments down while r increments up.
            r_prime = self._num_rings - 1 - r
            theta = 0
            y = middle + 2*r

            # Climb down the top-right.
            for i in range(r):
                # Add the universe.
                universe = universes[r_prime][theta]
                rows[y].append(id_form.format(universe._id))

                # Translate the indices.
                y -= 1
                theta += 1

            # Climb down the right.
            for i in range(r):
                # Add the universe.
                universe = universes[r_prime][theta]
                rows[y].append(id_form.format(universe._id))

                # Translate the indices.
                y -= 2
                theta += 1

            # Climb down the bottom-right.
            for i in range(r):
                # Add the universe.
                universe = universes[r_prime][theta]
                rows[y].append(id_form.format(universe._id))

                # Translate the indices.
                y -= 1
                theta += 1

            # Climb up the bottom-left.
            for i in range(r):
                # Add the universe.
                universe = universes[r_prime][theta]
                rows[y].insert(0, id_form.format(universe._id))

                # Translate the indices.
                y += 1
                theta += 1

            # Climb up the left.
            for i in range(r):
                # Add the universe.
                universe = universes[r_prime][theta]
                rows[y].insert(0, id_form.format(universe._id))

                # Translate the indices.
                y += 2
                theta += 1

            # Climb up the top-left.
            for i in range(r):
                # Add the universe.
                universe = universes[r_prime][theta]
                rows[y].insert(0, id_form.format(universe._id))

                # Translate the indices.
                y += 1
                theta += 1

        # Flip the rows and join each row into a single string.
        rows = [pad.join(x) for x in rows[::-1]]

        # Pad the beginning of the rows so they line up properly.
        for y in range(self._num_rings - 1):
            rows[y] = (self._num_rings - 1 - y)*pad + rows[y]
            rows[-1 - y] = (self._num_rings - 1 - y)*pad + rows[-1 - y]

        for y in range(self._num_rings % 2, self._num_rings, 2):
            rows[middle + y] = pad + rows[middle + y]
            if y != 0:
                rows[middle - y] = pad + rows[middle - y]

        # Join the rows together and return the string.
        universe_ids = '\n'.join(rows)
        return universe_ids

    @staticmethod
    def _show_indices_y(num_rings):
        """Return a diagram of the hexagonal lattice layout with indices.

        This method can be used to show the proper indices to be used when
        setting the :attr:`HexLattice.universes` property. For example, running
        this method with num_rings=3 will return the following diagram::

                      (0, 0)
                (0,11)      (0, 1)
          (0,10)      (1, 0)      (0, 2)
                (1, 5)      (1, 1)
          (0, 9)      (2, 0)      (0, 3)
                (1, 4)      (1, 2)
          (0, 8)      (1, 3)      (0, 4)
                (0, 7)      (0, 5)
                      (0, 6)

        Parameters
        ----------
        num_rings : int
            Number of rings in the hexagonal lattice

        Returns
        -------
        str
            Diagram of the hexagonal lattice showing indices

        """

        # Find the largest string and count the number of digits so we can
        # properly pad the output string later
        largest_index = 6*(num_rings - 1)
        n_digits_index = len(str(largest_index))
        n_digits_ring = len(str(num_rings - 1))
        str_form = '({{:{}}},{{:{}}})'.format(n_digits_ring, n_digits_index)
        pad = ' '*(n_digits_index + n_digits_ring + 3)

        # Initialize the list for each row.
        rows = [[] for i in range(1 + 4 * (num_rings-1))]
        middle = 2 * (num_rings - 1)

        # Start with the degenerate first ring.
        rows[middle] = [str_form.format(num_rings - 1, 0)]

        # Add universes one ring at a time.
        for r in range(1, num_rings):
            # r_prime increments down while r increments up.
            r_prime = num_rings - 1 - r
            theta = 0
            y = middle + 2*r

            for i in range(r):
                # Climb down the top-right.
                rows[y].append(str_form.format(r_prime, theta))
                y -= 1
                theta += 1

            for i in range(r):
                # Climb down the right.
                rows[y].append(str_form.format(r_prime, theta))
                y -= 2
                theta += 1

            for i in range(r):
                # Climb down the bottom-right.
                rows[y].append(str_form.format(r_prime, theta))
                y -= 1
                theta += 1

            for i in range(r):
                # Climb up the bottom-left.
                rows[y].insert(0, str_form.format(r_prime, theta))
                y += 1
                theta += 1

            for i in range(r):
                # Climb up the left.
                rows[y].insert(0, str_form.format(r_prime, theta))
                y += 2
                theta += 1

            for i in range(r):
                # Climb up the top-left.
                rows[y].insert(0, str_form.format(r_prime, theta))
                y += 1
                theta += 1

        # Flip the rows and join each row into a single string.
        rows = [pad.join(x) for x in rows[::-1]]

        # Pad the beginning of the rows so they line up properly.
        for y in range(num_rings - 1):
            rows[y] = (num_rings - 1 - y)*pad + rows[y]
            rows[-1 - y] = (num_rings - 1 - y)*pad + rows[-1 - y]

        for y in range(num_rings % 2, num_rings, 2):
            rows[middle + y] = pad + rows[middle + y]
            if y != 0:
                rows[middle - y] = pad + rows[middle - y]

        # Join the rows together and return the string.
        return '\n'.join(rows)

    @staticmethod
    def _show_indices_x(num_rings):
        """Return a diagram of the hexagonal lattice with x orientation
        layout with indices.

        This method can be used to show the proper indices to be used when
        setting the :attr:`HexLattice.universes` property. For example,running
        this method with num_rings=3 will return the similar diagram::

                      (0, 8)      (0, 9)      (0,10)

                (0, 7)      (1, 4)      (1, 5)      (0,11)

          (0, 6)      (1, 3)      (2, 0)      (1, 0)      (0, 0)

                (0, 5)      (1, 2)      (1, 1)      (0, 1)

                      (0, 4)      (0, 3)      (0, 2)

        Parameters
        ----------
        num_rings : int
            Number of rings in the hexagonal lattice

        Returns
        -------
        str
            Diagram of the hexagonal lattice showing indices in OX orientation

        """

        # Find the largest string and count the number of digits so we can
        # properly pad the output string later
        largest_index = 6*(num_rings - 1)
        n_digits_index = len(str(largest_index))
        n_digits_ring = len(str(num_rings - 1))
        str_form = '({{:{}}},{{:{}}})'.format(n_digits_ring, n_digits_index)
        pad = ' '*(n_digits_index + n_digits_ring + 3)

        # Initialize the list for each row.
        rows = [[] for i in range(2*num_rings - 1)]
        middle = num_rings - 1

        # Start with the degenerate first ring.
        rows[middle] = [str_form.format(num_rings - 1, 0)]

        # Add universes one ring at a time.
        for r in range(1, num_rings):
            # r_prime increments down while r increments up.
            r_prime = num_rings - 1 - r
            theta = 0
            y = middle

            for i in range(r):
                # Climb down the bottom-right
                rows[y].append(str_form.format(r_prime, theta))
                y += 1
                theta += 1

            for i in range(r):
                # Climb left across the bottom
                rows[y].insert(0, str_form.format(r_prime, theta))
                theta += 1

            for i in range(r):
                # Climb up the bottom-left
                rows[y].insert(0, str_form.format(r_prime, theta))
                y -= 1
                theta += 1

            for i in range(r):
                # Climb up the top-left
                rows[y].insert(0, str_form.format(r_prime, theta))
                y -= 1
                theta += 1

            for i in range(r):
                # Climb right across the top
                rows[y].append(str_form.format(r_prime, theta))
                theta += 1

            for i in range(r):
                # Climb down the top-right
                rows[y].append(str_form.format(r_prime, theta))
                y += 1
                theta += 1

        # Flip the rows and join each row into a single string.
        rows = [pad.join(x) for x in rows]

        # Pad the beginning of the rows so they line up properly.
        for y in range(num_rings - 1):
            rows[y] = (num_rings - 1 - y)*pad + rows[y]
            rows[-1 - y] = (num_rings - 1 - y)*pad + rows[-1 - y]

        # Join the rows together and return the string.
        return '\n\n'.join(rows)

    @staticmethod
    def show_indices(num_rings, orientation="y"):
        """Return a diagram of the hexagonal lattice layout with indices.

        Parameters
        ----------
        num_rings : int
            Number of rings in the hexagonal lattice
        orientation : {"x", "y"}
            Orientation of the hexagonal lattice

        Returns
        -------
        str
            Diagram of the hexagonal lattice showing indices

        """

        if orientation == 'x':
            return HexLattice._show_indices_x(num_rings)
        else:
            return HexLattice._show_indices_y(num_rings)

    @classmethod
    def from_hdf5(cls, group, universes):
        """Create rectangular lattice from HDF5 group

        Parameters
        ----------
        group : h5py.Group
            Group in HDF5 file
        universes : dict
            Dictionary mapping universe IDs to instances of
            :class:`openmc.Universe`.

        Returns
        -------
        openmc.RectLattice
            Rectangular lattice

        """
        n_rings = group['n_rings'][()]
        n_axial = group['n_axial'][()]
        center = group['center'][()]
        pitch = group['pitch'][()]
        outer = group['outer'][()]
        if 'orientation' in group:
            orientation = group['orientation'][()].decode()
        else:
            orientation = "y"
        universe_ids = group['universes'][()]

        # Create the Lattice
        lattice_id = int(group.name.split('/')[-1].lstrip('lattice '))
        name = group['name'][()].decode() if 'name' in group else ''
        lattice = openmc.HexLattice(lattice_id, name)
        lattice.center = center
        lattice.pitch = pitch
        lattice.orientation = orientation
        # If the Universe specified outer the Lattice is not void
        if outer >= 0:
            lattice.outer = universes[outer]
        if orientation == "y":
            # Build array of Universe pointers for the Lattice.  Note that
            # we need to convert between the HDF5's square array of
            # (x, alpha, z) to the Python API's format of a ragged nested
            # list of (z, ring, theta).
            uarray = []
            for z in range(n_axial):
                # Add a list for this axial level.
                uarray.append([])
                x = n_rings - 1
                a = 2*n_rings - 2
                for r in range(n_rings - 1, 0, -1):
                    # Add a list for this ring.
                    uarray[-1].append([])

                    # Climb down the top-right.
                    for i in range(r):
                        uarray[-1][-1].append(universe_ids[z, a, x])
                        x += 1
                        a -= 1

                    # Climb down the right.
                    for i in range(r):
                        uarray[-1][-1].append(universe_ids[z, a, x])
                        a -= 1

                    # Climb down the bottom-right.
                    for i in range(r):
                        uarray[-1][-1].append(universe_ids[z, a, x])
                        x -= 1

                    # Climb up the bottom-left.
                    for i in range(r):
                        uarray[-1][-1].append(universe_ids[z, a, x])
                        x -= 1
                        a += 1

                    # Climb up the left.
                    for i in range(r):
                        uarray[-1][-1].append(universe_ids[z, a, x])
                        a += 1

                    # Climb up the top-left.
                    for i in range(r):
                        uarray[-1][-1].append(universe_ids[z, a, x])
                        x += 1

                    # Move down to the next ring.
                    a -= 1

                    # Convert the ids into Universe objects.
                    uarray[-1][-1] = [universes[u_id]
                                        for u_id in uarray[-1][-1]]

                # Handle the degenerate center ring separately.
                u_id = universe_ids[z, a, x]
                uarray[-1].append([universes[u_id]])
        else:
            # Build array of Universe pointers for the Lattice.  Note that
            # we need to convert between the HDF5's square array of
            # (alpha, y, z) to the Python API's format of a ragged nested
            # list of (z, ring, theta).
            uarray = []
            for z in range(n_axial):
                # Add a list for this axial level.
                uarray.append([])
                a = 2*n_rings - 2
                y = n_rings - 1
                for r in range(n_rings - 1, 0, -1):
                    # Add a list for this ring.
                    uarray[-1].append([])

                    # Climb down the bottom-right.
                    for i in range(r):
                        uarray[-1][-1].append(universe_ids[z, y, a])
                        y -= 1

                    # Climb across the bottom.
                    for i in range(r):
                        uarray[-1][-1].append(universe_ids[z, y, a])
                        a -= 1

                    # Climb up the bottom-left.
                    for i in range(r):
                        uarray[-1][-1].append(universe_ids[z, y, a])
                        a -= 1
                        y += 1

                    # Climb up the top-left.
                    for i in range(r):
                        uarray[-1][-1].append(universe_ids[z, y, a])
                        y += 1

                    # Climb across the top.
                    for i in range(r):
                        uarray[-1][-1].append(universe_ids[z, y, a])
                        a += 1

                    # Climb down the top-right.
                    for i in range(r):
                        uarray[-1][-1].append(universe_ids[z, y, a])
                        a += 1
                        y -= 1

                    # Move down to the next ring.
                    a -= 1

                    # Convert the ids into Universe objects.
                    uarray[-1][-1] = [universes[u_id]
                                        for u_id in uarray[-1][-1]]

                # Handle the degenerate center ring separately.
                u_id = universe_ids[z, y, a]
                uarray[-1].append([universes[u_id]])

        # Add the universes to the lattice.
        if len(pitch) == 2:
            # Lattice is 3D
            lattice.universes = uarray
        else:
            # Lattice is 2D; extract the only axial level
            lattice.universes = uarray[0]

        return lattice
