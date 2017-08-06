from collections import OrderedDict, Iterable
from copy import deepcopy
from math import cos, sin, pi
from numbers import Real, Integral
from xml.etree import ElementTree as ET
import sys
import warnings

from six import string_types
import numpy as np

import openmc
import openmc.checkvalue as cv
from openmc.surface import Halfspace
from openmc.region import Region, Intersection, Complement
from .mixin import IDManagerMixin


class Cell(IDManagerMixin):
    r"""A region of space defined as the intersection of half-space created by
    quadric surfaces.

    Parameters
    ----------
    cell_id : int, optional
        Unique identifier for the cell. If not specified, an identifier will
        automatically be assigned.
    name : str, optional
        Name of the cell. If not specified, the name is the empty string.
    fill : openmc.Material or openmc.Universe or openmc.Lattice or None or iterable of openmc.Material, optional
        Indicates what the region of space is filled with
    region : openmc.Region, optional
        Region of space that is assigned to the cell.

    Attributes
    ----------
    id : int
        Unique identifier for the cell
    name : str
        Name of the cell
    fill : openmc.Material or openmc.Universe or openmc.Lattice or None or iterable of openmc.Material
        Indicates what the region of space is filled with. If None, the cell is
        treated as a void. An iterable of materials is used to fill repeated
        instances of a cell with different materials.
    fill_type : {'material', 'universe', 'lattice', 'distribmat', 'void'}
        Indicates what the cell is filled with.
    region : openmc.Region or None
        Region of space that is assigned to the cell.
    rotation : Iterable of float
        If the cell is filled with a universe, this array specifies the angles
        in degrees about the x, y, and z axes that the filled universe should be
        rotated. The rotation applied is an intrinsic rotation with specified
        Tait-Bryan angles. That is to say, if the angles are :math:`(\phi,
        \theta, \psi)`, then the rotation matrix applied is :math:`R_z(\psi)
        R_y(\theta) R_x(\phi)` or

        .. math::

           \left [ \begin{array}{ccc} \cos\theta \cos\psi & -\cos\theta \sin\psi
           + \sin\phi \sin\theta \cos\psi & \sin\phi \sin\psi + \cos\phi
           \sin\theta \cos\psi \\ \cos\theta \sin\psi & \cos\phi \cos\psi +
           \sin\phi \sin\theta \sin\psi & -\sin\phi \cos\psi + \cos\phi
           \sin\theta \sin\psi \\ -\sin\theta & \sin\phi \cos\theta & \cos\phi
           \cos\theta \end{array} \right ]
    rotation_matrix : numpy.ndarray
        The rotation matrix defined by the angles specified in the
        :attr:`Cell.rotation` property.
    temperature : float or iterable of float
        Temperature of the cell in Kelvin.  Multiple temperatures can be given
        to give each distributed cell instance a unique temperature.
    translation : Iterable of float
        If the cell is filled with a universe, this array specifies a vector
        that is used to translate (shift) the universe.
    paths : list of str
        The paths traversed through the CSG tree to reach each cell
        instance. This property is initialized by calling the
        :meth:`Geometry.determine_paths` method.
    num_instances : int
        The number of instances of this cell throughout the geometry.
    volume : float
        Volume of the cell in cm^3. This can either be set manually or
        calculated in a stochastic volume calculation and added via the
        :meth:`Cell.add_volume_information` method.

    """

    next_id = 1
    used_ids = set()

    def __init__(self, cell_id=None, name='', fill=None, region=None):
        # Initialize Cell class attributes
        self.id = cell_id
        self.name = name
        self.fill = fill
        self.region = region
        self._rotation = None
        self._rotation_matrix = None
        self._temperature = None
        self._translation = None
        self._paths = None
        self._num_instances = None
        self._volume = None
        self._atoms = None

    def __contains__(self, point):
        if self.region is None:
            return True
        else:
            return point in self.region

    def __eq__(self, other):
        if not isinstance(other, Cell):
            return False
        elif self.id != other.id:
            return False
        elif self.name != other.name:
            return False
        elif self.fill != other.fill:
            return False
        elif self.region != other.region:
            return False
        elif self.rotation != other.rotation:
            return False
        elif self.temperature != other.temperature:
            return False
        elif self.translation != other.translation:
            return False
        else:
            return True

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(repr(self))

    def __repr__(self):
        string = 'Cell\n'
        string += '{: <16}=\t{}\n'.format('\tID', self.id)
        string += '{: <16}=\t{}\n'.format('\tName', self.name)

        if self.fill_type == 'material':
            string += '{: <16}=\tMaterial {}\n'.format('\tFill', self.fill.id)
        elif self.fill_type == 'void':
            string += '{: <16}=\tNone\n'.format('\tFill')
        elif self.fill_type == 'distribmat':
            string += '{: <16}=\t{}\n'.format('\tFill', list(map(
                lambda m: m if m is None else m.id, self.fill)))
        else:
            string += '{: <16}=\t{}\n'.format('\tFill', self.fill.id)

        string += '{: <16}=\t{}\n'.format('\tRegion', self.region)
        string += '{: <16}=\t{}\n'.format('\tRotation', self.rotation)
        if self.fill_type == 'material':
            string += '\t{0: <15}=\t{1}\n'.format('Temperature',
                                                  self.temperature)
        string += '{: <16}=\t{}\n'.format('\tTranslation', self.translation)

        return string

    @property
    def name(self):
        return self._name

    @property
    def fill(self):
        return self._fill

    @property
    def fill_type(self):
        if isinstance(self.fill, openmc.Material):
            return 'material'
        elif isinstance(self.fill, openmc.Universe):
            return 'universe'
        elif isinstance(self.fill, openmc.Lattice):
            return 'lattice'
        elif isinstance(self.fill, Iterable):
            return 'distribmat'
        else:
            return 'void'

    @property
    def region(self):
        return self._region

    @property
    def rotation(self):
        return self._rotation

    @property
    def rotation_matrix(self):
        return self._rotation_matrix

    @property
    def temperature(self):
        return self._temperature

    @property
    def translation(self):
        return self._translation

    @property
    def volume(self):
        return self._volume

    @property
    def paths(self):
        if self._paths is None:
            raise ValueError('Cell instance paths have not been determined. '
                             'Call the Geometry.determine_paths() method.')
        return self._paths

    @property
    def bounding_box(self):
        if self.region is not None:
            return self.region.bounding_box
        else:
            return (np.array([-np.inf, -np.inf, -np.inf]),
                    np.array([np.inf, np.inf, np.inf]))

    @property
    def num_instances(self):
        if self._num_instances is None:
            raise ValueError(
                'Number of cell instances have not been determined. Call the '
                'Geometry.determine_paths() method.')
        return self._num_instances

    @name.setter
    def name(self, name):
        if name is not None:
            cv.check_type('cell name', name, string_types)
            self._name = name
        else:
            self._name = ''

    @fill.setter
    def fill(self, fill):
        if fill is not None:
            if isinstance(fill, string_types):
                if fill.strip().lower() != 'void':
                    msg = 'Unable to set Cell ID="{0}" to use a non-Material ' \
                          'or Universe fill "{1}"'.format(self._id, fill)
                    raise ValueError(msg)
                fill = None

            elif isinstance(fill, Iterable):
                for i, f in enumerate(fill):
                    if f is not None:
                        cv.check_type('cell.fill[i]', f, openmc.Material)

            elif not isinstance(fill, (openmc.Material, openmc.Lattice,
                                       openmc.Universe)):
                msg = 'Unable to set Cell ID="{0}" to use a non-Material or ' \
                      'Universe fill "{1}"'.format(self._id, fill)
                raise ValueError(msg)

        self._fill = fill

    @rotation.setter
    def rotation(self, rotation):
        if not isinstance(self.fill, openmc.Universe):
            raise TypeError('Cell rotation can only be applied if the cell '
                            'is filled with a Universe.')

        cv.check_type('cell rotation', rotation, Iterable, Real)
        cv.check_length('cell rotation', rotation, 3)
        self._rotation = np.asarray(rotation)

        # Save rotation matrix -- the reason we do this instead of having it be
        # automatically calculated when the rotation_matrix property is accessed
        # is so that plotting on a rotated geometry can be done faster.
        phi, theta, psi = self.rotation*(-pi/180.)
        c3, s3 = cos(phi), sin(phi)
        c2, s2 = cos(theta), sin(theta)
        c1, s1 = cos(psi), sin(psi)
        self._rotation_matrix = np.array([
            [c1*c2, c1*s2*s3 - c3*s1, s1*s3 + c1*c3*s2],
            [c2*s1, c1*c3 + s1*s2*s3, c3*s1*s2 - c1*s3],
            [-s2, c2*s3, c2*c3]])

    @translation.setter
    def translation(self, translation):
        cv.check_type('cell translation', translation, Iterable, Real)
        cv.check_length('cell translation', translation, 3)
        self._translation = np.asarray(translation)

    @temperature.setter
    def temperature(self, temperature):
        # Make sure temperatures are positive
        cv.check_type('cell temperature', temperature, (Iterable, Real))
        if isinstance(temperature, Iterable):
            cv.check_type('cell temperature', temperature, Iterable, Real)
            for T in temperature:
                cv.check_greater_than('cell temperature', T, 0.0, True)
        else:
            cv.check_greater_than('cell temperature', temperature, 0.0, True)

        # If this cell is filled with a universe or lattice, propagate
        # temperatures to all cells contained. Otherwise, simply assign it.
        if self.fill_type in ('universe', 'lattice'):
            for c in self.get_all_cells().values():
                if c.fill_type == 'material':
                    c._temperature = temperature
        else:
            self._temperature = temperature

    @region.setter
    def region(self, region):
        if region is not None:
            cv.check_type('cell region', region, Region)
        self._region = region

    @volume.setter
    def volume(self, volume):
        if volume is not None:
            cv.check_type('cell volume', volume, Real)
        self._volume = volume

    def add_surface(self, surface, halfspace):
        """Add a half-space to the list of half-spaces whose intersection defines the
        cell.

        .. deprecated:: 0.7.1
            Use the :attr:`Cell.region` property to directly specify a Region
            expression.

        Parameters
        ----------
        surface : openmc.Surface
            Quadric surface dividing space
        halfspace : {-1, 1}
            Indicate whether the negative or positive half-space is to be used

        """

        warnings.warn("Cell.add_surface(...) has been deprecated and may be "
                      "removed in a future version. The region for a Cell "
                      "should be defined using the region property directly.",
                      DeprecationWarning)

        if not isinstance(surface, openmc.Surface):
            msg = 'Unable to add Surface "{0}" to Cell ID="{1}" since it is ' \
                        'not a Surface object'.format(surface, self._id)
            raise ValueError(msg)

        if halfspace not in [-1, +1]:
            msg = 'Unable to add Surface "{0}" to Cell ID="{1}" with halfspace ' \
                  '"{2}" since it is not +/-1'.format(surface, self._id, halfspace)
            raise ValueError(msg)

        # If no region has been assigned, simply use the half-space. Otherwise,
        # take the intersection of the current region and the half-space
        # specified
        region = +surface if halfspace == 1 else -surface
        if self.region is None:
            self.region = region
        else:
            if isinstance(self.region, Intersection):
                self.region &= region
            else:
                self.region = Intersection(self.region, region)

    def add_volume_information(self, volume_calc):
        """Add volume information to a cell.

        Parameters
        ----------
        volume_calc : openmc.VolumeCalculation
            Results from a stochastic volume calculation

        """
        if volume_calc.domain_type == 'cell':
            if self.id in volume_calc.volumes:
                self._volume = volume_calc.volumes[self.id][0]
                self._atoms = volume_calc.atoms[self.id]
            else:
                raise ValueError('No volume information found for this cell.')
        else:
            raise ValueError('No volume information found for this cell.')

    def get_nuclides(self):
        """Returns all nuclides in the cell

        Returns
        -------
        nuclides : list of str
            List of nuclide names

        """
        return self.fill.get_nuclides() if self.fill_type != 'void' else []

    def get_nuclide_densities(self):
        """Return all nuclides contained in the cell and their densities

        Returns
        -------
        nuclides : collections.OrderedDict
            Dictionary whose keys are nuclide names and values are 2-tuples of
            (nuclide, density)

        """

        nuclides = OrderedDict()

        if self.fill_type == 'material':
            nuclides.update(self.fill.get_nuclide_densities())
        elif self.fill_type == 'void':
            pass
        else:
            if self._atoms is not None:
                volume = self.volume
                for name, atoms in self._atoms.items():
                    nuclide = openmc.Nuclide(name)
                    density = 1.0e-24 * atoms[0]/volume  # density in atoms/b-cm
                    nuclides[name] = (nuclide, density)
            else:
                raise RuntimeError(
                    'Volume information is needed to calculate microscopic cross '
                    'sections for cell {}. This can be done by running a '
                    'stochastic volume calculation via the '
                    'openmc.VolumeCalculation object'.format(self.id))

        return nuclides

    def get_all_cells(self):
        """Return all cells that are contained within this one if it is filled with a
        universe or lattice

        Returns
        -------
        cells : collections.orderedDict
            Dictionary whose keys are cell IDs and values are :class:`Cell`
            instances

        """

        cells = OrderedDict()

        if self.fill_type in ('universe', 'lattice'):
            cells.update(self.fill.get_all_cells())

        return cells

    def get_all_materials(self):
        """Return all materials that are contained within the cell

        Returns
        -------
        materials : collections.OrderedDict
            Dictionary whose keys are material IDs and values are
            :class:`Material` instances

        """
        materials = OrderedDict()
        if self.fill_type == 'material':
            materials[self.fill.id] = self.fill
        elif self.fill_type == 'distribmat':
            for m in self.fill:
                if m is not None:
                    materials[m.id] = m
        else:
            # Append all Cells in each Cell in the Universe to the dictionary
            cells = self.get_all_cells()
            for cell in cells.values():
                materials.update(cell.get_all_materials())

        return materials

    def get_all_universes(self):
        """Return all universes that are contained within this one if any of
        its cells are filled with a universe or lattice.

        Returns
        -------
        universes : collections.OrderedDict
            Dictionary whose keys are universe IDs and values are
            :class:`Universe` instances

        """

        universes = OrderedDict()

        if self.fill_type == 'universe':
            universes[self.fill.id] = self.fill
            universes.update(self.fill.get_all_universes())
        elif self.fill_type == 'lattice':
            universes.update(self.fill.get_all_universes())

        return universes

    def clone(self, memo=None):
        """Create a copy of this cell with a new unique ID, and clones
        the cell's region and fill.

        Parameters
        ----------
        memo : dict or None
            A nested dictionary of previously cloned objects. This parameter
            is used internally and should not be specified by the user.

        Returns
        -------
        clone : openmc.Cell
            The clone of this cell

        """

        if memo is None:
            memo = {}

        # If no nemoize'd clone exists, instantiate one
        if self not in memo:
            # Temporarily remove paths
            paths = self._paths
            self._paths = None

            clone = deepcopy(self)
            clone.id = None
            clone._num_instances = None

            # Restore paths on original instance
            self._paths = paths

            if self.region is not None:
                clone.region = self.region.clone(memo)
            if self.fill is not None:
                if self.fill_type == 'distribmat':
                    clone.fill = [fill.clone(memo) if fill is not None else None
                                  for fill in self.fill]
                else:
                    clone.fill = self.fill.clone(memo)

            # Memoize the clone
            memo[self] = clone

        return memo[self]

    def create_xml_subelement(self, xml_element):
        element = ET.Element("cell")
        element.set("id", str(self.id))

        if len(self._name) > 0:
            element.set("name", str(self.name))

        if self.fill_type == 'void':
            element.set("material", "void")

        elif self.fill_type == 'material':
            element.set("material", str(self.fill.id))

        elif self.fill_type == 'distribmat':
            element.set("material", ' '.join(['void' if m is None else str(m.id)
                                              for m in self.fill]))

        elif self.fill_type in ('universe', 'lattice'):
            element.set("fill", str(self.fill.id))
            self.fill.create_xml_subelement(xml_element)

        if self.region is not None:
            # Set the region attribute with the region specification
            region = str(self.region)
            if region.startswith('('):
                region = region[1:-1]
            if len(region) > 0:
                element.set("region", region)

            # Only surfaces that appear in a region are added to the geometry
            # file, so the appropriate check is performed here. First we create
            # a function which is called recursively to navigate through the CSG
            # tree. When it reaches a leaf (a Halfspace), it creates a <surface>
            # element for the corresponding surface if none has been created
            # thus far.
            def create_surface_elements(node, element):
                if isinstance(node, Halfspace):
                    path = "./surface[@id='{}']".format(node.surface.id)
                    if xml_element.find(path) is None:
                        xml_element.append(node.surface.to_xml_element())
                elif isinstance(node, Complement):
                    create_surface_elements(node.node, element)
                else:
                    for subnode in node:
                        create_surface_elements(subnode, element)

            # Call the recursive function from the top node
            create_surface_elements(self.region, xml_element)

        if self.temperature is not None:
            if isinstance(self.temperature, Iterable):
                element.set("temperature", ' '.join(
                    str(t) for t in self.temperature))
            else:
                element.set("temperature", str(self.temperature))

        if self.translation is not None:
            element.set("translation", ' '.join(map(str, self.translation)))

        if self.rotation is not None:
            element.set("rotation", ' '.join(map(str, self.rotation)))

        return element
