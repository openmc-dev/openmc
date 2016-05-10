from collections import OrderedDict, Iterable
from numbers import Real, Integral
from xml.etree import ElementTree as ET
import sys
import warnings

import openmc
import openmc.checkvalue as cv
from openmc.surface import Halfspace
from openmc.region import Region, Intersection, Complement

if sys.version_info[0] >= 3:
    basestring = str


# A static variable for auto-generated Cell IDs
AUTO_CELL_ID = 10000


def reset_auto_cell_id():
    global AUTO_CELL_ID
    AUTO_CELL_ID = 10000


class Cell(object):
    r"""A region of space defined as the intersection of half-space created by
    quadric surfaces.

    Parameters
    ----------
    cell_id : int, optional
        Unique identifier for the cell. If not specified, an identifier will
        automatically be assigned.
    name : str, optional
        Name of the cell. If not specified, the name is the empty string.
    fill : openmc.Material or openmc.Universe or openmc.Lattice or 'void' or iterable of openmc.Material, optional
        Indicates what the region of space is filled with
    region : openmc.Region, optional
        Region of space that is assigned to the cell.

    Attributes
    ----------
    id : int
        Unique identifier for the cell
    name : str
        Name of the cell
    fill : openmc.Material or openmc.Universe or openmc.Lattice or 'void' or iterable of openmc.Material
        Indicates what the region of space is filled with
    region : openmc.Region
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

    translation : Iterable of float
        If the cell is filled with a universe, this array specifies a vector
        that is used to translate (shift) the universe.
    offsets : ndarray
        Array of offsets used for distributed cell searches
    distribcell_index : int
        Index of this cell in distribcell arrays

    """

    def __init__(self, cell_id=None, name='', fill=None, region=None):
        # Initialize Cell class attributes
        self.id = cell_id
        self.name = name
        self._fill = None
        self._type = None
        self._region = None
        self._rotation = None
        self._translation = None
        self._offsets = None
        self._distribcell_index = None

        if fill is not None:
            self.fill = fill
        if region is not None:
            self.region = region

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
        string += '{0: <16}{1}{2}\n'.format('\tID', '=\t', self._id)
        string += '{0: <16}{1}{2}\n'.format('\tName', '=\t', self._name)

        if isinstance(self._fill, openmc.Material):
            string += '{0: <16}{1}{2}\n'.format('\tMaterial', '=\t',
                                                self._fill._id)
        elif isinstance(self._fill, Iterable):
            string += '{0: <16}{1}'.format('\tMaterial', '=\t')
            string += '['
            string += ', '.join(['void' if m == 'void' else str(m.id)
                                 for m in self.fill])
            string += ']\n'
        elif isinstance(self._fill, (openmc.Universe, openmc.Lattice)):
            string += '{0: <16}{1}{2}\n'.format('\tFill', '=\t',
                                                self._fill._id)
        else:
            string += '{0: <16}{1}{2}\n'.format('\tFill', '=\t', self._fill)

        string += '{0: <16}{1}{2}\n'.format('\tRegion', '=\t', self._region)

        string += '{0: <16}{1}{2}\n'.format('\tRotation', '=\t',
                                            self._rotation)
        string += '{0: <16}{1}{2}\n'.format('\tTranslation', '=\t',
                                            self._translation)
        string += '{0: <16}{1}{2}\n'.format('\tOffset', '=\t', self._offsets)
        string += '{0: <16}{1}{2}\n'.format('\tDistribcell index', '=\t',
                                            self._distribcell_index)

        return string

    @property
    def id(self):
        return self._id

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
        else:
            return None

    @property
    def region(self):
        return self._region

    @property
    def rotation(self):
        return self._rotation

    @property
    def translation(self):
        return self._translation

    @property
    def offsets(self):
        return self._offsets

    @property
    def distribcell_index(self):
        return self._distribcell_index

    @id.setter
    def id(self, cell_id):
        if cell_id is None:
            global AUTO_CELL_ID
            self._id = AUTO_CELL_ID
            AUTO_CELL_ID += 1
        else:
            cv.check_type('cell ID', cell_id, Integral)
            cv.check_greater_than('cell ID', cell_id, 0, equality=True)
            self._id = cell_id

    @name.setter
    def name(self, name):
        if name is not None:
            cv.check_type('cell name', name, basestring)
            self._name = name
        else:
            self._name = ''

    @fill.setter
    def fill(self, fill):
        if isinstance(fill, basestring):
            if fill.strip().lower() == 'void':
                self._type = 'void'
            else:
                msg = 'Unable to set Cell ID="{0}" to use a non-Material or ' \
                       'Universe fill "{1}"'.format(self._id, fill)
                raise ValueError(msg)

        elif isinstance(fill, openmc.Material):
            self._type = 'normal'

        elif isinstance(fill, Iterable):
            cv.check_type('cell.fill', fill, Iterable,
                          (openmc.Material, basestring))
            self._type = 'normal'

        elif isinstance(fill, openmc.Universe):
            self._type = 'fill'

        elif isinstance(fill, openmc.Lattice):
            self._type = 'lattice'

        else:
            msg = 'Unable to set Cell ID="{0}" to use a non-Material or ' \
                   'Universe fill "{1}"'.format(self._id, fill)
            raise ValueError(msg)

        self._fill = fill

    @rotation.setter
    def rotation(self, rotation):
        cv.check_type('cell rotation', rotation, Iterable, Real)
        cv.check_length('cell rotation', rotation, 3)
        self._rotation = rotation

    @translation.setter
    def translation(self, translation):
        cv.check_type('cell translation', translation, Iterable, Real)
        cv.check_length('cell translation', translation, 3)
        self._translation = translation

    @offsets.setter
    def offsets(self, offsets):
        cv.check_type('cell offsets', offsets, Iterable)
        self._offsets = offsets

    @region.setter
    def region(self, region):
        cv.check_type('cell region', region, Region)
        self._region = region

    @distribcell_index.setter
    def distribcell_index(self, ind):
        cv.check_type('distribcell index', ind, Integral)
        self._distribcell_index = ind

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
                self.region.nodes.append(region)
            else:
                self.region = Intersection(self.region, region)

    def get_cell_instance(self, path, distribcell_index):

        # If the Cell is filled by a Material
        if self._type == 'normal' or self._type == 'void':
            offset = 0

        # If the Cell is filled by a Universe
        elif self._type == 'fill':
            offset = self.offsets[distribcell_index-1]
            offset += self.fill.get_cell_instance(path, distribcell_index)

        # If the Cell is filled by a Lattice
        else:
            offset = self.fill.get_cell_instance(path, distribcell_index)

        return offset

    def get_all_nuclides(self):
        """Return all nuclides contained in the cell

        Returns
        -------
        nuclides : dict
            Dictionary whose keys are nuclide names and values are 2-tuples of
            (nuclide, density)

        """

        nuclides = OrderedDict()

        if self._type != 'void':
            nuclides.update(self._fill.get_all_nuclides())

        return nuclides

    def get_all_cells(self):
        """Return all cells that are contained within this one if it is filled with a
        universe or lattice

        Returns
        -------
        cells : dict
            Dictionary whose keys are cell IDs and values are :class:`Cell`
            instances

        """

        cells = OrderedDict()

        if self._type == 'fill' or self._type == 'lattice':
            cells.update(self._fill.get_all_cells())

        return cells

    def get_all_materials(self):
        """Return all materials that are contained within the cell

        Returns
        -------
        materials : dict
            Dictionary whose keys are material IDs and values are
            :class:`Material` instances

        """

        materials = OrderedDict()
        if self.fill_type == 'material':
            materials[self.fill.id] = self.fill

        # Append all Cells in each Cell in the Universe to the dictionary
        cells = self.get_all_cells()
        for cell_id, cell in cells.items():
            materials.update(cell.get_all_materials())

        return materials

    def get_all_universes(self):
        """Return all universes that are contained within this one if any of
        its cells are filled with a universe or lattice.

        Returns
        -------
        universes : dict
            Dictionary whose keys are universe IDs and values are
            :class:`Universe` instances

        """

        universes = OrderedDict()

        if self._type == 'fill':
            universes[self._fill._id] = self._fill
            universes.update(self._fill.get_all_universes())
        elif self._type == 'lattice':
            universes.update(self._fill.get_all_universes())

        return universes

    def create_xml_subelement(self, xml_element):
        element = ET.Element("cell")
        element.set("id", str(self.id))

        if len(self._name) > 0:
            element.set("name", str(self.name))

        if isinstance(self.fill, basestring):
            element.set("material", "void")

        elif isinstance(self.fill, openmc.Material):
            element.set("material", str(self.fill.id))

        elif isinstance(self.fill, Iterable):
            element.set("material", ' '.join([m if m == 'void' else str(m.id)
                                              for m in self.fill]))

        elif isinstance(self.fill, (openmc.Universe, openmc.Lattice)):
            element.set("fill", str(self.fill.id))
            self.fill.create_xml_subelement(xml_element)

        else:
            element.set("fill", str(self.fill))
            self.fill.create_xml_subelement(xml_element)

        if self.region is not None:
            # Set the region attribute with the region specification
            element.set("region", str(self.region))

            # Only surfaces that appear in a region are added to the geometry
            # file, so the appropriate check is performed here. First we create
            # a function which is called recursively to navigate through the CSG
            # tree. When it reaches a leaf (a Halfspace), it creates a <surface>
            # element for the corresponding surface if none has been created
            # thus far.
            def create_surface_elements(node, element):
                if isinstance(node, Halfspace):
                    path = './surface[@id=\'{0}\']'.format(node.surface.id)
                    if xml_element.find(path) is None:
                        surface_subelement = node.surface.create_xml_subelement()
                        xml_element.append(surface_subelement)
                elif isinstance(node, Complement):
                    create_surface_elements(node.node, element)
                else:
                    for subnode in node.nodes:
                        create_surface_elements(subnode, element)

            # Call the recursive function from the top node
            create_surface_elements(self.region, xml_element)

        if self.translation is not None:
            element.set("translation", ' '.join(map(str, self.translation)))

        if self.rotation is not None:
            element.set("rotation", ' '.join(map(str, self.rotation)))

        return element
