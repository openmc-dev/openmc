from collections import Iterable, OrderedDict
from copy import deepcopy
from numbers import Real, Integral
import warnings
from xml.etree import ElementTree as ET
import sys
if sys.version_info[0] >= 3:
    basestring = str

import openmc
import openmc.checkvalue as cv
from openmc.clean_xml import *
from openmc.data import natural_abundance


# A static variable for auto-generated Material IDs
AUTO_MATERIAL_ID = 10000


def reset_auto_material_id():
    global AUTO_MATERIAL_ID
    AUTO_MATERIAL_ID = 10000


# Units for density supported by OpenMC
DENSITY_UNITS = ['g/cm3', 'g/cc', 'kg/cm3', 'atom/b-cm', 'atom/cm3', 'sum',
                 'macro']


class Material(object):
    """A material composed of a collection of nuclides/elements that can be
    assigned to a region of space.

    Parameters
    ----------
    material_id : int, optional
        Unique identifier for the material. If not specified, an identifier will
        automatically be assigned.
    name : str, optional
        Name of the material. If not specified, the name will be the empty
        string.

    Attributes
    ----------
    id : int
        Unique identifier for the material
    density : float
        Density of the material (units defined separately)
    density_units : str
        Units used for `density`. Can be one of 'g/cm3', 'g/cc', 'kg/cm3',
        'atom/b-cm', 'atom/cm3', 'sum', or 'macro'.  The 'macro' unit only
        applies in the case of a multi-group calculation.
    elements : collections.OrderedDict
        Dictionary whose keys are element names and values are 3-tuples
        consisting of an :class:`openmc.Element` instance, the percent density,
        and the percent type (atom or weight fraction).
    nuclides : collections.OrderedDict
        Dictionary whose keys are nuclide names and values are 3-tuples
        consisting of an :class:`openmc.Nuclide` instance, the percent density,
        and the percent type (atom or weight fraction).

    """

    def __init__(self, material_id=None, name=''):
        # Initialize class attributes
        self.id = material_id
        self.name = name
        self._density = None
        self._density_units = ''

        # An ordered dictionary of Nuclides (order affects OpenMC results)
        # Keys         - Nuclide names
        # Values     - tuple (nuclide, percent, percent type)
        self._nuclides = OrderedDict()

        # The single instance of Macroscopic data present in this material
        # (only one is allowed, hence this is different than _nuclides, etc)
        self._macroscopic = None

        # An ordered dictionary of Elements (order affects OpenMC results)
        # Keys         - Element names
        # Values     - tuple (element, percent, percent type)
        self._elements = OrderedDict()

        # If specified, a list of tuples of (table name, xs identifier)
        self._sab = []

        # If true, the material will be initialized as distributed
        self._convert_to_distrib_comps = False

        # If specified, this file will be used instead of composition values
        self._distrib_otf_file = None

    def __eq__(self, other):
        if not isinstance(other, Material):
            return False
        elif self.id != other.id:
            return False
        elif self.name != other.name:
            return False
        # FIXME: We cannot compare densities since OpenMC outputs densities
        # in atom/b-cm in summary.h5 irregardless of input units, and we
        # cannot compute the sum percent in Python since we lack AWR
        #elif self.density != other.density:
        #    return False
        #elif self._nuclides != other._nuclides:
        #    return False
        #elif self._elements != other._elements:
        #   return False
        elif self._sab != other._sab:
            return False
        else:
            return True

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(repr(self))

    def __repr__(self):
        string = 'Material\n'
        string += '{0: <16}{1}{2}\n'.format('\tID', '=\t', self._id)
        string += '{0: <16}{1}{2}\n'.format('\tName', '=\t', self._name)

        string += '{0: <16}{1}{2}'.format('\tDensity', '=\t', self._density)
        string += ' [{0}]\n'.format(self._density_units)

        string += '{0: <16}\n'.format('\tS(a,b) Tables')

        for sab in self._sab:
            string += '{0: <16}{1}[{2}{3}]\n'.format('\tS(a,b)', '=\t',
                                                     sab[0], sab[1])

        string += '{0: <16}\n'.format('\tNuclides')

        for nuclide in self._nuclides:
            percent = self._nuclides[nuclide][1]
            percent_type = self._nuclides[nuclide][2]
            string += '{0: <16}'.format('\t{0}'.format(nuclide))
            string += '=\t{0: <12} [{1}]\n'.format(percent, percent_type)

        if self._macroscopic is not None:
            string += '{0: <16}\n'.format('\tMacroscopic Data')
            string += '{0: <16}'.format('\t{0}'.format(self._macroscopic))

        string += '{0: <16}\n'.format('\tElements')

        for element in self._elements:
            percent = self._elements[element][1]
            percent_type = self._elements[element][2]
            string += '{0: <16}'.format('\t{0}'.format(element))
            string += '=\t{0: <12} [{1}]\n'.format(percent, percent_type)

        return string

    def __deepcopy__(self, memo):
        existing = memo.get(id(self))

        if existing is None:
            # If this is the first time we have tried to copy this object, create a copy
            clone = type(self).__new__(type(self))
            clone._id = self._id
            clone._name = self._name
            clone._density = self._density
            clone._density_units = self._density_units
            clone._nuclides = deepcopy(self._nuclides, memo)
            clone._macroscopic = self._macroscopic
            clone._elements = deepcopy(self._elements, memo)
            clone._sab = deepcopy(self._sab, memo)
            clone._convert_to_distrib_comps = self._convert_to_distrib_comps
            clone._distrib_otf_file = self._distrib_otf_file

            memo[id(self)] = clone

            return clone

        else:
            # If this object has been copied before, return the first copy made
            return existing

    @property
    def id(self):
        return self._id

    @property
    def name(self):
        return self._name

    @property
    def density(self):
        return self._density

    @property
    def density_units(self):
        return self._density_units

    @property
    def elements(self):
        return self._elements

    @property
    def nuclides(self):
        return self._nuclides

    @property
    def convert_to_distrib_comps(self):
        return self._convert_to_distrib_comps

    @property
    def distrib_otf_file(self):
        return self._distrib_otf_file

    @id.setter
    def id(self, material_id):

        if material_id is None:
            global AUTO_MATERIAL_ID
            self._id = AUTO_MATERIAL_ID
            AUTO_MATERIAL_ID += 1
        else:
            cv.check_type('material ID', material_id, Integral)
            cv.check_greater_than('material ID', material_id, 0, equality=True)
            self._id = material_id

    @name.setter
    def name(self, name):
        if name is not None:
            cv.check_type('name for Material ID="{0}"'.format(self._id),
                          name, basestring)
            self._name = name
        else:
            self._name = ''

    def set_density(self, units, density=None):
        """Set the density of the material

        Parameters
        ----------
        units : {'g/cm3', 'g/cc', 'km/cm3', 'atom/b-cm', 'atom/cm3', 'sum', 'macro'}
            Physical units of density.
        density : float, optional
            Value of the density. Must be specified unless units is given as
            'sum'.

        """

        cv.check_value('density units', units, DENSITY_UNITS)
        self._density_units = units

        if units is 'sum':
            if density is not None:
                msg = 'Density "{0}" for Material ID="{1}" is ignored ' \
                      'because the unit is "sum"'.format(density, self.id)
                warnings.warn(msg)
        else:
            if density is None:
                msg = 'Unable to set the density for Material ID="{0}" ' \
                      'because a density value must be given when not using ' \
                      '"sum" unit'.format(self.id)
                raise ValueError(msg)

            cv.check_type('the density for Material ID="{0}"'.format(self.id),
                          density, Real)
            self._density = density

    @distrib_otf_file.setter
    def distrib_otf_file(self, filename):
        # TODO: remove this when distributed materials are merged
        warnings.warn('This feature is not yet implemented in a release '
                      'version of openmc')

        if not isinstance(filename, basestring) and filename is not None:
            msg = 'Unable to add OTF material file to Material ID="{0}" with a ' \
                  'non-string name "{1}"'.format(self._id, filename)
            raise ValueError(msg)

        self._distrib_otf_file = filename

    @convert_to_distrib_comps.setter
    def convert_to_distrib_comps(self):
        # TODO: remove this when distributed materials are merged
        warnings.warn('This feature is not yet implemented in a release '
                      'version of openmc')

        self._convert_to_distrib_comps = True

    def add_nuclide(self, nuclide, percent, percent_type='ao'):
        """Add a nuclide to the material

        Parameters
        ----------
        nuclide : str or openmc.Nuclide
            Nuclide to add
        percent : float
            Atom or weight percent
        percent_type : {'ao', 'wo'}
            'ao' for atom percent and 'wo' for weight percent

        """

        if self._macroscopic is not None:
            msg = 'Unable to add a Nuclide to Material ID="{0}" as a ' \
                  'macroscopic data-set has already been added'.format(self._id)
            raise ValueError(msg)

        if not isinstance(nuclide, (openmc.Nuclide, basestring)):
            msg = 'Unable to add a Nuclide to Material ID="{0}" with a ' \
                  'non-Nuclide value "{1}"'.format(self._id, nuclide)
            raise ValueError(msg)

        elif not isinstance(percent, Real):
            msg = 'Unable to add a Nuclide to Material ID="{0}" with a ' \
                  'non-floating point value "{1}"'.format(self._id, percent)
            raise ValueError(msg)

        elif percent_type not in ['ao', 'wo', 'at/g-cm']:
            msg = 'Unable to add a Nuclide to Material ID="{0}" with a ' \
                  'percent type "{1}"'.format(self._id, percent_type)
            raise ValueError(msg)

        if isinstance(nuclide, openmc.Nuclide):
            # Copy this Nuclide to separate it from the Nuclide in
            # other Materials
            nuclide = deepcopy(nuclide)
        else:
            nuclide = openmc.Nuclide(nuclide)

        self._nuclides[nuclide._name] = (nuclide, percent, percent_type)

    def remove_nuclide(self, nuclide):
        """Remove a nuclide from the material

        Parameters
        ----------
        nuclide : openmc.Nuclide
            Nuclide to remove

        """

        if not isinstance(nuclide, openmc.Nuclide):
            msg = 'Unable to remove a Nuclide "{0}" in Material ID="{1}" ' \
                  'since it is not a Nuclide'.format(self._id, nuclide)
            raise ValueError(msg)

        # If the Material contains the Nuclide, delete it
        if nuclide._name in self._nuclides:
            del self._nuclides[nuclide._name]

    def add_macroscopic(self, macroscopic):
        """Add a macroscopic to the material.  This will also set the
        density of the material to 1.0, unless it has been otherwise set,
        as a default for Macroscopic cross sections.

        Parameters
        ----------
        macroscopic : str or openmc.Macroscopic
            Macroscopic to add

        """

        # Ensure no nuclides, elements, or sab are added since these would be
        # incompatible with macroscopics
        if self._nuclides or self._elements or self._sab:
            msg = 'Unable to add a Macroscopic data set to Material ID="{0}" ' \
                  'with a macroscopic value "{1}" as an incompatible data ' \
                  'member (i.e., nuclide, element, or S(a,b) table) ' \
                  'has already been added'.format(self._id, macroscopic)
            raise ValueError(msg)

        if not isinstance(macroscopic, (openmc.Macroscopic, basestring)):
            msg = 'Unable to add a Macroscopic to Material ID="{0}" with a ' \
                  'non-Macroscopic value "{1}"'.format(self._id, macroscopic)
            raise ValueError(msg)

        if isinstance(macroscopic, openmc.Macroscopic):
            # Copy this Macroscopic to separate it from the Macroscopic in
            # other Materials
            macroscopic = deepcopy(macroscopic)
        else:
            macroscopic = openmc.Macroscopic(macroscopic)

        if self._macroscopic is None:
            self._macroscopic = macroscopic
        else:
            msg = 'Unable to add a Macroscopic to Material ID="{0}", ' \
                  'Only One Macroscopic allowed per ' \
                  'Material!'.format(self._id, macroscopic)
            raise ValueError(msg)

        # Generally speaking, the density for a macroscopic object will
        # be 1.0.  Therefore, lets set density to 1.0 so that the user
        # doesnt need to set it unless its needed.
        # Of course, if the user has already set a value of density,
        # then we will not override it.
        if self._density is None:
            self.set_density('macro', 1.0)

    def remove_macroscopic(self, macroscopic):
        """Remove a macroscopic from the material

        Parameters
        ----------
        macroscopic : openmc.Macroscopic
            Macroscopic to remove

        """

        if not isinstance(macroscopic, openmc.Macroscopic):
            msg = 'Unable to remove a Macroscopic "{0}" in Material ID="{1}" ' \
                  'since it is not a Macroscopic'.format(self._id, macroscopic)
            raise ValueError(msg)

        # If the Material contains the Macroscopic, delete it
        if macroscopic._name == self._macroscopic.name:
            self._macroscopic = None

    def add_element(self, element, percent, percent_type='ao', expand=False):
        """Add a natural element to the material

        Parameters
        ----------
        element : openmc.Element or str
            Element to add
        percent : float
            Atom or weight percent
        percent_type : {'ao', 'wo'}, optional
            'ao' for atom percent and 'wo' for weight percent. Defaults to atom
            percent.
        expand : bool, optional
            Whether to expand the natural element into its naturally-occurring
            isotopes. Defaults to False.

        """

        if self._macroscopic is not None:
            msg = 'Unable to add an Element to Material ID="{0}" as a ' \
                  'macroscopic data-set has already been added'.format(self._id)
            raise ValueError(msg)

        if not isinstance(element, (openmc.Element, basestring)):
            msg = 'Unable to add an Element to Material ID="{0}" with a ' \
                  'non-Element value "{1}"'.format(self._id, element)
            raise ValueError(msg)

        if not isinstance(percent, Real):
            msg = 'Unable to add an Element to Material ID="{0}" with a ' \
                  'non-floating point value "{1}"'.format(self._id, percent)
            raise ValueError(msg)

        if percent_type not in ['ao', 'wo']:
            msg = 'Unable to add an Element to Material ID="{0}" with a ' \
                  'percent type "{1}"'.format(self._id, percent_type)
            raise ValueError(msg)

        # Copy this Element to separate it from same Element in other Materials
        if isinstance(element, openmc.Element):
            element = deepcopy(element)
        else:
            element = openmc.Element(element)

        if expand:
            if percent_type == 'wo':
                raise NotImplementedError('Expanding natural element based on '
                                          'weight percent is not yet supported.')
            for isotope, abundance in element.expand():
                self._nuclides[isotope.name] = (
                    isotope, percent*abundance, percent_type)
        else:
            self._elements[element.name] = (element, percent, percent_type)

    def remove_element(self, element):
        """Remove a natural element from the material

        Parameters
        ----------
        element : openmc.Element
            Element to remove

        """

        # If the Material contains the Element, delete it
        if element._name in self._elements:
            del self._elements[element._name]

    def add_s_alpha_beta(self, name, xs):
        r"""Add an :math:`S(\alpha,\beta)` table to the material

        Parameters
        ----------
        name : str
            Name of the :math:`S(\alpha,\beta)` table
        xs : str
            Cross section identifier, e.g. '71t'

        """

        if self._macroscopic is not None:
            msg = 'Unable to add an S(a,b) table to Material ID="{0}" as a ' \
                  'macroscopic data-set has already been added'.format(self._id)
            raise ValueError(msg)

        if not isinstance(name, basestring):
            msg = 'Unable to add an S(a,b) table to Material ID="{0}" with a ' \
                        'non-string table name "{1}"'.format(self._id, name)
            raise ValueError(msg)

        if not isinstance(xs, basestring):
            msg = 'Unable to add an S(a,b) table to Material ID="{0}" with a ' \
                  'non-string cross-section identifier "{1}"'.format(self._id, xs)
            raise ValueError(msg)

        self._sab.append((name, xs))

    def make_isotropic_in_lab(self):
        for nuclide_name in self._nuclides:
            self._nuclides[nuclide_name][0].scattering = 'iso-in-lab'
        for element_name in self._elements:
            self._elements[element_name][0].scattering = 'iso-in-lab'

    def get_all_nuclides(self):
        """Returns all nuclides in the material

        Returns
        -------
        nuclides : dict
            Dictionary whose keys are nuclide names and values are 2-tuples of
            (nuclide, density)

        """

        nuclides = OrderedDict()

        for nuclide_name, nuclide_tuple in self._nuclides.items():
            nuclide = nuclide_tuple[0]
            density = nuclide_tuple[1]
            nuclides[nuclide._name] = (nuclide, density)

        for element_name, element_tuple in self._elements.items():
            element = element_tuple[0]
            density = element_tuple[1]

            # Expand natural element into isotopes
            for isotope, abundance in element.expand():
                nuclides[isotope.name] = (isotope, density*abundance)

        return nuclides

    def _get_nuclide_xml(self, nuclide, distrib=False):
        xml_element = ET.Element("nuclide")
        xml_element.set("name", nuclide[0]._name)

        if not distrib:
            if nuclide[2] == 'ao':
                xml_element.set("ao", str(nuclide[1]))
            else:
                xml_element.set("wo", str(nuclide[1]))

        if nuclide[0].xs is not None:
            xml_element.set("xs", nuclide[0].xs)

        if not nuclide[0].scattering is None:
            xml_element.set("scattering", nuclide[0].scattering)

        return xml_element

    def _get_macroscopic_xml(self, macroscopic):
        xml_element = ET.Element("macroscopic")
        xml_element.set("name", macroscopic._name)

        if macroscopic.xs is not None:
            xml_element.set("xs", macroscopic.xs)

        return xml_element

    def _get_element_xml(self, element, distrib=False):
        xml_element = ET.Element("element")
        xml_element.set("name", str(element[0]._name))

        if not distrib:
            if element[2] == 'ao':
                xml_element.set("ao", str(element[1]))
            else:
                xml_element.set("wo", str(element[1]))

        if element[0].xs is not None:
            xml_element.set("xs", element[0].xs)

        if not element[0].scattering is None:
            xml_element.set("scattering", element[0].scattering)

        return xml_element

    def _get_nuclides_xml(self, nuclides, distrib=False):
        xml_elements = []

        for nuclide in nuclides.values():
            xml_elements.append(self._get_nuclide_xml(nuclide, distrib))

        return xml_elements

    def _get_elements_xml(self, elements, distrib=False):
        xml_elements = []

        for element in elements.values():
            xml_elements.append(self._get_element_xml(element, distrib))

        return xml_elements

    def get_material_xml(self):
        """Return XML representation of the material

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing material data

        """

        # Create Material XML element
        element = ET.Element("material")
        element.set("id", str(self._id))

        if len(self._name) > 0:
            element.set("name", str(self._name))

        # Create density XML subelement
        subelement = ET.SubElement(element, "density")
        if self._density_units is not 'sum':
            subelement.set("value", str(self._density))
        subelement.set("units", self._density_units)

        if not self._convert_to_distrib_comps:
            if self._macroscopic is None:
                # Create nuclide XML subelements
                subelements = self._get_nuclides_xml(self._nuclides)
                for subelement in subelements:
                    element.append(subelement)

                # Create element XML subelements
                subelements = self._get_elements_xml(self._elements)
                for subelement in subelements:
                    element.append(subelement)
            else:
                # Create macroscopic XML subelements
                subelement = self._get_macroscopic_xml(self._macroscopic)
                element.append(subelement)

        else:
            subelement = ET.SubElement(element, "compositions")

            comps = []
            allnucs = self._nuclides.values() + self._elements.values()
            dist_per_type = allnucs[0][2]
            for nuc, per, typ in allnucs:
                if not typ == dist_per_type:
                    msg = 'All nuclides and elements in a distributed ' \
                          'material must have the same type, either ao or wo'
                    raise ValueError(msg)
                comps.append(per)

            if self._distrib_otf_file is None:
                # Create values and units subelements
                subsubelement = ET.SubElement(subelement, "values")
                subsubelement.text = ' '.join([str(c) for c in comps])
                subsubelement = ET.SubElement(subelement, "units")
                subsubelement.text = dist_per_type
            else:
                # Specify the materials file
                subsubelement = ET.SubElement(subelement, "otf_file_path")
                subsubelement.text = self._distrib_otf_file

            if self._macroscopic is None:
                # Create nuclide XML subelements
                subelements = self.get_nuclides_xml(self._nuclides, distrib=True)
                for subelement_nuc in subelements:
                    subelement.append(subelement_nuc)

                # Create element XML subelements
                subelements = self._get_elements_xml(self._elements, distrib=True)
                for subsubelement in subelements:
                    subelement.append(subsubelement)
            else:
                # Create macroscopic XML subelements
                subsubelement = self._get_macroscopic_xml(self._macroscopic,
                                                          distrib=True)
                subelement.append(subsubelement)

        if len(self._sab) > 0:
            for sab in self._sab:
                subelement = ET.SubElement(element, "sab")
                subelement.set("name", sab[0])
                subelement.set("xs", sab[1])

        return element


class Materials(cv.CheckedList):
    """Collection of Materials used for an OpenMC simulation.

    This class corresponds directly to the materials.xml input file. It can be
    thought of as a normal Python list where each member is a
    :class:`Material`. It behaves like a list as the following example
    demonstrates:

    >>> fuel = openmc.Material()
    >>> clad = openmc.Material()
    >>> water = openmc.Material()
    >>> m = openmc.Materials([fuel])
    >>> m.append(water)
    >>> m += [clad]

    Parameters
    ----------
    materials : Iterable of openmc.Material
        Materials to add to the collection

    Attributes
    ----------
    default_xs : str
        The default cross section identifier applied to a nuclide when none is
        specified

    """

    def __init__(self, materials=None):
        super(Materials, self).__init__(Material, 'materials collection')
        self._default_xs = None
        self._materials_file = ET.Element("materials")
        if materials is not None:
            self += materials

    @property
    def default_xs(self):
        return self._default_xs

    @default_xs.setter
    def default_xs(self, xs):
        cv.check_type('default xs', xs, basestring)
        self._default_xs = xs

    def add_material(self, material):
        """Append material to collection

        .. deprecated:: 0.8
            Use :meth:`Materials.append` instead.

        Parameters
        ----------
        material : openmc.Material
            Material to add

        """
        warnings.warn("Materials.add_material(...) has been deprecated and may be "
                      "removed in a future version. Use Material.append(...) "
                      "instead.", DeprecationWarning)
        self.append(material)

    def add_materials(self, materials):
        """Add multiple materials to the collection

        .. deprecated:: 0.8
            Use compound assignment instead.

        Parameters
        ----------
        materials : Iterable of openmc.Material
            Materials to add

        """
        warnings.warn("Materials.add_materials(...) has been deprecated and may be "
                      "removed in a future version. Use compound assignment "
                      "instead.", DeprecationWarning)
        for material in materials:
            self.append(material)

    def append(self, material):
        """Append material to collection

        Parameters
        ----------
        material : openmc.Material
            Material to append

        """
        super(Materials, self).append(material)

    def insert(self, index, material):
        """Insert material before index

        Parameters
        ----------
        index : int
            Index in list
        material : openmc.Material
            Material to insert

        """
        super(Materials, self).insert(index, material)

    def remove_material(self, material):
        """Remove a material from the file

        .. deprecated:: 0.8
            Use :meth:`Materials.remove` instead.

        Parameters
        ----------
        material : openmc.Material
            Material to remove

        """
        warnings.warn("Materials.remove_material(...) has been deprecated and "
                      "may be removed in a future version. Use "
                      "Materials.remove(...) instead.", DeprecationWarning)
        self.remove(material)

    def make_isotropic_in_lab(self):
        for material in self:
            material.make_isotropic_in_lab()

    def _create_material_subelements(self):
        if self._default_xs is not None:
            subelement = ET.SubElement(self._materials_file, "default_xs")
            subelement.text = self._default_xs

        for material in self:
            xml_element = material.get_material_xml()
            self._materials_file.append(xml_element)

    def export_to_xml(self):
        """Create a materials.xml file that can be used for a simulation.

        """

        # Reset xml element tree
        self._materials_file.clear()

        self._create_material_subelements()

        # Clean the indentation in the file to be user-readable
        sort_xml_elements(self._materials_file)
        clean_xml_indentation(self._materials_file)

        # Write the XML Tree to the materials.xml file
        tree = ET.ElementTree(self._materials_file)
        tree.write("materials.xml", xml_declaration=True,
                             encoding='utf-8', method="xml")
