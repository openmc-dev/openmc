from collections import OrderedDict
from copy import deepcopy
from numbers import Real, Integral
import warnings
from xml.etree import ElementTree as ET
import sys

from six import string_types

import openmc
import openmc.data
import openmc.checkvalue as cv
from openmc.clean_xml import sort_xml_elements, clean_xml_indentation


# A static variable for auto-generated Material IDs
AUTO_MATERIAL_ID = 10000


def reset_auto_material_id():
    """Reset counter for auto-generated material IDs."""
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
    temperature : float, optional
        Temperature of the material in Kelvin. If not specified, the material
        inherits the default temperature applied to the model.

    Attributes
    ----------
    id : int
        Unique identifier for the material
    temperature : float
        Temperature of the material in Kelvin.
    density : float
        Density of the material (units defined separately)
    density_units : str
        Units used for `density`. Can be one of 'g/cm3', 'g/cc', 'kg/cm3',
        'atom/b-cm', 'atom/cm3', 'sum', or 'macro'.  The 'macro' unit only
        applies in the case of a multi-group calculation.
    elements : list of tuple
        List in which each item is a 4-tuple consisting of an
        :class:`openmc.Element` instance, the percent density, the percent
        type ('ao' or 'wo'), and enrichment.
    nuclides : list of tuple
        List in which each item is a 3-tuple consisting of an
        :class:`openmc.Nuclide` instance, the percent density, and the percent
        type ('ao' or 'wo').

    """

    def __init__(self, material_id=None, name='', temperature=None):
        # Initialize class attributes
        self.id = material_id
        self.name = name
        self.temperature = temperature
        self._density = None
        self._density_units = ''

        # A list of tuples (nuclide, percent, percent type)
        self._nuclides = []

        # The single instance of Macroscopic data present in this material
        # (only one is allowed, hence this is different than _nuclides, etc)
        self._macroscopic = None

        # A list of tuples (element, percent, percent type, enrichment)
        self._elements = []

        # If specified, a list of table names
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
        string += '{0: <16}{1}{2}\n'.format('\Temperature', '=\t',
                                            self._temperature)

        string += '{0: <16}{1}{2}'.format('\tDensity', '=\t', self._density)
        string += ' [{0}]\n'.format(self._density_units)

        string += '{0: <16}\n'.format('\tS(a,b) Tables')

        for sab in self._sab:
            string += '{0: <16}{1}{2}\n'.format('\tS(a,b)', '=\t', sab)

        string += '{0: <16}\n'.format('\tNuclides')

        for nuclide, percent, percent_type in self._nuclides:
            string += '{0: <16}'.format('\t{0.name}'.format(nuclide))
            string += '=\t{0: <12} [{1}]\n'.format(percent, percent_type)

        if self._macroscopic is not None:
            string += '{0: <16}\n'.format('\tMacroscopic Data')
            string += '{0: <16}'.format('\t{0}'.format(self._macroscopic))

        string += '{0: <16}\n'.format('\tElements')

        for element, percent, percent_type, enr in self._elements:
            string += '{0: <16}'.format('\t{0.name}'.format(element))
            if enr is None:
                string += '=\t{0: <12} [{1}]\n'.format(percent, percent_type)
            else:
                string += '=\t{0: <12} [{1}] @ {2} w/o enrichment\n'\
                          .format(percent, percent_type, enr)

        return string

    @property
    def id(self):
        return self._id

    @property
    def name(self):
        return self._name

    @property
    def temperature(self):
        return self._temperature

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
                          name, string_types)
            self._name = name
        else:
            self._name = ''

    @temperature.setter
    def temperature(self, temperature):
        cv.check_type('Temperature for Material ID="{0}"'.format(self._id),
                      temperature, (Real, type(None)))
        self._temperature = temperature

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

        if not isinstance(filename, string_types) and filename is not None:
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

        if not isinstance(nuclide, string_types + (openmc.Nuclide,)):
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

        self._nuclides.append((nuclide, percent, percent_type))

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
        for nuc in self._nuclides:
            if nuclide == nuc:
                self._nuclides.remove(nuc)

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

        if not isinstance(macroscopic, string_types + (openmc.Macroscopic,)):
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
            msg = 'Unable to add a Macroscopic to Material ID="{0}". ' \
                  'Only one Macroscopic allowed per ' \
                  'Material.'.format(self._id)
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
        if macroscopic.name == self._macroscopic.name:
            self._macroscopic = None

    def add_element(self, element, percent, percent_type='ao', enrichment=None):
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
        enrichment : float, optional
            Enrichment for U235 in weight percent. For example, input 4.95 for
            4.95 weight percent enriched U. Default is None
            (natural composition).

        """

        if self._macroscopic is not None:
            msg = 'Unable to add an Element to Material ID="{0}" as a ' \
                  'macroscopic data-set has already been added'.format(self._id)
            raise ValueError(msg)

        if not isinstance(element, string_types + (openmc.Element,)):
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

        if enrichment is not None:
            if not isinstance(enrichment, Real):
                msg = 'Unable to add an Element to Material ID="{0}" with a ' \
                      'non-floating point enrichment value "{1}"'\
                      .format(self._id, enrichment)
                raise ValueError(msg)

            elif element.name != 'U':
                msg = 'Unable to use enrichment for element {0} which is not ' \
                      'uranium for Material ID="{1}"'.format(element.name,
                                                             self._id)
                raise ValueError(msg)

            # Check that the enrichment is in the valid range
            cv.check_less_than('enrichment', enrichment, 100./1.008)
            cv.check_greater_than('enrichment', enrichment, 0., equality=True)

            if enrichment > 5.0:
                msg = 'A uranium enrichment of {0} was given for Material ID='\
                      '"{1}". OpenMC assumes the U234/U235 mass ratio is '\
                      'constant at 0.008, which is only valid at low ' \
                      'enrichments. Consider setting the isotopic ' \
                      'composition manually for enrichments over 5%.'.\
                      format(enrichment, self._id)
                warnings.warn(msg)

        self._elements.append((element, percent, percent_type, enrichment))

    def remove_element(self, element):
        """Remove a natural element from the material

        Parameters
        ----------
        element : openmc.Element
            Element to remove

        """

        if not isinstance(element, openmc.Element):
            msg = 'Unable to remove "{0}" in Material ID="{1}" ' \
                  'since it is not an Element'.format(self.id, element)
            raise ValueError(msg)

        # If the Material contains the Element, delete it
        for elm in self._elements:
            if element == elm:
                self._elements.remove(elm)

    def add_s_alpha_beta(self, name):
        r"""Add an :math:`S(\alpha,\beta)` table to the material

        Parameters
        ----------
        name : str
            Name of the :math:`S(\alpha,\beta)` table

        """

        if self._macroscopic is not None:
            msg = 'Unable to add an S(a,b) table to Material ID="{0}" as a ' \
                  'macroscopic data-set has already been added'.format(self._id)
            raise ValueError(msg)

        if not isinstance(name, string_types):
            msg = 'Unable to add an S(a,b) table to Material ID="{0}" with a ' \
                        'non-string table name "{1}"'.format(self._id, name)
            raise ValueError(msg)

        new_name = openmc.data.get_thermal_name(name)
        if new_name != name:
            msg = 'OpenMC S(a,b) tables follow the GND naming convention. ' \
                  'Table "{}" is being renamed as "{}".'.format(name, new_name)
            warnings.warn(msg)

        self._sab.append(new_name)

    def make_isotropic_in_lab(self):
        for nuclide, percent, percent_type in self._nuclides:
            nuclide.scattering = 'iso-in-lab'
        for element, percent, percent_type in self._elements:
            element.scattering = 'iso-in-lab'

    def get_nuclides(self):
        """Returns all nuclides in the material

        Returns
        -------
        nuclides : list of str
            List of nuclide names

        """

        nuclides = []

        for nuclide, percent, percent_type in self._nuclides:
            nuclides.append(nuclide.name)

        for ele, ele_pct, ele_pct_type, enr in self._elements:

            # Expand natural element into isotopes
            isotopes = ele.expand(ele_pct, ele_pct_type, enr)
            for iso, iso_pct, iso_pct_type in isotopes:
                nuclides.append(iso.name)

        return nuclides

    def get_nuclide_densities(self):
        """Returns all nuclides in the material and their densities

        Returns
        -------
        nuclides : dict
            Dictionary whose keys are nuclide names and values are 3-tuples of
            (nuclide, density percent, density percent type)

        """

        nuclides = OrderedDict()

        for nuclide, density, density_type in self._nuclides:
            nuclides[nuclide.name] = (nuclide, density)

        for ele, ele_pct, ele_pct_type, enr in self._elements:

            # Expand natural element into isotopes
            isotopes = ele.expand(ele_pct, ele_pct_type, enr)
            for iso, iso_pct, iso_pct_type in isotopes:
                nuclides[iso.name] = (iso, iso_pct, iso_pct_type)

        return nuclides

    def _get_nuclide_xml(self, nuclide, distrib=False):
        xml_element = ET.Element("nuclide")
        xml_element.set("name", nuclide[0].name)

        if not distrib:
            if nuclide[2] == 'ao':
                xml_element.set("ao", str(nuclide[1]))
            else:
                xml_element.set("wo", str(nuclide[1]))

        if not nuclide[0].scattering is None:
            xml_element.set("scattering", nuclide[0].scattering)

        return xml_element

    def _get_macroscopic_xml(self, macroscopic):
        xml_element = ET.Element("macroscopic")
        xml_element.set("name", macroscopic.name)

        return xml_element

    def _get_element_xml(self, element, cross_sections, distrib=False):

        # Get the nuclides in this element
        nuclides = element[0].expand(element[1], element[2], element[3],
                                     cross_sections)

        xml_elements = []
        for nuclide in nuclides:
            xml_elements.append(self._get_nuclide_xml(nuclide, distrib))

        return xml_elements

    def _get_nuclides_xml(self, nuclides, distrib=False):

        xml_elements = []

        for nuclide in nuclides:
            xml_elements.append(self._get_nuclide_xml(nuclide, distrib))

        return xml_elements

    def _get_elements_xml(self, elements, cross_sections, distrib=False):

        xml_elements = []

        for element in elements:
            nuclide_elements = self._get_element_xml(element, cross_sections,
                                                     distrib)
            for nuclide_element in nuclide_elements:
                xml_elements.append(nuclide_element)

        return xml_elements

    def get_material_xml(self, cross_sections):
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

        # Create temperature XML subelement
        if self.temperature is not None:
            subelement = ET.SubElement(element, "temperature")
            subelement.text = str(self.temperature)

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
                subelements = self._get_elements_xml(self._elements,
                                                     cross_sections)
                for subelement in subelements:
                    element.append(subelement)
            else:
                # Create macroscopic XML subelements
                subelement = self._get_macroscopic_xml(self._macroscopic)
                element.append(subelement)

        else:
            subelement = ET.SubElement(element, "compositions")

            comps = []
            allnucs = self._nuclides + self._elements
            dist_per_type = allnucs[0][2]
            for nuc in allnucs:
                if nuc[2] != dist_per_type:
                    msg = 'All nuclides and elements in a distributed ' \
                          'material must have the same type, either ao or wo'
                    raise ValueError(msg)
                comps.append(nuc[1])

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
                subelements = self._get_nuclides_xml(self._nuclides,
                                                     distrib=True)
                for subelement_nuc in subelements:
                    subelement.append(subelement_nuc)

                # Create element XML subelements
                subelements = self._get_elements_xml(self._elements,
                                                     cross_sections,
                                                     distrib=True)
                for subsubelement in subelements:
                    subelement.append(subsubelement)
            else:
                # Create macroscopic XML subelements
                subsubelement = self._get_macroscopic_xml(self._macroscopic)
                subelement.append(subsubelement)

        if len(self._sab) > 0:
            for sab in self._sab:
                subelement = ET.SubElement(element, "sab")
                subelement.set("name", sab)

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
    cross_sections : str
        Indicates the path to an XML cross section listing file (usually named
        cross_sections.xml). If it is not set, the
        :envvar:`OPENMC_CROSS_SECTIONS` environment variable will be used for
        continuous-energy calculations and
        :envvar:`OPENMC_MG_CROSS_SECTIONS` will be used for multi-group
        calculations to find the path to the XML cross section file.
    multipole_library : str
        Indicates the path to a directory containing a windowed multipole
        cross section library. If it is not set, the
        :envvar:`OPENMC_MULTIPOLE_LIBRARY` environment variable will be used. A
        multipole library is optional.

    """

    def __init__(self, materials=None):
        super(Materials, self).__init__(Material, 'materials collection')

        self._materials_file = ET.Element("materials")
        self._cross_sections = None
        self._multipole_library = None

        if materials is not None:
            self += materials

    @property
    def cross_sections(self):
        return self._cross_sections

    @property
    def multipole_library(self):
        return self._multipole_library

    @cross_sections.setter
    def cross_sections(self, cross_sections):
        cv.check_type('cross sections', cross_sections, string_types)
        self._cross_sections = cross_sections

    @multipole_library.setter
    def multipole_library(self, multipole_library):
        cv.check_type('cross sections', multipole_library, string_types)
        self._multipole_library = multipole_library

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
        for material in self:
            xml_element = material.get_material_xml(self.cross_sections)
            self._materials_file.append(xml_element)

    def _create_cross_sections_subelement(self):
        if self._cross_sections is not None:
            element = ET.SubElement(self._materials_file, "cross_sections")
            element.text = str(self._cross_sections)

    def _create_multipole_library_subelement(self):
        if self._multipole_library is not None:
            element = ET.SubElement(self._materials_file, "multipole_library")
            element.text = str(self._multipole_library)

    def export_to_xml(self, path='materials.xml'):
        """Export material collection to an XML file.

        Parameters
        ----------
        path : str
            Path to file to write. Defaults to 'materials.xml'.

        """

        # Reset xml element tree
        self._materials_file.clear()

        self._create_material_subelements()
        self._create_cross_sections_subelement()
        self._create_multipole_library_subelement()

        # Clean the indentation in the file to be user-readable
        sort_xml_elements(self._materials_file)
        clean_xml_indentation(self._materials_file)

        # Write the XML Tree to the materials.xml file
        tree = ET.ElementTree(self._materials_file)
        tree.write(path, xml_declaration=True, encoding='utf-8', method="xml")
