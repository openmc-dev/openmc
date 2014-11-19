#!/usr/bin/env python

import warnings

import openmc
from openmc.checkvalue import *
from openmc.clean_xml import *
from xml.etree import ElementTree as ET
from collections import MappingView
from copy import deepcopy
import numpy as np


# A list of all IDs for all Materials created
MATERIAL_IDS = list()

# A static variable for auto-generated Material IDs
AUTO_MATERIAL_ID = 10000

def reset_auto_material_id():
    global AUTO_MATERIAL_ID, MATERIAL_IDS
    AUTO_MATERIAL_ID = 10000
    MATERIAL_IDS = list()


# Units for density supported by OpenMC
DENSITY_UNITS = ['g/cm3', 'g/cc', 'kg/cm3', 'at/b-cm', 'at/cm3', 'sum']

# ENDF temperatures
ENDF_TEMPS = np.array([300, 600, 700, 900, 1200, 1500])

# ENDF ZAIDs
ENDF_ZAIDS = np.array(['70c', '71c', '72c', '73c', '74c'])

# Constant for density when not needed
NO_DENSITY = 99999.


class Material(object):

    def __init__(self, material_id=None, name=''):

        # Initialize class attributes
        self._id = None
        self._name = ''
        self._density = None
        self._density_units = ''

        # A dictionary of Nuclides
        # Keys         - Nuclide names
        # Values     - tuple (nuclide, percent, percent type)
        self._nuclides = dict()

        # A dictionary of Elements
        # Keys         - Element names
        # Values     - tuple (element, percent, percent type)
        self._elements = dict()

        # If specified, a list of tuples of (table name, xs identifier)
        self._sab = list()

        # If true, the material will be initialized as distributed
        self._convert_to_distrib_comps = False

        # If specified, this file will be used instead of composition values
        self._distrib_otf_file = None

        # Set the Material class attributes
        self.set_id(material_id)
        self.set_name(name)


    def set_id(self, material_id=None):

        global MATERIAL_IDS

        # If the Material already has an ID, remove it from global list
        if not self._id is None:
            MATERIAL_IDS.remove(self._id)

        if material_id is None:
            global AUTO_MATERIAL_ID
            self._id = AUTO_MATERIAL_ID
            MATERIAL_IDS.append(AUTO_MATERIAL_ID)
            AUTO_MATERIAL_ID += 1

        # Check that the ID is an integer and wasn't already used
        elif not is_integer(material_id):
            msg = 'Unable to set a non-integer Material ' \
                  'ID {0}'.format(material_id)
            raise ValueError(msg)

        elif material_id in MATERIAL_IDS:
            msg = 'Unable to set Material ID to {0} since a Material with ' \
                  'this ID was already initialized'.format(material_id)
            raise ValueError(msg)

        elif material_id < 0:
            msg = 'Unable to set Material ID to {0} since it must be a ' \
                  'non-negative integer'.format(material_id)
            raise ValueError(msg)

        else:
            self._id = material_id
            MATERIAL_IDS.append(material_id)


    def set_name(self, name):

        if not is_string(name):
            msg = 'Unable to set name for Material ID={0} with a non-string ' \
                  'value {1}'.format(self._id, name)
            raise ValueError(msg)

        else:
            self._name = name


    def set_density(self, units, density=NO_DENSITY):

        if not is_float(density):
            msg = 'Unable to set the density for Material ID={0} to a ' \
                  'non-floating point value {1}'.format(self._id, density)
            raise ValueError(msg)

        elif not units in DENSITY_UNITS:
            msg = 'Unable to set the density for Material ID={0} with ' \
                  'units {1}'.format(self._id, units)
            raise ValueError(msg)

        if density == NO_DENSITY and units is not 'sum':
            msg = 'Unable to set the density Material ID={0} ' \
                  'because a density must be set when not using ' \
                  'sum unit'.format(self._id)
            raise ValueError(msg)

        self._density = density
        self._density_units = units


    def add_nuclide(self, nuclide, percent, percent_type='ao'):

        if not isinstance(nuclide, openmc.Nuclide):
            msg = 'Unable to add an Nuclide to Material ID={0} with a ' \
                  'non-Nuclide value {1}'.format(self._id, nuclide)
            raise ValueError(msg)

        elif not is_float(percent):
            msg = 'Unable to add an Nuclide to Material ID={0} with a ' \
                  'non-floating point value {1}'.format(self._id, percent)
            raise ValueError(msg)

        elif not percent_type in ['ao', 'wo', 'at/g-cm']:
            msg = 'Unable to add an Nuclide to Material ID={0} with a ' \
                  'percent type {1}'.format(self._id, percent_type)
            raise ValueError(msg)

        # Copy this Nuclide to separate it from the Nuclide in other Materials
        nuclide = deepcopy(nuclide)

        self._nuclides[nuclide._name] = (nuclide, percent, percent_type)


    def remove_nuclide(self, nuclide):

        if not isinstance(nuclide, openmc.Nuclide):
            msg = 'Unable to remove a Nuclide {0} in Material ID={1} ' \
                  'since it is not a Nuclide'.format(self._id, nuclide)
            raise ValueError(msg)

        # If the Material contains the Nuclide, delete it
        if nuclide._name in self._nuclides:
            del self._nuclides[nuclide._name]


    def add_element(self, element, percent, percent_type='ao'):

        if not isinstance(element, openmc.Element):
            msg = 'Unable to add an Element to Material ID={0} with a ' \
                  'non-Element value {1}'.format(self._id, element)
            raise ValueError(msg)

        if not is_float(percent):
            msg = 'Unable to add an Element to Material ID={0} with a ' \
                  'non-floating point value {1}'.format(self._id, percent)
            raise ValueError(msg)

        if not percent_type in ['ao', 'wo']:
            msg = 'Unable to add an Element to Material ID={0} with a ' \
                  'percent type {1}'.format(self._id, percent_type)
            raise ValueError(msg)

        # Copy this Element to separate it from same Element in other Materials
        element = deepcopy(element)

        self._elements[element._name] = (element, percent, percent_type)


    def remove_element(self, element):

        # If the Material contains the Element, delete it
        if element._name in self._elements:
            del self._elements[element._name]


    def add_s_alpha_beta(self, name, xs):

        if not is_string(name):
            msg = 'Unable to add an S(a,b) table to Material ID={0} with a ' \
                        'non-string table name {1}'.format(self._id, name)
            raise ValueError(msg)

        if not is_string(xs):
            msg = 'Unable to add an S(a,b) table to Material ID={0} with a ' \
                  'non-string cross-section identifier {1}'.format(self._id, xs)
            raise ValueError(msg)

        self._sab.append((name, xs))


    def set_otf_mat_file(self, name):

        # TODO: remove this when distributed materials are merged
        warnings.warn('This feature is not yet implemented in a release ' \
                      'version of openmc')

        if not is_string(name):
            msg = 'Unable to add OTF material file to Material ID={0} with a ' \
                        'non-string name {1}'.format(self._id, name)
            raise ValueError(msg)

        self._distrib_otf_file = name


    def set_as_distrib_comp(self):

        # TODO: remove this when distributed materials are merged
        warnings.warn('This feature is not yet implemented in a release ' \
                      'version of openmc')

        self._convert_to_distrib_comps = True


    def get_all_nuclides(self):

        nuclides = dict()

        for nuclide_name, nuclide_tuple in self._nuclides.items():
            nuclide = nuclide_tuple[0]
            density = nuclide_tuple[1]
            nuclides[nuclide._name] = (nuclide, density)

        return nuclides


    def _repr__(self):

        string = 'Material\n'
        string += '{0: <16}{1}{2}\n'.format('\tID', '=\t', self._id)
        string += '{0: <16}{1}{2}\n'.format('\tName', '=\t', self._name)

        string += '{0: <16}{1}{2}'.format('\tDensity', '=\t', self._density)
        string += ' [{0}]\n'.format(self._density_units)

        string += '{0: <16}'.format('\tS(a,b) Tables') + '\n'

        for sab in self._sab:
            string += '{0: <16}{1}[{2}{3}]\n'.format('\tS(a,b)', '=\t',
                                                     sab[0], sab[1])

        string += '{0: <16}'.format('\tNuclides') + '\n'

        for nuclide in self._nuclides:
            percent = self._nuclides[nuclide][1]
            percent_type = self._nuclides[nuclide][2]
            string += '{0: <16}'.format('\t{0}'.format(nuclide))
            string += '=\t{0: <12} [{1}]\n'.format(percent, percent_type)

        string += '{0: <16}\n'.format('\tElements')

        for element in self._elements:
            percent = self._nuclides[element][1]
            percent_type = self._nuclides[element][2]
            string += '{0: >16}'.format('\t{0}'.format(element))
            string += '=\t{0: <12} [{1}]\n'.format(percent, percent_type)

        return string


    def get_nuclide_xml(self, nuclide, distrib=False):

        xml_element = ET.Element("nuclide")
        xml_element.set("name", nuclide[0]._name)

        if not distrib:
            if nuclide[2] is 'ao':
                xml_element.set("ao", str(nuclide[1]))
            else:
                xml_element.set("wo", str(nuclide[1]))

        if not nuclide[0]._xs is None:
            xml_element.set("xs", nuclide[0]._xs)

        return xml_element


    def get_element_xml(self, element, distrib=False):

        xml_element = ET.Element("element")
        xml_element.set("name", str(element[0]._name))

        if not distrib:
            if element[2] is 'ao':
                xml_element.set("ao", str(element[1]))
            else:
                xml_element.set("wo", str(element[1]))

        return xml_element


    def get_nuclides_xml(self, nuclides, distrib=False):

        xml_elements = list()

        for nuclide in nuclides.values():
            xml_elements.append(self.get_nuclide_xml(nuclide, distrib))

        return xml_elements


    def get_elements_xml(self, elements, distrib=False):

        xml_elements = list()

        for element in elements.values():
            xml_elements.append(self.get_element_xml(element, distrib))

        return xml_elements


    def get_material_xml(self):

        # Create Material XML element
        element = ET.Element("material")
        element.set("id", str(self._id))

        # Create density XML subelement
        subelement = ET.SubElement(element, "density")
        if self._density_units is not 'sum':
            subelement.set("value", str(self._density))
        subelement.set("units", self._density_units)

        if not self._convert_to_distrib_comps:

            # Create nuclide XML subelements
            subelements = self.get_nuclides_xml(self._nuclides)
            for subelement in subelements:
                element.append(subelement)

            # Create element XML subelements
            subelements = self.get_elements_xml(self._elements)
            for subelement in subelements:
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


            # Create nuclide XML subelements
            subelements = self.get_nuclides_xml(self._nuclides, distrib=True)
            for subelement_nuc in subelements:
                subelement.append(subelement_nuc)

            # Create element XML subelements
            subelements = self.get_elements_xml(self._elements, distrib=True)
            for subelement_ele in subelements:
                subelement.append(subelement_ele)

        if len(self._sab) > 0:
            for sab in self._sab:
                subelement = ET.SubElement(element, "sab")
                subelement.set("name", sab[0])
                subelement.set("xs", sab[1])

        return element


class MaterialsFile(object):

    def __init__(self):

        # Initialize MaterialsFile class attributes
        self._materials = list()
        self._default_xs = None
        self._materials_file = ET.Element("materials")


    def add_material(self, material):

        if not isinstance(material, Material):
            msg = 'Unable to add a non-Material {0} to the ' \
                  'MaterialsFile'.format(material)
            raise ValueError(msg)

        self._materials.append(material)


    def add_materials(self, materials):

        if not isinstance(materials, (tuple, list, MappingView)):
            msg = 'Unable to create OpenMC materials.xml file from {0} which ' \
                  'is not a Python tuple/list'.format(materials)
            raise ValueError(msg)

        for material in materials:
            self.add_material(material)


    def remove_materials(self, material):

        if not isinstance(material, Material):
            msg = 'Unable to remove a non-Material {0} from the ' \
                  'MaterialsFile'.format(material)
            raise ValueError(msg)

        self._materials.remove(material)


    def set_default_xs(self, xs):

        if not is_string(xs):
            msg = 'Unable to set default xs to a non-string value'.format(xs)
            raise ValueError(msg)

        self._default_xs = xs


    def create_material_subelements(self):

        subelement = ET.SubElement(self._materials_file, "default_xs")

        if not self._default_xs is None:
            subelement.text = self._default_xs

        for material in self._materials:
            xml_element = material.get_material_xml()

            if len(material._name) > 0:
                self._materials_file.append(ET.Comment(material._name))

            self._materials_file.append(xml_element)


    def export_to_xml(self):

        self.create_material_subelements()

        # Clean the indentation in the file to be user-readable
        sort_xml_elements(self._materials_file)
        clean_xml_indentation(self._materials_file)

        # Write the XML Tree to the materials.xml file
        tree = ET.ElementTree(self._materials_file)
        tree.write("materials.xml", xml_declaration=True,
                             encoding='utf-8', method="xml")
