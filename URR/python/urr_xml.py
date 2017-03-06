#!/usr/bin/env python

"""
API for generating OpenMC XML input files for the URR package

public classes:
URRElement               -- top-level object
SettingsElement          -- object containing required sub-element of <urr>
IsotopesElement          -- object containing required sub-element of <urr>
ProbabilityTablesElement -- object containing optional sub-element of <urr>
OnTheFlyElement          -- object containing optional sub-element of <urr>
PointwiseElement         -- object containing optional sub-element of <urr>
"""

import argparse
import lxml.etree as et
import copy as cp
import logging as log


class URRElement(object):
    """
    required top-level object for generating urr.xml containing <urr> element

    data:
    xml_element -- the actual <urr> XML element
    methods:
    __init__
    add_element
    """

    def __init__(self):
        self.xml_element_ = et.Element('urr')

    def add_element(self, element):
        """add an XML sub-element to the <urr> XML element"""

        self.xml_element_.append(cp.deepcopy(element.xml_element_))

    def write_xml(self):
        urr_xml = open('urr.xml', 'w')
        urr_xml.write(et.tostring(self.xml_element_, pretty_print=True))
        urr_xml.close()


class SettingsElement(object):
    """
    required object containing <settings> sub-element of <urr>
    """

    def __init__(self):
        self.xml_element_ = et.Element('settings')

    def endf_6_filepath(self, endf_6_filepath_str):
        self.endf_6_filepath_ = et.SubElement(
            self.xml_element_, 'endf_6_filepath')
        self.endf_6_filepath_.text = endf_6_filepath_str

    def avg_xs_filepath(self, avg_xs_filepath_str):
        self.avg_xs_filepath_ = et.SubElement(
            self.xml_element_, 'avg_xs_filepath')
        self.avg_xs_filepath_.text = avg_xs_filepath_str

    def formalism(self, formalism_str):
        self.formalism_ = et.SubElement(
            self.xml_element_, 'formalism')
        self.formalism_.text = formalism_str

    def w_function_implementation(self, w_function_implementation_str):
        self.w_function_implementation_ = et.SubElement(
            self.xml_element_, 'w_function_implementation')
        self.w_function_implementation_.text = w_function_implementation_str

    def num_s_wave(self, num_s_wave_int):
        self.num_s_wave_ = et.SubElement(
            self.xml_element_, 'num_s_wave')
        self.num_s_wave_.text = str(num_s_wave_int)

    def num_p_wave(self, num_p_wave_int):
        self.num_p_wave_ = et.SubElement(
            self.xml_element_, 'num_p_wave')
        self.num_p_wave_.text = str(num_p_wave_int)

    def num_d_wave(self, num_d_wave_int):
        self.num_d_wave_ = et.SubElement(
            self.xml_element_, 'num_d_wave')
        self.num_d_wave_.text = str(num_d_wave_int)

    def num_f_wave(self, num_f_wave_int):
        self.num_f_wave_ = et.SubElement(
            self.xml_element_, 'num_f_wave')
        self.num_f_wave_.text = str(num_f_wave_int)

    def parameter_energy_dependence(self, parameter_energy_dependence_str):
        self.parameter_energy_dependence_ = et.SubElement(
            self.xml_element_, 'parameter_energy_dependence')
        self.parameter_energy_dependence_.text \
            = parameter_energy_dependence_str

    def competitive_structure(self, competitive_structure_bool):
        self.competitive_structure_ = et.SubElement(
            self.xml_element_, 'competitive_structure')
        self.competitive_structure_.text = str(competitive_structure_bool)


class IsotopesListElement(object):
    """
    required object containing <isotopes> sub-element of <urr>
    """

    def __init__(self):
        self.xml_element_ = et.Element('isotopes')

    def add_isotope(self, isotope):
        """add an isotope sub-element to the <isotopes> XML element"""

        self.xml_element_.append(cp.deepcopy(isotope.xml_element_))


class IsotopeElement(object):
    """
    required object containing <isotope> sub-element of <isotopes>
    """

    def __init__(self):
        self.xml_element_ = et.Element('isotope')

    def zaid(self, zaid_int):
        self.xml_element_.set('zaid', str(zaid_int))

    def endf_6_file(self, endf_6_file_str):
        self.xml_element_.set('endf_6_file', endf_6_file_str)

    def max_energy(self, max_energy_float):
        self.xml_element_.set('max_energy', str(max_energy_float))


class OnTheFlyElement(object):
    """
    optional object containing <on_the_fly> sub-element of <urr>
    """

    def __init__(self):
        self.xml_element_ = et.Element('on_the_fly')

    def frequency(self, frequency_str):
        self.frequency_ = et.SubElement(
            self.xml_element_, 'frequency')
        self.frequency_.text = frequency_str

    def num_realizations(self, num_realizations_int):
        self.num_realizations_ = et.SubElement(
            self.xml_element_, 'num_realizations')
        self.num_realizations_.text = str(num_realizations_int)

    def realization_index(self, realization_index_int):
        self.realization_index_ = et.SubElement(
            self.xml_element_, 'realization_index')
        self.realization_index_.text = str(realization_index_int)


class ProbabilityTablesElement(object):
    """
    optional object containing <probability_tables> sub-element of <urr>
    """

    def __init__(self):
        self.xml_element_ = et.Element('probability_tables')

    def read_tables(self, read_tables_bool):
        self.read_tables_ = et.SubElement(
            self.xml_element_, 'read_tables')
        self.read_tables_.text = str(read_tables_bool)

    def filepath(self, filepath_str):
        self.filepath_ = et.SubElement(
            self.xml_element_, 'filepath')
        self.filepath_.text = filepath_str

    def temperature_interpolation_method(
            self, temperature_interpolation_method_str):
        self.temperature_interpolation_method_ = et.SubElement(
            self.xml_element_, 'temperature_interpolation_method')
        self.temperature_interpolation_method_.text \
            = temperature_interpolation_method_str

    def num_bands(self, num_bands_int):
        self.num_bands_ = et.SubElement(
            self.xml_element_, 'num_bands')
        self.num_bands_.text = str(num_bands_int)

    def temperature_grid(self, temperature_grid_list):
        temperature_grid = '\n      '
        if len(temperature_grid_list) == 1:
            temperature_grid = str(temperature_grid_list[0]) + '\n    '
        else:
            for i in range(len(temperature_grid_list)-1):
                temperature_grid += str(temperature_grid_list[i]) + '\n      '
            temperature_grid += str(temperature_grid_list[i+1]) + '\n    '
        self.temperature_grid_ = et.SubElement(
            self.xml_element_, 'temperature_grid')
        self.temperature_grid_.text = temperature_grid

    def min_num_batches(self, min_num_batches_int):
        self.min_num_batches_ = et.SubElement(
            self.xml_element_, 'min_num_batches')
        self.min_num_batches_.text = str(min_num_batches_int)

    def max_num_batches(self, max_num_batches_int):
        self.max_num_batches_ = et.SubElement(
            self.xml_element_, 'max_num_batches')
        self.max_num_batches_.text = str(max_num_batches_int)

    def num_histories(self, num_histories_int):
        self.num_histories_ = et.SubElement(
            self.xml_element_, 'num_histories')
        self.num_histories_.text = str(num_histories_int)

    def tolerance(self, tolerance_float):
        self.tolerance_ = et.SubElement(
            self.xml_element_, 'tolerance')
        self.tolerance_.text = str(tolerance_float)

    def energy_spacing(self, energy_spacing_str):
        self.energy_spacing_ = et.SubElement(
            self.xml_element_, 'energy_spacing')
        self.energy_spacing_.text = energy_spacing_str

    def num_energies(self, num_energies_int):
        self.num_energies_ = et.SubElement(
            self.xml_element_, 'num_energies')
        self.num_energies_.text = str(num_energies_int)

    def energy_grid(self, energy_grid_list):
        energy_grid = '\n      '
        if len(energy_grid_list) == 1:
            energy_grid = str(energy_grid_list[0]) + '\n    '
        else:
            for i in range(len(energy_grid_list)-1):
                energy_grid += str(energy_grid_list[i]) + '\n      '
            energy_grid += str(energy_grid_list[i+1]) + '\n    '
        self.energy_grid_ = et.SubElement(
            self.xml_element_, 'energy_grid')
        self.energy_grid_.text = energy_grid

    def background_xs(self, background_xs_str):
        self.background_xs_ = et.SubElement(
            self.xml_element_, 'backgound_xs')
        self.background_xs_.text = background_xs_str

    def write_tables(self, write_tables_bool):
        self.write_tables_ = et.SubElement(
            self.xml_element_, 'write_tables')
        self.write_tables_.text = str(write_tables_bool)

    def write_avg_xs(self, write_avg_xs_bool):
        self.write_avg_xs_ = et.SubElement(
            self.xml_element_, 'write_avg_xs')
        self.write_avg_xs_.text = str(write_avg_xs_bool)


class PointwiseElement(object):
    """
    optional object containing <pointwise> sub-element of <urr>
    """

    def __init__(self):
        self.xml_element_ = et.Element('pointwise')

    def min_energy_spacing(self, min_energy_spacing_float):
        self.min_energy_spacing_ = et.SubElement(
            self.xml_element_, 'min_energy_spacing')
        self.min_energy_spacing_.text = str(min_energy_spacing_float)

    def tolerance(self, tolerance_float):
        self.tolerance_ = et.SubElement(
            self.xml_element_, 'tolerance')
        self.tolerance_.text = str(tolerance_float)

def urr_materials_xml(symbols):
    """Write materials.xml XML elements for URR isotopes"""
    matfile = open('materials.xml', 'w')
    for i in range(len(symbols)):
        matfile.write('<nuclide name="'+symbols[i]+'" ao="1.0" />'+'\n')
    matfile.close()

def urr_get_args():
    """Return list of isotopes for URR processing"""
    parser = argparse.ArgumentParser()
    parser.add_argument('isotopes', nargs='+',
        help='space-separated list of isotopes for URR processing (e.g., U-235 U-238 Pu-239 Pu-240)')
    args = parser.parse_args()

    return args.isotopes
