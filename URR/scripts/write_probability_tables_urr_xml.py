#!/usr/bin/env python

"""
Script for generating probability table URR cross sections urr.xml input file
"""

import os
import urr_xml as xml
import parse_endf_6 as endf6

urr_isotopes = xml.urr_get_args()

# URR settings
settings = xml.SettingsElement()
endf_6_filepath = os.environ.get('ENDF_6_FILEPATH')
if endf_6_filepath:
    if endf_6_filepath[-1] != '/': endf_6_filepath += '/'
    settings.endf_6_filepath(endf_6_filepath)
else:
    raise ValueError('Set ENDF_6_FILEPATH environment variable')
avg_xs_filepath = os.environ.get('AVG_XS_FILEPATH')
if avg_xs_filepath:
    if avg_xs_filepath[-1] != '/': avg_xs_filepath += '/'
    settings.avg_xs_filepath(avg_xs_filepath)
else:
    raise ValueError('Set AVG_XS_FILEPATH environment variable')
settings.formalism('SLBW')
settings.w_function_implementation('quick_w')
settings.num_s_wave(64)
settings.num_p_wave(64)
settings.num_d_wave(64)
settings.num_f_wave(64)
settings.parameter_energy_dependence('neutron')
settings.competitive_structure(True)

# URR probability tables input
probability_tables = xml.ProbabilityTablesElement()
prob_tables_filepath = os.environ.get('PROB_TABLES_FILEPATH')
if prob_tables_filepath:
    if prob_tables_filepath[-1] != '/': prob_tables_filepath += '/'
    probability_tables.filepath(prob_tables_filepath)
else:
    raise ValueError('Set PROB_TABLES_FILEPATH environment variable')
probability_tables.read_tables(True)
probability_tables.temperature_interpolation_method('log-log')

# URR isotopes input
isotopes = xml.IsotopesListElement()
urr_files = endf6.urr_filenames(endf_6_filepath)
urr_files.sort()
zaids, symbols = endf6.zaids_symbols(urr_files)
for i in range(len(urr_files)):
    include_isotope = False
    if ('all' in urr_isotopes) or (symbols[i] in urr_isotopes):
        include_isotope = True
    if include_isotope:
        isotope = xml.IsotopeElement()
        isotope.zaid(zaids[i])
        isotope.endf_6_file(urr_files[i])
        isotopes.add_isotope(isotope)

# Generate urr.xml
urr = xml.URRElement()
urr.add_element(settings)
urr.add_element(isotopes)
urr.add_element(probability_tables)
urr.write_xml()
