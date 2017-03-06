#!/usr/bin/env python

"""
Script for generating URR data generation urr.xml input files
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
probability_tables.read_tables(False)
probability_tables.temperature_interpolation_method('log-log')
probability_tables.num_bands(1)
probability_tables.temperature_grid(['293.6'])
probability_tables.min_num_batches(50)
probability_tables.max_num_batches(1000)
probability_tables.num_histories(128)
probability_tables.tolerance(5.0e-3)
probability_tables.energy_spacing('logarithmic')
probability_tables.num_energies(64)
probability_tables.write_tables(True)
probability_tables.write_avg_xs(True)

# URR isotopes input
urr_files = endf6.urr_filenames(endf_6_filepath)
urr_files.sort()
zaids, symbols = endf6.zaids_symbols(urr_files, 'JENDL')
for i in range(len(urr_files)):
    include_isotope = False
    if ('all' in urr_isotopes) or (symbols[i] in urr_isotopes):
        include_isotope = True
    if include_isotope:
        isotopes = xml.IsotopesListElement()
        isotope = xml.IsotopeElement()
        isotope.zaid(zaids[i])
        isotope.endf_6_file(urr_files[i])
        isotopes.add_isotope(isotope)

        # Generate urr.xml
        os.makedirs(settings.avg_xs_filepath_.text+symbols[i])
        os.chdir(settings.avg_xs_filepath_.text+symbols[i])
        urr = xml.URRElement()
        urr.add_element(settings)
        urr.add_element(isotopes)
        urr.add_element(probability_tables)
        urr.write_xml()
