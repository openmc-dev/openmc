#/usr/bin/env python

"""
Script for generating urr.xml input file
"""

import os
import lxml.etree as et
import urr_xml as xml
import parse_endf_6 as endf6

urr_files = endf6.urr_filenames('/g/g19/walsh23/data/evaluations/endf-7.1/neutrons/')
zaids, symbols = endf6.zaids_symbols(urr_files)
for i in range(len(urr_files)):
    urr = xml.URRElement()
    settings = xml.SettingsElement()
    isotopes = xml.IsotopesListElement()
    probability_tables = xml.ProbabilityTablesElement()

    # URR settings
    settings.endf_6_filepath(
        '/g/g19/walsh23/data/evaluations/endf-7.1/neutrons/')
    settings.avg_xs_filepath(
        '/g/g19/walsh23/dev/shutmc/src/URR/data/processed/PURXS/avg_xs/')
    settings.formalism('SLBW')
    settings.w_function_implementation('quick_w')
    settings.num_s_wave(64)
    settings.num_p_wave(64)
    settings.num_d_wave(64)
    settings.num_f_wave(64)
    settings.parameter_energy_dependence('neutron')
    settings.competitive_structure(True)

    # URR isotopes
    isotope = xml.IsotopeElement()
    isotope.zaid(zaids[i])
    isotope.endf_6_file(urr_files[i])
    isotopes.add_isotope(isotope)

    # URR probability tables input
    probability_tables.read_tables(False)
    probability_tables.filepath('/g/g19/walsh23/dev/shutmc/src/URR/data/processed/PURXS/prob_tables/')
    probability_tables.temperature_interpolation_method('log-log')
    probability_tables.num_bands(16)
    probability_tables.temperature_grid([0.0])
    probability_tables.min_num_batches(64)
    probability_tables.max_num_batches(1024)
    probability_tables.num_histories(128)
    probability_tables.tolerance(5.0e-3)
    probability_tables.energy_spacing('logarithmic')
    probability_tables.num_energies(16)
    probability_tables.write_tables(True)
    probability_tables.write_avg_xs(True)

    # Construct urr.xml and materials.xml
    os.makedirs(settings.avg_xs_filepath_.text+symbols[i])
    os.chdir(settings.avg_xs_filepath_.text+symbols[i])
    urr.add_element(settings)
    urr.add_element(isotopes)
    urr.add_element(probability_tables)
    urr.write_xml()
    xml.urr_materials_xml(symbols[i])
