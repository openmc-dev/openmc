#/usr/bin/env python

"""
Script for generating urr.xml input file
"""

import urr_xml as xml
import parse_endf_6 as endf6

urr = xml.URRElement()
settings = xml.SettingsElement()
isotopes = xml.IsotopesListElement()
probability_tables = xml.ProbabilityTablesElement()

settings.endf_6_filepath(
    '/Users/walshjon/dev/shutmc/src/URR/data/evaluations/endf-7.1/neutrons/')
settings.avg_xs_filepath(
    '/Users/walshjon/dev/shutmc/src/URR/data/avg_xs/')
settings.formalism('SLBW')
settings.w_function_implementation('quick_w')
settings.num_s_wave(64)
settings.num_p_wave(64)
settings.num_d_wave(64)
settings.num_f_wave(64)
settings.parameter_energy_dependence('neutron')
settings.competitive_structure(True)

urr_files = endf6.urr_filenames('/Users/walshjon/dev/shutmc/src/URR/data/evaluations/endf-7.1/neutrons/')
zaids, symbols = endf6.zaids_symbols(urr_files)
for i in range(len(urr_files)):
    isotope = xml.IsotopeElement()
    isotope.zaid(zaids[i])
    isotope.endf_6_file(urr_files[i])
    isotopes.add_isotope(isotope)

probability_tables.read_tables(False)
probability_tables.filepath('/Users/walshjon/dev/shutmc/src/URR/data/probability_tables/')
probability_tables.temperature_interpolation_method('log-log')
probability_tables.num_bands(16)
probability_tables.temperature_grid([0.0,
                                     293.6,
                                     600.,
                                     900.,
                                     1200.,
                                     2500.,
                                     11604.5,
                                     116045.2,
                                     1160452.2,
                                     11604522.1,
                                     116045221.1,
                                     348135663.3,
                                     1160452211.1,
                                     11604522110.5])
probability_tables.min_num_batches(64)
probability_tables.max_num_batches(1024)
probability_tables.num_histories(64)
probability_tables.tolerance(1.0e-3)
probability_tables.energy_spacing('logarithmic')
probability_tables.num_energies(16)
probability_tables.write_tables(True)
probability_tables.write_avg_xs(True)

urr.add_element(settings)
urr.add_element(isotopes)
urr.add_element(probability_tables)

urr.write_xml()
xml.urr_materials_xml(symbols)
