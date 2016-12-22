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
    '/Users/walshjon/Dropbox (Personal)/transfer/kilkenny/walshjon/data/evaluations/endf-7.1/neutrons')
settings.avg_xs_filepath(
    '/Users/walshjon/Dropbox (Personal)/transfer/kilkenny/walshjon/data/avg-urr-xs/endf-7.1/slbw')
settings.formalism('SLBW')
settings.w_function_implementation('quick_w')
settings.num_s_wave(64)
settings.num_p_wave(64)
settings.num_d_wave(64)
settings.num_f_wave(64)
settings.parameter_energy_dependence('neutron')
settings.competitive_structure(True)

urr_files = endf6.urr_filenames()
zaids, symbols = endf6.zaids_symbols(urr_files)
for i in range(len(urr_files)):
    isotope = xml.IsotopeElement()
    isotope.zaid(zaids[i])
    isotope.endf_6_file(urr_files[i])
    isotopes.add_isotope(isotope)
#    isotope.max_energy(20.0e6) # defaults to infinity

probability_tables.load_tables(False)
#probability_tables.filepath('/path/to/prob_table/data')
#probability_tables.temperature_interpolation_method('loglog')
probability_tables.num_bands(1)
probability_tables.temperature_grid([293.6])
probability_tables.min_num_batches(64)
probability_tables.max_num_batches(1024)
probability_tables.num_histories(64)
probability_tables.tolerance(1.0e-3)
probability_tables.energy_spacing('logarithmic')
probability_tables.num_energies(16)
#probability_tables.energy_grid()
#probability_tables.background_xs('ENDF_6') # defaults to no background component
probability_tables.generate_tables(False)
probability_tables.generate_avg_xs(True)

urr.add_element(settings)
urr.add_element(isotopes)
urr.add_element(probability_tables)

urr.write_xml()
