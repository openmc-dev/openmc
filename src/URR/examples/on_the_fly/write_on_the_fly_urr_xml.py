#/usr/bin/env python

"""
Script for generating urr.xml input file
"""

import urr_xml as xml
import parse_endf_6 as endf6

urr = xml.URRElement()
settings = xml.SettingsElement()
isotopes = xml.IsotopesListElement()
on_the_fly = xml.OnTheFlyElement()

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

on_the_fly.frequency('event')

urr.add_element(settings)
urr.add_element(isotopes)
urr.add_element(on_the_fly)

urr.write_xml()
xml.urr_materials_xml(symbols)
