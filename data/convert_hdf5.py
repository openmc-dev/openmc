#!/usr/bin/env python

import glob
import os
from xml.dom.minidom import getDOMImplementation

import openmc.data.ace


if not os.path.isdir('nndc_hdf5'):
    os.mkdir('nndc_hdf5')

nndc_files = glob.glob('nndc/293.6K/*.ace')
nndc_thermal_files = glob.glob('nndc/tsl/*.acer')

thermal_names = {'al': 'c_Al27',
                 'be': 'c_Be',
                 'bebeo': 'c_Be_in_BeO',
                 'benzine': 'c_Benzine',
                 'dd2o': 'c_D_in_D2O',
                 'fe': 'c_Fe56',
                 'graphite': 'c_Graphite',
                 'hch2': 'c_H_in_CH2',
                 'hh2o': 'c_H_in_H2O',
                 'hzrh': 'c_H_in_ZrH',
                 'lch4': 'c_liquid_CH4',
                 'obeo': 'c_O_in_BeO',
                 'orthod': 'c_ortho_D',
                 'orthoh': 'c_ortho_H',
                 'ouo2': 'c_O_in_UO2',
                 'parad': 'c_para_D',
                 'parah': 'c_para_H',
                 'sch4': 'c_solid_CH4',
                 'uuo2': 'c_U_in_UO2',
                 'zrzrh': 'c_Zr_in_ZrH'}

impl = getDOMImplementation()
doc = impl.createDocument(None, "cross_sections", None)
doc_root = doc.documentElement

for f in sorted(nndc_files):
    print('Converting {}...'.format(f))

    # Deterine output file name
    dirname, basename = os.path.split(f)
    root, ext = os.path.splitext(basename)
    outfile = os.path.join('nndc_hdf5', root + '.h5')
    if os.path.exists(outfile):
        os.remove(outfile)

    # Determine elemental symbol, mass number and metastable state
    element, mass_number, temp = basename.split('_')
    metastable = int(mass_number[-1]) if 'm' in mass_number else 0
    mass_number = int(mass_number[:3])

    # Parse ACE file, create HDF5 file
    t = openmc.data.ace.get_table(f)
    t.export_to_hdf5(outfile, element, mass_number, metastable)
    xs = t.name.split('.')[1]
    if metastable > 0:
        name = "{}{}_m{}.{}".format(element, mass_number, metastable, xs)
    else:
        name = "{}{}.{}".format(element, mass_number, xs)

    # Add entry to XML listing
    libraryNode = doc.createElement("library")
    libraryNode.setAttribute("path", root + '.h5')
    libraryNode.setAttribute("materials", name)
    libraryNode.setAttribute("type", "neutron")
    doc_root.appendChild(libraryNode)

for f in sorted(nndc_thermal_files):
    print('Converting {}...'.format(f))

    # Deterine output file name
    dirname, basename = os.path.split(f)
    root, ext = os.path.splitext(basename)
    outfile = os.path.join('nndc_hdf5', root + '.h5')
    if os.path.exists(outfile):
        os.remove(outfile)

    # Parse ACE file, create HDF5 file
    t = openmc.data.ace.get_table(f)
    t.export_to_hdf5(outfile, thermal_names[root])
    xs = t.name.split('.')[1]

    # Add entry to XML listing
    libraryNode = doc.createElement("library")
    libraryNode.setAttribute("path", root + '.h5')
    libraryNode.setAttribute("materials", thermal_names[root] + '.' + xs)
    libraryNode.setAttribute("type", "thermal")
    doc_root.appendChild(libraryNode)

# Write cross_sections.xml
lines = doc.toprettyxml(indent='  ')
open(os.path.join('nndc_hdf5', 'cross_sections.xml'), 'w').write(lines)
