#!/usr/bin/env python3

import argparse
from collections import defaultdict
import glob
import os

import openmc.data


description = """
Convert ENDF/B-VII.0 ACE data from the MCNP5/6 distribution into an HDF5 library
that can be used by OpenMC. This assumes that you have a directory containing
files named endf70a, endf70b, ..., endf70k, and endf70sab.

"""

class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

parser = argparse.ArgumentParser(
    description=description,
    formatter_class=CustomFormatter
)
parser.add_argument('-d', '--destination', default='mcnp_endfb70',
                    help='Directory to create new library in')
parser.add_argument('--libver', choices=['earliest', 'latest'],
                    default='earliest', help="Output HDF5 versioning. Use "
                    "'earliest' for backwards compatibility or 'latest' for "
                    "performance")
parser.add_argument('mcnpdata', help='Directory containing endf70[a-k] and endf70sab')
args = parser.parse_args()
assert os.path.isdir(args.mcnpdata)

# Get a list of all neutron ACE files
endf70 = glob.glob(os.path.join(args.mcnpdata, 'endf70[a-k]'))

# Create output directory if it doesn't exist
if not os.path.isdir(args.destination):
    os.mkdir(args.destination)

library = openmc.data.DataLibrary()

for path in sorted(endf70):
    print('Loading data from {}...'.format(path))
    lib = openmc.data.ace.Library(path)

    # Group together tables for the same nuclide
    tables = defaultdict(list)
    for table in lib.tables:
        zaid, xs = table.name.split('.')
        tables[zaid].append(table)

    for zaid, tables in sorted(tables.items()):
        # Convert first temperature for the table
        print('Converting: ' + tables[0].name)
        data = openmc.data.IncidentNeutron.from_ace(tables[0], 'mcnp')

        # For each higher temperature, add cross sections to the existing table
        for table in tables[1:]:
            print('Adding: ' + table.name)
            data.add_temperature_from_ace(table, 'mcnp')

        # Export HDF5 file
        h5_file = os.path.join(args.destination, data.name + '.h5')
        print('Writing {}...'.format(h5_file))
        data.export_to_hdf5(h5_file, 'w', libver=args.libver)

        # Register with library
        library.register_file(h5_file)

# Handle S(a,b) tables
endf70sab = os.path.join(args.mcnpdata, 'endf70sab')
if os.path.exists(endf70sab):
    lib = openmc.data.ace.Library(endf70sab)

    # Group together tables for the same nuclide
    tables = defaultdict(list)
    for table in lib.tables:
        name, xs = table.name.split('.')
        tables[name].append(table)

    for zaid, tables in sorted(tables.items()):
        # Convert first temperature for the table
        print('Converting: ' + tables[0].name)
        data = openmc.data.ThermalScattering.from_ace(tables[0])

        # For each higher temperature, add cross sections to the existing table
        for table in tables[1:]:
            print('Adding: ' + table.name)
            data.add_temperature_from_ace(table)

        # Export HDF5 file
        h5_file = os.path.join(args.destination, data.name + '.h5')
        print('Writing {}...'.format(h5_file))
        data.export_to_hdf5(h5_file, 'w', libver=args.libver)

        # Register with library
        library.register_file(h5_file)

# Write cross_sections.xml
libpath = os.path.join(args.destination, 'cross_sections.xml')
library.export_to_xml(libpath)
