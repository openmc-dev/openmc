#!/usr/bin/env python3

import argparse
from collections import defaultdict
import glob
import os

import openmc.data


description = """
Convert ENDF/B-VIII.0 ACE data from LANL into an HDF5 library
that can be used by OpenMC. This assumes that you have a directory containing
subdirectories 'Lib80x' and 'ENDF80SaB'.

"""


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass


parser = argparse.ArgumentParser(
    description=description,
    formatter_class=CustomFormatter
)
parser.add_argument('-o', '--output_dir', default='lib80x_hdf5',
                    help='Directory to create new library in')
parser.add_argument('--libver', choices=['earliest', 'latest'],
                    default='earliest', help="Output HDF5 versioning. Use "
                    "'earliest' for backwards compatibility or 'latest' for "
                    "performance")
parser.add_argument('--datadir', help='Directory containing Lib80x and ENDF80SaB',
                    default=os.curdir)
args = parser.parse_args()
assert os.path.isdir(args.datadir)

# Get a list of all ACE files
lib80x = glob.glob(os.path.join(args.datadir, 'Lib80x', '**', '*.80?nc'), recursive=True)
lib80sab = glob.glob(os.path.join(args.datadir, 'ENDF80SaB', '**', '*.??t'), recursive=True)

# Find and fix B10 ACE files
b10files = glob.glob(os.path.join(args.datadir, 'Lib80x', '**', '5010.80?nc'), recursive=True)
nxs1_position = 523
for filename in b10files:
    with open(filename, 'r+') as fh:
        # Read NXS(1)
        fh.seek(nxs1_position)
        nxs1 = int(fh.read(5))

        # Increase length to match actual length of XSS, but make sure this
        # isn't done twice by checking the current length
        if nxs1 < 86870:
            fh.seek(nxs1_position)
            fh.write(str(nxs1 + 53))

# Group together tables for the same nuclide
suffixes = defaultdict(list)
for filename in sorted(lib80x + lib80sab):
    dirname, basename = os.path.split(filename)
    zaid, xs = basename.split('.')
    suffixes[os.path.join(dirname, zaid)].append(xs)

# Create output directory if it doesn't exist
if not os.path.isdir(args.output_dir):
    os.mkdir(args.output_dir)

library = openmc.data.DataLibrary()

for basename, xs_list in sorted(suffixes.items()):
    # Convert first temperature for the table
    filename = '.'.join((basename, xs_list[0]))
    print('Converting: ' + filename)
    if filename.endswith('t'):
        data = openmc.data.ThermalScattering.from_ace(filename)
    else:
        data = openmc.data.IncidentNeutron.from_ace(filename, 'mcnp')

    # For each higher temperature, add cross sections to the existing table
    for xs in xs_list[1:]:
        filename = '.'.join((basename, xs))
        print('Adding: ' + filename)
        if filename.endswith('t'):
            data.add_temperature_from_ace(filename)
        else:
            data.add_temperature_from_ace(filename, 'mcnp')

    # Export HDF5 file
    h5_file = os.path.join(args.output_dir, data.name + '.h5')
    print('Writing {}...'.format(h5_file))
    data.export_to_hdf5(h5_file, 'w', libver=args.libver)

    # Register with library
    library.register_file(h5_file)

# Write cross_sections.xml
libpath = os.path.join(args.output_dir, 'cross_sections.xml')
library.export_to_xml(libpath)
