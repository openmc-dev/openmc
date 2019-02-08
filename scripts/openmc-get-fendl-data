#!/usr/bin/env python3

import os
from collections import defaultdict
import sys
import tarfile
import zipfile
import glob
import argparse
from string import digits
from urllib.request import urlopen, Request

import openmc.data


description = """
Download FENDL 3.1d or FENDL 3.1c ACE data from the IAEA and convert it to a HDF5 library for 
use with OpenMC.

"""


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

parser = argparse.ArgumentParser(
    description=description,
    formatter_class=CustomFormatter
)
parser.add_argument('-b', '--batch', action='store_true',
                    help='supresses standard in')
parser.add_argument('-d', '--destination', default=None,
                    help='Directory to create new library in')
parser.add_argument('--libver', choices=['earliest', 'latest'],
                    default='latest', help="Output HDF5 versioning. Use "
                    "'earliest' for backwards compatibility or 'latest' for "
                    "performance")
parser.add_argument('-r', '--release', choices=['3.1a', '3.1d'],
                    default='3.1d', help="The nuclear data library release version. "
                    "The currently supported options are 3.1a and 3.1d")
args = parser.parse_args()

# this could be added as an argument to allow different libraries to be downloaded
library = 'fendl'
ace_files_dir = '-'.join([library, args.release, 'ace'])
# the destination is decided after the release is know to avoid putting the release in a folder with a misleading name
if args.destination == None:
    args.destination = '-'.join([library, args.release, 'hdf5'])

# This dictionary contains all the unique information about each release. This can be exstened to accommodated new releases
release_details = {'3.1a': {'base_url': 'https://www-nds.iaea.org/fendl31/data/neutron/',
                            'files': ['fendl31a-neutron-ace.zip'],
                            'neutron_files': os.path.join(ace_files_dir, '*'),
                            'compressed_file_size': '0.4 GB',
                            'uncompressed_file_size': '3 GB'
                            },
                   '3.1d': {'base_url': 'https://www-nds.iaea.org/fendl/data/neutron/',
                            'files': ['fendl31d-neutron-ace.zip'],
                            'neutron_files': os.path.join(ace_files_dir, 'fendl31d_ACE', '*'),
                            'compressed_file_size': '0.5 GB',
                            'uncompressed_file_size': '3 GB'
                            }
                   }

download_warning = """
WARNING: This script will download {} of data. 
Extracting and processing the data requires {} of additional free disk space.

Are you sure you want to continue? ([y]/n)
""".format(release_details[args.release]['compressed_file_size'],
           release_details[args.release]['uncompressed_file_size'])

response = input(download_warning) if not args.batch else 'y'
if response.lower().startswith('n'):
    sys.exit()

block_size = 16384

# ==============================================================================
# DOWNLOAD FILES FROM IAEA SITE

# The fendl website requires the web browser to be mocked
user_agent = 'Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.9.0.7) Gecko/2009021910 Firefox/3.0.7'
headers = {'User-Agent': user_agent, }

files_complete = []
for f in release_details[args.release]['files']:
    # Establish connection to URL
    url = release_details[args.release]['base_url'] + f
    req = urlopen(Request(url, None, headers))

    # Get file size from header
    if sys.version_info[0] < 3:
        file_size = int(req.info().getheaders('Content-Length')[0])
    else:
        file_size = req.length
    downloaded = 0

    # Check if file already downloaded
    if os.path.exists(f):
        if os.path.getsize(f) == file_size:
            print('Skipping {}, already downloaded'.format(f))
            files_complete.append(f)
            continue
        else:
            overwrite = input('Overwrite {}? ([y]/n) '.format(f))
            if overwrite.lower().startswith('n'):
                continue

    # Copy file to disk
    print('Downloading {}... '.format(f), end='')
    with open(f, 'wb') as fh:
        while True:
            chunk = req.read(block_size)
            if not chunk: break
            fh.write(chunk)
            downloaded += len(chunk)
            status = '{:10}  [{:3.2f}%]'.format(downloaded, downloaded * 100. / file_size)
            print(status + chr(8)*len(status), end='')
        print('')
        files_complete.append(f)

# ==============================================================================
# EXTRACT FILES FROM TGZ

for f in release_details[args.release]['files']:
    if f not in files_complete:
        continue

    # Extract files, the fendl release was compressed using type 9 zip format
    # unfortunatly which is incompatible with the standard python zipfile library
    # therefore the following system command is used

    os.system('unzip -o ' + f + ' -d '+ ace_files_dir) 


# # ==============================================================================
# # GENERATE HDF5 LIBRARY -- NEUTRON FILES

# Get a list of all ACE files, excluding files ending with _ which are old incorrect files kept in the release for backwards compatability
neutron_files = [f for f in glob.glob(release_details[args.release]['neutron_files']) if not f.endswith('_') and not f.endswith('.xsd')]

# Create output directory if it doesn't exist
if not os.path.isdir(args.destination):
    os.mkdir(args.destination)

library = openmc.data.DataLibrary()

for filename in sorted(neutron_files):
    
    print('Converting: ' + filename)
    data = openmc.data.IncidentNeutron.from_ace(filename)

    # Export HDF5 file
    h5_file = os.path.join(args.destination, data.name + '.h5')
    print('Writing {}...'.format(h5_file))
    data.export_to_hdf5(h5_file, 'w', libver=args.libver)

    # Register with library
    library.register_file(h5_file)

# Write cross_sections.xml
libpath = os.path.join(args.destination, 'cross_sections.xml')
library.export_to_xml(libpath)
