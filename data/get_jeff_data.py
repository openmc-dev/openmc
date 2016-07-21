#!/usr/bin/env python

from __future__ import print_function
import os
import shutil
import subprocess
import sys
import tarfile
import zipfile
import glob
import hashlib
import argparse

import openmc.data

try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen


thermal_suffix = {20: '01t', 100: '02t', 293: '03t', 296: '03t', 323: '04t',
                  350: '05t', 373: '06t', 400: '07t', 423: '08t', 473: '09t',
                  500: '10t', 523: '11t', 573: '12t', 600: '13t', 623: '14t',
                  643: '15t', 647: '15t', 700: '16t', 773: '17t', 800: '18t',
                  1000: '19t', 1200: '20t', 1600: '21t', 2000: '22t',
                  3000: '23t'}



parser = argparse.ArgumentParser()
parser.add_argument('-b', '--batch', action='store_true',
                    help='supresses standard in')
args = parser.parse_args()

base_url = 'https://www.oecd-nea.org/dbforms/data/eva/evatapes/jeff_32/Processed/'
files = ['JEFF32-ACE-293K.tar.gz',
         'JEFF32-ACE-400K.tar.gz',
         'JEFF32-ACE-500K.tar.gz',
         'JEFF32-ACE-600K.tar.gz',
         'JEFF32-ACE-700K.tar.gz',
         'JEFF32-ACE-800K.zip',
         'JEFF32-ACE-900K.tar.gz',
         'JEFF32-ACE-1000K.tar.gz',
         'JEFF32-ACE-1200K.tar.gz',
         'JEFF32-ACE-1500K.tar.gz',
         'JEFF32-ACE-1800K.tar.gz',
         'TSLs.tar.gz']

block_size = 16384

# ==============================================================================
# DOWNLOAD FILES FROM OECD SITE

files_complete = []
for f in files:
    # Establish connection to URL
    url = base_url + f
    req = urlopen(url)

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
            if sys.version_info[0] < 3:
                overwrite = raw_input('Overwrite {}? ([y]/n) '.format(f))
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

for f in files:
    if f not in files_complete:
        continue

    # Extract files
    if f.endswith('.zip'):
        with zipfile.ZipFile(f, 'r') as zipf:
            print('Extracting {}...'.format(f))
            zipf.extractall('jeff-3.2')

    else:
        suffix = 'ACEs_293K' if '293' in f else ''
        with tarfile.open(f, 'r') as tgz:
            print('Extracting {}...'.format(f))
            tgz.extractall(os.path.join('jeff-3.2', suffix))

        # Remove thermal scattering tables from 293K data since they are
        # redundant
        if '293' in f:
            for path in glob.glob(os.path.join('jeff-3.2', 'ACEs_293K', '*-293.ACE')):
                os.remove(path)

# ==============================================================================
# FIX ERRORS

# A few nuclides at 400K has 03c instead of 04c
print('Assigning new cross section identifiers...')
wrong_nuclides = ['Mn55', 'Mo95', 'Nb93', 'Pd105', 'Pu239', 'Pu240', 'U235',
                  'U238', 'Y89']
for nuc in wrong_nuclides:
    path = os.path.join('jeff-3.2', 'ACEs_400K', nuc + '.ACE')
    print('    Fixing {} (03c --> 04c)...'.format(path))
    if os.path.isfile(path):
        text = open(path, 'r').read()
        text = text[:7] + '04c' + text[10:]
        open(path, 'w').write(text)

# ==============================================================================
# CHANGE ZAID FOR METASTABLES

metastables = glob.glob(os.path.join('jeff-3.2', '**', '*M.ACE'))
for path in metastables:
    print('    Fixing {} (ensure metastable)...'.format(path))
    text = open(path, 'r').read()
    mass_first_digit = int(text[3])
    if mass_first_digit <= 2:
        text = text[:3] + str(mass_first_digit + 4) + text[4:]
        open(path, 'w').write(text)

# ==============================================================================
# CHANGE IDENTIFIER FOR S(A,B) TABLES

thermals = glob.glob(os.path.join('jeff-3.2', 'ANNEX_6_3_STLs', '**', '*.ace'))
for path in thermals:
    print('    Fixing {} (unique suffix)...'.format(path))
    basename = os.path.basename(path)
    temperature = int(basename.split('-')[1][:-4])
    text = open(path, 'r').read()
    text = text[:7] + thermal_suffix[temperature] + text[10:]
    open(path, 'w').write(text)

# ==============================================================================
# CONVERT TO BINARY TO SAVE DISK SPACE

# get a list of all ACE files
ace_files = (glob.glob(os.path.join('jeff-3.2', '**', '*.ACE')) +
             glob.glob(os.path.join('jeff-3.2', 'ANNEX_6_3_STLs', '**', '*.ace')))

# Ask user to convert
if not args.batch:
    if sys.version_info[0] < 3:
        response = raw_input('Convert ACE files to binary? ([y]/n) ')
    else:
        response = input('Convert ACE files to binary? ([y]/n) ')
else:
    response = 'y'

# Convert files if requested
if not response or response.lower().startswith('y'):
    for f in ace_files:
        print('    Converting {}...'.format(f))
        openmc.data.ace.ascii_to_binary(f, f)

# ==============================================================================
# PROMPT USER TO GENERATE HDF5 LIBRARY

# Ask user to convert
if not args.batch:
    if sys.version_info[0] < 3:
        response = raw_input('Generate HDF5 library? ([y]/n) ')
    else:
        response = input('Generate HDF5 library? ([y]/n) ')
else:
    response = 'y'

# Convert files if requested
if not response or response.lower().startswith('y'):
    # Ensure 'import openmc.data' works in the openmc-ace-to-xml script
    env = os.environ.copy()
    env['PYTHONPATH'] = os.path.join(os.getcwd(), os.pardir)

    subprocess.call(['../scripts/openmc-ace-to-hdf5', '-d', 'jeff-3.2-hdf5']
                    + sorted(ace_files), env=env)
