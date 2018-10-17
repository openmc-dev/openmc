#!/usr/bin/env python

"""
Download ENDF/B-VII.1 incident neutron ACE data and incident photon ENDF data
from NNDC and convert it to an HDF5 library for use with OpenMC. This data is
used for OpenMC's regression test suite.
"""

import os
import shutil
import subprocess
import sys
import tarfile
import glob
import hashlib
import argparse
from urllib.request import urlopen

import openmc.data


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=CustomFormatter
)
parser.add_argument('-b', '--batch', action='store_true',
                    help='supresses standard in')
parser.add_argument('-n', '--neutron-only', action='store_true',
                    help='Whether to exclude photon interaction/atomic data')
parser.add_argument('--libver', choices=['earliest', 'latest'],
                    default='earliest', help="Output HDF5 versioning. Use "
                    "'earliest' for backwards compatibility or 'latest' for "
                    "performance")
args = parser.parse_args()

base_url = 'http://www.nndc.bnl.gov/endf/b7.1/aceFiles/'
files = ['ENDF-B-VII.1-neutron-293.6K.tar.gz',
         'ENDF-B-VII.1-tsl.tar.gz']
checksums = ['9729a17eb62b75f285d8a7628ace1449',
             'e17d827c92940a30f22f096d910ea186']
block_size = 16384

# ==============================================================================
# DOWNLOAD FILES FROM NNDC SITE

files_complete = []
for f in files:
    # Establish connection to URL
    url = base_url + f
    req = urlopen(url)

    # Get file size from header
    file_size = req.length
    downloaded = 0

    # Check if file already downloaded
    if os.path.exists(f):
        if os.path.getsize(f) == file_size:
            print('Skipping ' + f)
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
            status = '{0:10}  [{1:3.2f}%]'.format(
                downloaded, downloaded * 100. / file_size)
            print(status + chr(8)*len(status), end='')
        print('')
        files_complete.append(f)

# ==============================================================================
# VERIFY MD5 CHECKSUMS

print('Verifying MD5 checksums...')
for f, checksum in zip(files, checksums):
    downloadsum = hashlib.md5(open(f, 'rb').read()).hexdigest()
    if downloadsum != checksum:
        raise IOError("MD5 checksum for {} does not match. If this is your first "
                      "time receiving this message, please re-run the script. "
                      "Otherwise, please contact OpenMC developers by emailing "
                      "openmc-users@googlegroups.com.".format(f))

# ==============================================================================
# EXTRACT FILES FROM TGZ

for f in files:
    if f not in files_complete:
        continue

    # Extract files
    suffix = f[f.rindex('-') + 1:].rstrip('.tar.gz')
    with tarfile.open(f, 'r') as tgz:
        print('Extracting {}...'.format(f))
        tgz.extractall(path='nndc/' + suffix)

# Move ACE files down one level
for filename in glob.glob('nndc/293.6K/ENDF-B-VII.1-neutron-293.6K/*'):
    shutil.move(filename, 'nndc/293.6K/' + os.path.basename(filename))

# ==============================================================================
# FIX ZAID ASSIGNMENTS FOR VARIOUS S(A,B) TABLES

def fix_zaid(table, old, new):
    filename = os.path.join('nndc', 'tsl', table)
    with open(filename, 'r') as fh:
        text = fh.read()
    text = text.replace(old, new, 1)
    with open(filename, 'w') as fh:
        fh.write(text)

print('Fixing ZAIDs for S(a,b) tables')
fix_zaid('bebeo.acer', '8016', '   0')
fix_zaid('obeo.acer', '4009', '   0')

# ==============================================================================
# PROMPT USER TO DELETE .TAR.GZ FILES

# Ask user to delete
if not args.batch:
    response = input('Delete *.tar.gz files? ([y]/n) ')
else:
    response = 'y'

# Delete files if requested
if not response or response.lower().startswith('y'):
    for f in files:
        if os.path.exists(f):
            print('Removing {}...'.format(f))
            os.remove(f)

# ==============================================================================
# GENERATE HDF5 LIBRARY

# get a list of all ACE files
ace_files = sorted(glob.glob(os.path.join('nndc', '**', '*.ace*')))

# Call the ace-to-hdf5 conversion script
pwd = os.path.dirname(os.path.realpath(__file__))
ace2hdf5 = os.path.join(pwd, 'openmc-ace-to-hdf5')
subprocess.call([ace2hdf5,
                 '-d', 'nndc_hdf5',
                 '--libver', args.libver] + ace_files)

# Generate photo interaction library files
if not args.neutron_only:
    pwd = os.path.dirname(os.path.realpath(__file__))
    photo_endf = os.path.join(pwd, 'openmc-get-photon-data')
    subprocess.call([photo_endf, '-c', 'cross_sections.xml'],
                    cwd='nndc_hdf5')
