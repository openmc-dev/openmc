#!/usr/bin/env python

from __future__ import print_function
import os
import shutil
import subprocess
import sys
import tarfile
import glob
import hashlib
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-b', '--batch', action='store_true',
                    help='supresses standard in')
args = parser.parse_args()

try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen

baseUrl = 'http://www.nndc.bnl.gov/endf/b7.1/aceFiles/'
files = ['ENDF-B-VII.1-neutron-293.6K.tar.gz',
         'ENDF-B-VII.1-tsl.tar.gz']
checksums = ['9729a17eb62b75f285d8a7628ace1449',
             'e17d827c92940a30f22f096d910ea186']
block_size = 16384

# ==============================================================================
# DOWNLOAD FILES FROM NNDC SITE

filesComplete = []
for f in files:
    # Establish connection to URL
    url = baseUrl + f
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
            print('Skipping ' + f)
            filesComplete.append(f)
            continue
        else:
            if sys.version_info[0] < 3:
                overwrite = raw_input('Overwrite {0}? ([y]/n) '.format(f))
            else:
                overwrite = input('Overwrite {0}? ([y]/n) '.format(f))
            if overwrite.lower().startswith('n'):
                continue

    # Copy file to disk
    print('Downloading {0}... '.format(f), end='')
    with open(f, 'wb') as fh:
        while True:
            chunk = req.read(block_size)
            if not chunk: break
            fh.write(chunk)
            downloaded += len(chunk)
            status = '{0:10}  [{1:3.2f}%]'.format(downloaded, downloaded * 100. / file_size)
            print(status + chr(8)*len(status), end='')
        print('')
        filesComplete.append(f)

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
    if f not in filesComplete:
        continue

    # Extract files
    suffix = f[f.rindex('-') + 1:].rstrip('.tar.gz')
    with tarfile.open(f, 'r') as tgz:
        print('Extracting {0}...'.format(f))
        tgz.extractall(path='nndc/' + suffix)

# Move ACE files down one level
for filename in glob.glob('nndc/293.6K/ENDF-B-VII.1-neutron-293.6K/*'):
    shutil.move(filename, 'nndc/293.6K/')

#===============================================================================
# EDIT GRAPHITE ZAID (6012 to 6000)

print('Changing graphite ZAID from 6012 to 6000')
graphite = os.path.join('nndc', 'tsl', 'graphite.acer')
with open(graphite) as fh:
    text = fh.read()
text = text.replace('6012', '6000', 1)
with open(graphite, 'w') as fh:
    fh.write(text)

# ==============================================================================
# PROMPT USER TO DELETE .TAR.GZ FILES

# Ask user to delete
if not args.batch:
    if sys.version_info[0] < 3:
        response = raw_input('Delete *.tar.gz files? ([y]/n) ')
    else:
        response = input('Delete *.tar.gz files? ([y]/n) ')
else:
    response = 'y'

# Delete files if requested
if not response or response.lower().startswith('y'):
    for f in files:
        if os.path.exists(f):
            print('Removing {0}...'.format(f))
            os.remove(f)

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
    # get a list of all ACE files
    ace_files = sorted(glob.glob(os.path.join('nndc', '**', '*.ace*')))

    # Ensure 'import openmc.data' works in the openmc-ace-to-xml script
    cwd = os.getcwd()
    env = os.environ.copy()
    env['PYTHONPATH'] = os.path.join(cwd, '..')

    subprocess.call(['../scripts/openmc-ace-to-hdf5', '-d', 'nndc_hdf5',
                     '--fission_energy_release', 'fission_Q_data_endfb71.h5']
                    + ace_files, env=env)
