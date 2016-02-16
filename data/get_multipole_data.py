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
parser.add_argument('-b', '--batch', action = 'store_true',
                    help = 'supresses standard in')
args = parser.parse_args()

try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen

cwd = os.getcwd()
sys.path.insert(0, os.path.join(cwd, '..'))

baseUrl = 'http://web.mit.edu/smharper/Public/'
files = ['multipole_lib.tar.gz']
checksums = ['9f0307132fe5beca78b8fc7a01fb401c']
block_size = 16384

# ==============================================================================
# DOWNLOAD FILES FROM ATHENA LOCKER

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
    if not f in filesComplete:
        continue

    # Extract files
    with tarfile.open(f, 'r') as tgz:
        print('Extracting {0}...'.format(f))
        tgz.extractall(path='wmp/')

# Move data files down one level
for filename in glob.glob('wmp/multipole_lib/*'):
    shutil.move(filename, 'wmp/')
os.rmdir('wmp/multipole_lib')

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
