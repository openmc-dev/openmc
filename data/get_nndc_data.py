#!/usr/bin/env python

from __future__ import print_function
import os
import shutil
import subprocess
import sys
import tarfile
import glob

try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen

cwd = os.getcwd()
sys.path.append(os.path.join(cwd, '..', 'src', 'utils'))
from convert_binary import ascii_to_binary

baseUrl = 'http://www.nndc.bnl.gov/endf/b7.1/aceFiles/'
files = ['ENDF-B-VII.1-neutron-293.6K.tar.gz',
         'ENDF-B-VII.1-tsl.tar.gz']
block_size = 16384

# ==============================================================================
# DOWNLOAD FILES FROM NNDC SITE

filesComplete = []
for f in files:
    # Establish connection to URL
    url = baseUrl + f
    req = urlopen(url)

    # Get file size from header
    file_size = int(req.info().getheaders('Content-Length')[0])
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
# EXTRACT FILES FROM TGZ

for f in files:
    if not f in filesComplete:
        continue

    # Extract files
    suffix = f[f.rindex('-') + 1:].rstrip('.tar.gz')
    with tarfile.open(f, 'r') as tgz:
        print('Extracting {0}...'.format(f))
        tgz.extractall(path='nndc/' + suffix)

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
# COPY CROSS_SECTIONS.XML

print('Copying cross_sections_nndc.xml...')
shutil.copyfile('cross_sections_nndc.xml', 'nndc/cross_sections.xml')

# ==============================================================================
# PROMPT USER TO DELETE .TAR.GZ FILES

# Ask user to delete
if sys.version_info[0] < 3:
    response = raw_input('Delete *.tar.gz files? ([y]/n) ')
else:
    response = input('Delete *.tar.gz files? ([y]/n) ')

# Delete files if requested
if not response or response.lower().startswith('y'):
    for f in files:
        if os.path.exists(f):
            print('Removing {0}...'.format(f))
            os.remove(f)

# ==============================================================================
# PROMPT USER TO CONVERT ASCII TO BINARY

# Ask user to convert
if sys.version_info[0] < 3:
    response = raw_input('Convert ACE files to binary? ([y]/n) ')
else:
    response = input('Convert ACE files to binary? ([y]/n) ')

# Convert files if requested
if not response or response.lower().startswith('y'):

    # get a list of directories
    ace_dirs = glob.glob(os.path.join('nndc', '*K'))
    ace_dirs += glob.glob(os.path.join('nndc', 'tsl'))

    # loop around ace directories
    for d in ace_dirs:
        print('Coverting {0}...'.format(d))

        # get a list of files to convert
        ace_files = glob.glob(os.path.join(d, '*.ace*'))

        # convert files
        for f in ace_files:
            print('    Coverting {0}...'.format(os.path.split(f)[1]))
            ascii_to_binary(f, f)

    # Change cross_sections.xml file
    xs_file = os.path.join('nndc', 'cross_sections.xml')
    asc_str = "<filetype>ascii</filetype>"
    bin_str = "<filetype> binary </filetype>\n "
    bin_str += "<record_length> 4096 </record_length>\n "
    bin_str += "<entries> 512 </entries>"
    with open(xs_file) as fh:
        text = fh.read()
    text = text.replace(asc_str, bin_str)
    with open(xs_file, 'w') as fh:
        fh.write(text)
