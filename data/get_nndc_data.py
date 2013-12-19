#!/usr/bin/env python

from __future__ import print_function
import os
import os.path
import tarfile
import urllib2

baseUrl = 'http://www.nndc.bnl.gov/endf/b7.1/aceFiles/'
files = ['ENDF-B-VII.1-neutron-293.6K.tar.gz',
         'ENDF-B-VII.1-neutron-300K.tar.gz',
         'ENDF-B-VII.1-neutron-900K.tar.gz',
         'ENDF-B-VII.1-neutron-1500K.tar.gz',
         'ENDF-B-VII.1-tsl.tar.gz']
block_size = 16384

# ==============================================================================
# DOWNLOAD FILES FROM NNDC SITE

filesComplete = []
for f in files:
    # Establish connection to URL
    url = baseUrl + f
    req = urllib2.urlopen(url)

    # Get file size from header
    file_size = int(req.info().getheaders('Content-Length')[0])
    downloaded = 0

    # Check if file already downloaded
    if os.path.exists(f):
        if os.path.getsize(f) == file_size:
            print('Skipping ' + f)
            continue
        else:
            overwrite = raw_input('Overwrite {0}? ([y]/n) '.format(f))
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
    with tarfile.open(f, 'r') as tgz:
        tgz.extractall(path='nndc')

    # Give xsdir a unique name
    xsdir = 'nndc/xsdir'
    if os.path.exists(xsdir):
        os.rename(xsdir, xsdir + '_' + f.strip('.tar.gz'))
