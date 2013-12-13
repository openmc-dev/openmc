#!/usr/bin/env python

from __future__ import print_function
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
    print('Downloading {0}... '.format(f), end='')
    req = urllib2.urlopen(url)

    # Get file size from header
    file_size = int(req.info().getheaders("Content-Length")[0])
    downloaded = 0

    # Copy file to disk
    with open(f, 'wb') as fh:
        while True:
            chunk = req.read(block_size)
            if not chunk: break
            fh.write(chunk)
            downloaded += len(chunk)
            status = "{0:10}  [{1:3.2f}%]".format(downloaded, downloaded * 100. / file_size)
            print(status + chr(8)*len(status), end='')
        filesComplete.append(f)

# ==============================================================================
# EXTRACT FILES FROM TGZ

for f in files:
    if not f in filesComplete:
        continue

    with tarfile.open(f, 'r') as tgz:
        tgz.extractall(path='nndc')
