#!/usr/bin/env python

import glob
import os
from zipfile import ZipFile

import requests
from tqdm import tqdm
import openmc.deplete


urls = [
    'http://www.nndc.bnl.gov/endf/b7.1/zips/ENDF-B-VII.1-neutrons.zip',
    'http://www.nndc.bnl.gov/endf/b7.1/zips/ENDF-B-VII.1-decay.zip',
    'http://www.nndc.bnl.gov/endf/b7.1/zips/ENDF-B-VII.1-nfy.zip'
]


def download_file(url):
    response = requests.get(url, stream=True)
    filesize = int(response.headers.get('content-length'))

    # Check if file already downloaded
    basename = url.split('/')[-1]
    if os.path.exists(basename):
        if os.path.getsize(basename) == filesize:
            return basename
        else:
            overwrite = input('Overwrite {}? ([y]/n) '.format(basename))
            if overwrite.lower().startswith('n'):
                return basename

    with open(basename, 'wb') as f:
        with tqdm(desc='Downloading {}'.format(basename),
                  total=filesize, unit='B', unit_scale=True) as pbar:
            for i, chunk in enumerate(response.iter_content(chunk_size=4096)):
                pbar.update(4096)
                if chunk:
                    f.write(chunk)

    return basename


def main():
    for url in urls:
        basename = download_file(url)
        with ZipFile(basename, 'r') as zf:
            print('Extracting {}...'.format(basename))
            zf.extractall()

    decay_files = glob.glob(os.path.join('decay', '*.endf'))
    nfy_files = glob.glob(os.path.join('nfy', '*.endf'))
    neutron_files = glob.glob(os.path.join('neutrons', '*.endf'))

    chain = openmc.deplete.DepletionChain.from_endf(decay_files, nfy_files, neutron_files)
    chain.xml_write('chain_endfb71.xml')


if __name__ == '__main__':
    main()
