#!/usr/bin/env python3

import os
from pathlib import Path
from zipfile import ZipFile

from openmc._utils import download
import openmc.deplete


URLS = [
    'https://www.nndc.bnl.gov/endf/b7.1/zips/ENDF-B-VII.1-neutrons.zip',
    'https://www.nndc.bnl.gov/endf/b7.1/zips/ENDF-B-VII.1-decay.zip',
    'https://www.nndc.bnl.gov/endf/b7.1/zips/ENDF-B-VII.1-nfy.zip'
]

def main():
    endf_dir = os.environ.get("OPENMC_ENDF_DATA")
    if endf_dir is not None:
        endf_dir = Path(endf_dir)
    elif all(os.path.isdir(lib) for lib in ("neutrons", "decay", "nfy")):
        endf_dir = Path(".")
    else:
        for url in URLS:
            basename = download(url)
            with ZipFile(basename, 'r') as zf:
                print('Extracting {}...'.format(basename))
                zf.extractall()
        endf_dir = Path(".")

    decay_files = tuple((endf_dir / "decay").glob("*endf"))
    neutron_files = tuple((endf_dir / "neutrons").glob("*endf"))
    nfy_files = tuple((endf_dir / "nfy").glob("*endf"))

    # check files exist
    for flist, ftype in [(decay_files, "decay"), (neutron_files, "neutron"),
                         (nfy_files, "neutron fission product yield")]:
        if not flist:
            raise IOError("No {} endf files found in {}".format(ftype, endf_dir))

    chain = openmc.deplete.Chain.from_endf(decay_files, nfy_files, neutron_files)
    chain.export_to_xml('chain_endfb71.xml')


if __name__ == '__main__':
    main()
