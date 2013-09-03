#!/usr/bin/env python

import os
import glob

dirs = glob.glob('test_*')

for adir in dirs:

    os.chdir(adir)

    files = glob.glob('results.py')

    if len(files) > 0:

        files = files[0]
        with open(files, 'r') as fh:
            intxt = fh.read()
        intxt = intxt.replace('14.8E', '12.6E')
        with open(files, 'w') as fh:
            fh.write(intxt)

    os.chdir('..')
