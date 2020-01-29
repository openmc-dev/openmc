#!/usr/bin/env python

from urllib.parse import urlencode
from urllib.request import urlopen
from lxml import html

import numpy as np
import h5py
from openmc.data import ATOMIC_SYMBOL


base_url = 'https://physics.nist.gov/cgi-bin/Star/e_table-t.pl'
energies = np.logspace(-3, 3, 200)
data = {'matno': '', 'Energies': '\n'.join(str(x) for x in energies)}
columns = {1: 's_collision', 2: 's_radiative'}

# ==============================================================================
# SCRAPE DATA FROM ESTAR SITE AND GENERATE STOPPING POWER HDF5 FILE

print('Generating stopping_powers.h5...')

with h5py.File('stopping_powers.h5', 'w') as f:

    # Write energies
    f.create_dataset('energy', data=energies)

    # Look over atomic number; ESTAR only goes up to Z=98 (Californium)
    for Z in range(1, 99):
        print('Processing {} data...'.format(ATOMIC_SYMBOL[Z]))

        # Update form-encoded data to send in POST request for this element
        data['matno'] = '{:03}'.format(Z)
        payload = urlencode(data).encode("utf-8")

        # Retrieve data from ESTAR site
        with urlopen(url=base_url, data=payload) as response:
            r = response.read()

        # Remove text and reformat data -- omit first 12 and last 5 lines to get
        # only data in table
        r = html.fromstring(r).xpath('//pre//text()')
        values = np.fromstring(' '.join(r[12:-5]), sep=' ').reshape((-1, 5)).T

        # Create group for this element
        group = f.create_group('{:03}'.format(Z))

        # Write the mean excitation energy
        attributes = np.fromstring(r[3], sep=' ')
        group.attrs['I'] = attributes[2]

        # Write collision and radiative stopping powers
        for i in columns:
            group.create_dataset(columns[i], data=values[i])
