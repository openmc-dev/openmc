import warnings

import openmc.checkvalue as cv
from openmc.data import ATOMIC_NUMBER, ATOMIC_SYMBOL, zam

class Nuclide(str):
    """A nuclide that can be used in a material.

    Parameters
    ----------
    name : str
        Name of the nuclide, e.g. 'U235'

    Attributes
    ----------
    name : str
        Name of the nuclide, e.g. 'U235'
    element_name : str
        Name of the element, e.g. 'U'
    protons : int
        Number of the protons, e.g. 92
    nucleons : int
        Number of the nucleons, e.g. 235
    neutrons : int
        Number of neutrons, e.g. 143
    """

    def __new__(cls, name):
        # Initialize class attributes
        orig_name = name

        if '-' in name:
            name = name.replace('-', '')
            name = name.replace('Nat', '0')
            if name.endswith('m'):
                name = name[:-1] + '_m1'

            msg = 'OpenMC nuclides follow the GND naming convention. Nuclide ' \
                  '"{}" is being renamed as "{}".'.format(orig_name, name)
            warnings.warn(msg)

        return super().__new__(cls, name)

    @property
    def name(self):
        return self

    @property
    def element_name(self):
        return ATOMIC_SYMBOL[zam(self)[0]]

    @property
    def protons(self):
        return ATOMIC_NUMBER[self.element_name]

    @property
    def nucleons(self):
        return zam(self)[1]

    @property
    def neutrons(self):
        return zam(self)[1]-zam(self)[0]