import warnings

import openmc.checkvalue as cv


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
