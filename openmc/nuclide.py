import warnings

from six import string_types

import openmc.checkvalue as cv


class Nuclide(object):
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

    def __init__(self, name=''):
        # Initialize class attributes
        self._name = ''

        # Set the Material class attributes
        self.name = name

    def __eq__(self, other):
        if isinstance(other, Nuclide):
            if self.name != other.name:
                return False
            else:
                return True
        elif isinstance(other, string_types) and other == self.name:
            return True
        else:
            return False

    def __ne__(self, other):
        return not self == other

    def __gt__(self, other):
        return repr(self) > repr(other)

    def __lt__(self, other):
        return not self > other

    def __hash__(self):
        return hash(repr(self))

    def __repr__(self):
        return 'Nuclide    -    {0}\n'.format(self._name)

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        cv.check_type('name', name, string_types)
        self._name = name

        if '-' in name:
            self._name = name.replace('-', '')
            self._name = self._name.replace('Nat', '0')
            if self._name.endswith('m'):
                self._name = self._name[:-1] + '_m1'

            msg = 'OpenMC nuclides follow the GND naming convention. Nuclide ' \
                  '"{}" is being renamed as "{}".'.format(name, self._name)
            warnings.warn(msg)
