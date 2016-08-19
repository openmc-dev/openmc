import sys

from openmc.checkvalue import check_type

if sys.version_info[0] >= 3:
    basestring = str


class Macroscopic(object):
    """A Macroscopic object that can be used in a material.

    Parameters
    ----------
    name : str
        Name of the macroscopic data, e.g. UO2
    xs : str
        Cross section identifier, e.g. 71c

    Attributes
    ----------
    name : str
        Name of the nuclide, e.g. UO2
    xs : str
        Cross section identifier, e.g. 71c

    """

    def __init__(self, name='', xs=None):
        # Initialize class attributes
        self._name = ''
        self._xs = None

        # Set the Macroscopic class attributes
        self.name = name

        if xs is not None:
            self.xs = xs

    def __eq__(self, other):
        if isinstance(other, Macroscopic):
            if self.name != other.name:
                return False
            elif self.xs != other.xs:
                return False
            else:
                return True
        elif isinstance(other, basestring) and other == self.name:
            return True
        else:
            return False

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash((self._name, self._xs))

    def __repr__(self):
        string = 'Nuclide    -    {0}\n'.format(self._name)
        string += '{0: <16}{1}{2}\n'.format('\tXS', '=\t', self._xs)
        return string

    @property
    def name(self):
        return self._name

    @property
    def xs(self):
        return self._xs

    @name.setter
    def name(self, name):
        check_type('name', name, basestring)
        self._name = name

    @xs.setter
    def xs(self, xs):
        check_type('cross-section identifier', xs, basestring)
        self._xs = xs
