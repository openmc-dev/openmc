from numbers import Integral
import sys

from openmc.checkvalue import check_type

if sys.version_info[0] >= 3:
    basestring = str


class Nuclide(object):
    """A nuclide that can be used in a material.

    Parameters
    ----------
    name : str
        Name of the nuclide, e.g. U-235
    xs : str
        Cross section identifier, e.g. 71c

    Attributes
    ----------
    name : str
        Name of the nuclide, e.g. U-235
    xs : str
        Cross section identifier, e.g. 71c
    zaid : int
        1000*(atomic number) + mass number. As an example, the zaid of U-235
        would be 92235.
    scattering : 'data' or 'iso-in-lab' or None
        The type of angular scattering distribution to use

    """

    def __init__(self, name='', xs=None):
        # Initialize class attributes
        self._name = ''
        self._xs = None
        self._zaid = None
        self._scattering = None

        # Set the Material class attributes
        self.name = name

        if xs is not None:
            self.xs = xs

    def __eq__(self, other):
        if isinstance(other, Nuclide):
            if self._name != other._name:
                return False
            elif self._xs != other._xs:
                return False
            else:
                return True
        elif isinstance(other, basestring) and other == self.name:
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
        string = 'Nuclide    -    {0}\n'.format(self._name)
        string += '{0: <16}{1}{2}\n'.format('\tXS', '=\t', self._xs)
        if self._zaid is not None:
            string += '{0: <16}{1}{2}\n'.format('\tZAID', '=\t', self._zaid)
        return string

    @property
    def name(self):
        return self._name

    @property
    def xs(self):
        return self._xs

    @property
    def zaid(self):
        return self._zaid

    @property
    def scattering(self):
        return self._scattering

    @name.setter
    def name(self, name):
        check_type('name', name, basestring)
        self._name = name

    @xs.setter
    def xs(self, xs):
        check_type('cross-section identifier', xs, basestring)
        self._xs = xs

    @zaid.setter
    def zaid(self, zaid):
        check_type('zaid', zaid, Integral)
        self._zaid = zaid

    @scattering.setter
    def scattering(self, scattering):

        if not scattering in ['data', 'iso-in-lab']:
            msg = 'Unable to set scattering for Nuclide to {0} ' \
                  'which is not "data" or "iso-in-lab"'.format(scattering)
            raise ValueError(msg)

        self._scattering = scattering

    def __repr__(self):
        string = 'Nuclide    -    {0}\n'.format(self._name)
        string += '{0: <16}{1}{2}\n'.format('\tXS', '=\t', self.xs)
        if self.zaid is not None:
            string += '{0: <16}{1}{2}\n'.format('\tZAID', '=\t', self.zaid)
        if self.scattering is not None:
            string += '{0: <16}{1}{2}\n'.format('\tscattering', '=\t', 
                                                self.scattering)
        return string
