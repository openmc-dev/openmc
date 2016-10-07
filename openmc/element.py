import re
import sys

from six import string_types

import openmc
from openmc.checkvalue import check_type, check_length
from openmc.data import NATURAL_ABUNDANCE


class Element(object):
    """A natural element used in a material via <element>. Internally, OpenMC will
    expand the natural element into isotopes based on the known natural
    abundances.

    Parameters
    ----------
    name : str
        Chemical symbol of the element, e.g. Pu

    Attributes
    ----------
    name : str
        Chemical symbol of the element, e.g. Pu
    scattering : {'data', 'iso-in-lab', None}
        The type of angular scattering distribution to use

    """

    def __init__(self, name=''):
        # Initialize class attributes
        self._name = ''
        self._scattering = None

        # Set class attributes
        self.name = name

    def __eq__(self, other):
        if isinstance(other, Element):
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
        string = 'Element    -    {0}\n'.format(self._name)
        if self.scattering is not None:
            string += '{0: <16}{1}{2}\n'.format('\tscattering', '=\t',
                                                self.scattering)

        return string

    @property
    def name(self):
        return self._name

    @property
    def scattering(self):
        return self._scattering

    @name.setter
    def name(self, name):
        check_type('element name', name, string_types)
        check_length('element name', name, 1, 2)
        self._name = name

    @scattering.setter
    def scattering(self, scattering):

        if not scattering in ['data', 'iso-in-lab']:
            msg = 'Unable to set scattering for Element to {0} ' \
                  'which is not "data" or "iso-in-lab"'.format(scattering)
            raise ValueError(msg)

        self._scattering = scattering

    def expand(self):
        """Expand natural element into its naturally-occurring isotopes.

        Returns
        -------
        isotopes : list
            Naturally-occurring isotopes of the element. Each item of the list
            is a tuple consisting of an openmc.Nuclide instance and the natural
            abundance of the isotope.

        """

        isotopes = []
        for isotope, abundance in sorted(NATURAL_ABUNDANCE.items()):
            if re.match(r'{}\d+'.format(self.name), isotope):
                nuc = openmc.Nuclide(isotope)
                isotopes.append((nuc, abundance))
        return isotopes
