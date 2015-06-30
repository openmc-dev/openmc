import sys

from openmc.checkvalue import check_type

if sys.version_info[0] >= 3:
    basestring = str


class Element(object):
    """A natural element used in a material via <element>. Internally, OpenMC will
    expand the natural element into isotopes based on the known natural
    abundances.

    Parameters
    ----------
    name : str
        Chemical symbol of the element, e.g. Pu
    xs : str
        Cross section identifier, e.g. 71c

    Attributes
    ----------
    name : str
        Chemical symbol of the element, e.g. Pu
    xs : str
        Cross section identifier, e.g. 71c

    """

    def __init__(self, name='', xs=None):
        # Initialize class attributes
        self._name = ''
        self._xs = None

        # Set class attributes
        self.name = name

        if xs is not None:
            self.xs = xs

    def __eq__(self, element2):
        # Check type
        if not isinstance(element2, Element):
            return False

        # Check name and xs
        if self._name != element2._name:
            return False
        elif self._xs != element2._xs:
            return False
        else:
            return True

    def __hash__(self):
        return hash((self._name, self._xs))

    @property
    def xs(self):
        return self._xs

    @property
    def name(self):
        return self._name

    @xs.setter
    def xs(self, xs):
        check_type('cross section identifier', xs, basestring)
        self._xs = xs

    @name.setter
    def name(self, name):
        check_type('name', name, basestring)
        self._name = name

    def __repr__(self):
        string = 'Element    -    {0}\n'.format(self._name)
        string += '{0: <16}{1}{2}\n'.format('\tXS', '=\t', self._xs)
        return string
