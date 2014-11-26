#!/usr/bin/env python

from openmc.checkvalue import *


class Nuclide(object):

    def __init__(self, name='', xs=None):

        # Initialize class attributes
        self._name = ''
        self._xs = None
        self._zaid = None

        # Set the Material class attributes
        self.set_name(name)

        if not xs is None:
            self.set_xs(xs)


    def __eq__(self, nuclide2):

        # Check type
        if not isinstance(nuclide2, Nuclide):
            return False

        # Check name
        elif self._name != nuclide2._name:
            return False

        # Check xs
        elif self._xs != nuclide2._xs:
            return False

        else:
            return True


    def __hash__(self):
        hashable = list()
        hashable.append(self._name)
        hashable.append(self._xs)
        return hash(tuple(hashable))


    def set_name(self, name):

        if not is_string(name):
            msg = 'Unable to set name for Nuclide with a non-string ' \
                        'value {0}'.format(name)
            raise ValueError(msg)

        self._name = name


    def set_xs(self, xs):

        if not is_string(xs):
            msg = 'Unable to set cross-section identifier xs for Nuclide ' \
                  'with a non-string value {0}'.format(xs)
            raise ValueError(msg)

        self._xs = xs


    def set_zaid(self, zaid):

        if not is_integer(zaid):
            msg = 'Unable to set zaid for Nuclide ' \
                  'with a non-integer {0}'.format(zaid)
            raise ValueError(msg)

        self._zaid = zaid


    def __repr__(self):

        string = 'Nuclide    -    {0}\n'.format(self._name)
        string += '{0: <16}{1}{2}\n'.format('\tXS', '=\t', self._xs)
        if self._zaid is not None:
          string += '{0: <16}{1}{2}\n'.format('\tZAID', '=\t', self._zaid)
        return string