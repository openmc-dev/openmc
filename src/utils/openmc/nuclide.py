#!/usr/bin/env python

from openmc.checkvalue import *


class Nuclide(object):

    def __init__(self, name='', xs=None):

        # Initialize class attributes
        self._name = ''
        self._xs = None

        # Set the Material class attributes
        self.setName(name)

        if not xs is None:
            self.setXS(xs)


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


    def setName(self, name):

        if not is_string(name):
            msg = 'Unable to set name for Nuclide with a non-string ' \
                        'value {0}'.format(name)
            raise ValueError(msg)

        self._name = name


    def setXS(self, xs):

        if not is_string(xs):
            msg = 'Unable to set cross-section identifier xs for Nuclide ' \
                  'with a non-string value {0}'.format(xs)
            raise ValueError(msg)

        self._xs = xs


    def __repr__(self):

        string = 'Nuclide    -    {0}\n'.format(self._name)
        string += '{0: <16}{1}{2}\n'.format('\tXS', '=\t', self._xs)
        return string
