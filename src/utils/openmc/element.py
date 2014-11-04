#!/usr/bin/env python

from openmc.checkvalue import *

class Element(object):

    def __init__(self, name='', xs=None):

        # Initialize class attributes
        self._name = ''
        self._xs = None

        # Set the Material class attributes
        self.set_name(name)

        if not xs is None:
            self.set_xs(xs)


    def __eq__(self, element2):

        # Check type
        if not isinstance(element2, Element):
            return False

        # Check name
        if self._name != element2._name:
            return False

        # Check xs
        elif self._xs != element2._xs:
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
            msg = 'Unable to set name for Element with a non-string ' \
                        'value {0}'.format(name)
            raise ValueError(msg)

        self._name = name


    def set_xs(self, xs):

        if not is_string(xs):
            msg = 'Unable to set cross-section identifier xs for Element ' \
                  'with a non-string value {0}'.format(xs)
            raise ValueError(msg)

        self._xs = xs


    def __repr__(self):

        string = 'Element    -    {0}\n'.format(self._name)
        string += '{0: <16}{1}{2}\n'.format('\tXS', '=\t', self._xs)
        return string
