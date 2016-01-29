from __future__ import division

#from collections import Iterable, defaultdict
#import copy
#import os
#import pickle
#import itertools
from numbers import Integral
from xml.etree import ElementTree as ET
#import sys
#
#import numpy as np
#
#from openmc import Mesh, Filter, Trigger, Nuclide
#from openmc.cross import CrossScore, CrossNuclide, CrossFilter
#from openmc.aggregate import AggregateScore, AggregateNuclide, AggregateFilter
#from openmc.filter import _FILTER_TYPES
import openmc.checkvalue as cv
#from openmc.clean_xml import *
#
#
#if sys.version_info[0] >= 3:
#    basestring = str

# "Static" variable for auto-generated TallyDerivative IDs
AUTO_TALLY_DERIV_ID = 10000

def reset_auto_tally_deriv_id():
    global AUTO_TALLY_ID
    AUTO_TALLY_DERIV_ID = 10000


class TallyDerivative(object):
    """A material perturbation derivative to apply to a tally.

    Parameters
    ----------
    derivative_id : Integral, optional
        Unique identifier for the tally derivative. If none is specified, an
        identifier will automatically be assigned

    Attributes
    ----------
    id : Integral
        Unique identifier for the tally derivative
    variable : str
        Either 'density' or 'nuclide_density'
    material : Integral
        The perutrubed material
    nuclide : str
        The perturbed nuclide. Only needed for 'nuclide_density' derivatives.
        Ex: 'Xe-135'

    """

    def __init__(self, derivative_id=None):
        # Initialize Tally class attributes
        self.id = derivative_id
        self._variable = None
        self._material = None
        self._nuclide = None

    def __deepcopy__(self, memo):
        existing = memo.get(id(self))

        # If this is the first time we have tried to copy this object, create a
        # copy
        if existing is None:
            clone = type(self).__new__(type(self))
            clone.id = self.id
            clone.variable = self.variable
            clone.material = self.material
            clone.nuclide = self.nuclide

            memo[id(self)] = clone

            return clone

        # If this object has been copied before, return the first copy made
        else:
            return existing

    def __eq__(self, other):
        if not isinstance(other, TallyDerivative):
            return False

        return (self.variable == other.variable and
                self.material == other.material and
                self.nuclide == other.nuclide)

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(repr(self))

    def __repr__(self):
        string = 'Tally Derivative\n'
        string += '{0: <16}{1}{2}\n'.format('\tID', '=\t', self.id)
        string += '{0: <16}{1}{2}\n'.format('\tVariable', '=\t',
                                            self.variable)

        if self.variable == 'density':
            string += '{0: <16}{1}{2}\n'.format('\tMaterial', '=\t',
                                                self.material)
        elif self.variable == 'nuclide_density':
            string += '{0: <16}{1}{2}\n'.format('\tMaterial', '=\t',
                                                self.material)
            string += '{0: <16}{1}{2}\n'.format('\tNuclide', '=\t',
                                                self.nuclide)
        else:
            raise RuntimeError("Encountered unrecognized differential "
                               "variable in a tally.")

        return string

    @property
    def id(self):
        return self._id

    @property
    def variable(self):
        return self._variable

    @property
    def material(self):
        return self._material

    @property
    def nuclide(self):
        return self._nuclide

    @id.setter
    def id(self, deriv_id):
        if deriv_id is None:
            global AUTO_TALLY_DERIV_ID
            self._id = AUTO_TALLY_DERIV_ID
            AUTO_TALLY_DERIV_ID += 1
        else:
            cv.check_type('tally derivative ID', deriv_id, Integral)
            cv.check_greater_than('tally derivative ID', deriv_id, 0,
                                  equality=True)
            self._id = deriv_id

    @variable.setter
    def variable(self, var):
        if var is not None:
            cv.check_type('derivative variable', var, basestring)
            if var not in ('density', 'nuclide_density'):
                raise ValueError("A tally differential variable must be either "
                     "'density' or 'nuclide_density'")
            self._variable = var

    @material.setter
    def material(self, mat):
        if mat is not None:
            cv.check_type('derivative material', mat, Integral)
            self._material = mat

    @nuclide.setter
    def nuclide(self, nuc):
        if nuc is not None:
            cv.check_type('derivative nuclide', nuc, basestring)
            self._nuclide = nuc

    def get_derivative_xml(self):
        """Return XML representation of the tally derivative

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing derivative data

        """

        element = ET.Element("derivative")
        element.set("id", str(self.id))
        element.set("variable", self.variable)
        element.set("material", str(self.material))
        if self.variable == 'nuclide_density':
            element.set("nuclide", self.nuclide)
        return element
