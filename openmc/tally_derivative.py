from __future__ import division

import sys
from numbers import Integral
from xml.etree import ElementTree as ET

from six import string_types

import openmc.checkvalue as cv
from openmc.mixin import EqualityMixin


# "Static" variable for auto-generated TallyDerivative IDs
AUTO_TALLY_DERIV_ID = 10000

def reset_auto_tally_deriv_id():
    global AUTO_TALLY_ID
    AUTO_TALLY_DERIV_ID = 10000


class TallyDerivative(EqualityMixin):
    """A material perturbation derivative to apply to a tally.

    Parameters
    ----------
    derivative_id : Integral, optional
        Unique identifier for the tally derivative. If none is specified, an
        identifier will automatically be assigned
    variable : str, optional
        Accepted values are 'density', 'nuclide_density', and 'temperature'
    material : Integral, optional
        The perturubed material ID
    nuclide : str, optional
        The perturbed nuclide. Only needed for 'nuclide_density' derivatives.
        Ex: 'Xe135'

    Attributes
    ----------
    id : Integral
        Unique identifier for the tally derivative
    variable : str
        Accepted values are 'density', 'nuclide_density', and 'temperature'
    material : Integral
        The perturubed material ID
    nuclide : str
        The perturbed nuclide. Only needed for 'nuclide_density' derivatives.
        Ex: 'Xe135'

    """

    def __init__(self, derivative_id=None, variable=None, material=None,
                 nuclide=None):
        # Initialize Tally class attributes
        self.id = derivative_id
        self.variable = variable
        self.material = material
        self.nuclide = nuclide

    def __hash__(self):
        return hash(repr(self))

    def __repr__(self):
        string = 'Tally Derivative\n'
        string += '{: <16}=\t{}\n'.format('\tID', self.id)
        string += '{: <16}=\t{}\n'.format('\tVariable', self.variable)

        if self.variable == 'density':
            string += '{: <16}=\t{}\n'.format('\tMaterial', self.material)
        elif self.variable == 'nuclide_density':
            string += '{: <16}=\t{}\n'.format('\tMaterial', self.material)
            string += '{: <16}=\t{}\n'.format('\tNuclide', self.nuclide)
        elif self.variable == 'temperature':
            string += '{: <16}=\t{}\n'.format('\tMaterial', self.material)

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
            cv.check_type('derivative variable', var, string_types)
            cv.check_value('derivative variable', var,
                           ('density', 'nuclide_density', 'temperature'))
        self._variable = var

    @material.setter
    def material(self, mat):
        if mat is not None:
            cv.check_type('derivative material', mat, Integral)
        self._material = mat

    @nuclide.setter
    def nuclide(self, nuc):
        if nuc is not None:
            cv.check_type('derivative nuclide', nuc, string_types)
        self._nuclide = nuc

    def to_xml_element(self):
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
