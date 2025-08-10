from numbers import Integral

import lxml.etree as ET

import openmc.checkvalue as cv
from .mixin import EqualityMixin, IDManagerMixin
from ._xml import get_text


class TallyDerivative(EqualityMixin, IDManagerMixin):
    """A material perturbation derivative to apply to a tally.

    Parameters
    ----------
    derivative_id : int, optional
        Unique identifier for the tally derivative. If none is specified, an
        identifier will automatically be assigned
    variable : str, optional
        Accepted values are 'density', 'nuclide_density', and 'temperature'
    material : int, optional
        The perturbed material ID
    nuclide : str, optional
        The perturbed nuclide. Only needed for 'nuclide_density' derivatives.
        Ex: 'Xe135'

    Attributes
    ----------
    id : int
        Unique identifier for the tally derivative
    variable : str
        Accepted values are 'density', 'nuclide_density', and 'temperature'
    material : int
        The perturubed material ID
    nuclide : str
        The perturbed nuclide. Only needed for 'nuclide_density' derivatives.
        Ex: 'Xe135'

    """

    next_id = 1
    used_ids = set()

    def __init__(self, derivative_id=None, variable=None, material=None,
                 nuclide=None):
        # Initialize Tally class attributes
        self.id = derivative_id
        self.variable = variable
        self.material = material
        self.nuclide = nuclide

    def __repr__(self):
        string = 'Tally Derivative\n'
        string += '{: <16}=\t{}\n'.format('\tID', self.id)
        string += '{: <16}=\t{}\n'.format('\tVariable', self.variable)
        string += '{: <16}=\t{}\n'.format('\tMaterial', self.material)
        if self.variable == 'nuclide_density':
            string += '{: <16}=\t{}\n'.format('\tNuclide', self.nuclide)

        return string

    @property
    def variable(self):
        return self._variable

    @variable.setter
    def variable(self, var):
        if var is not None:
            cv.check_type('derivative variable', var, str)
            cv.check_value('derivative variable', var,
                           ('density', 'nuclide_density', 'temperature'))
        self._variable = var

    @property
    def material(self):
        return self._material

    @material.setter
    def material(self, mat):
        if mat is not None:
            cv.check_type('derivative material', mat, Integral)
        self._material = mat

    @property
    def nuclide(self):
        return self._nuclide

    @nuclide.setter
    def nuclide(self, nuc):
        if nuc is not None:
            cv.check_type('derivative nuclide', nuc, str)
        self._nuclide = nuc

    def to_xml_element(self):
        """Return XML representation of the tally derivative

        Returns
        -------
        element : lxml.etree._Element
            XML element containing derivative data

        """

        element = ET.Element("derivative")
        element.set("id", str(self.id))
        element.set("variable", self.variable)
        element.set("material", str(self.material))
        if self.variable == 'nuclide_density':
            element.set("nuclide", self.nuclide)
        return element

    @classmethod
    def from_xml_element(cls, elem):
        """Generate tally derivative from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element

        Returns
        -------
        openmc.TallyDerivative
            Tally derivative object

        """
        derivative_id = int(get_text(elem, "id"))
        variable = get_text(elem, "variable")
        material = int(get_text(elem, "material"))
        nuclide = get_text(elem, "nuclide") if variable == "nuclide_density" else None
        return cls(derivative_id, variable, material, nuclide)
