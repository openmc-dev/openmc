from numbers import Integral
import numpy as np
import lxml.etree as ET

import openmc.checkvalue as cv
from .mixin import EqualityMixin, IDManagerMixin
from ._xml import get_text


class Sensitivity(EqualityMixin, IDManagerMixin):
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

    def __init__(self, sensitivity_id=None, variable=None,
                 nuclide=None, reaction=None,energy=None):
        # Initialize Tally class attributes
        self.id = sensitivity_id
        self.variable = variable
        self.nuclide = nuclide
        self.reaction = reaction
        self.energy = np.array(energy)

    def __repr__(self):
        string = 'Sensitivity\n'
        string += '{: <16}=\t{}\n'.format('\tID', self.id)
        string += '{: <16}=\t{}\n'.format('\tVariable', self.variable)
        string += '{: <16}=\t{}\n'.format('\tNuclide', self.nuclide)
        if self.variable == 'cross_section':
            string += '{: <16}=\t{}\n'.format('\tReaction', self.reaction)
            string += '{: <16}=\t{}\n'.format('\tEnergy', self.energy)

        return string

    @property
    def variable(self):
        return self._variable

    @property
    def nuclide(self):
        return self._nuclide
        
    @property
    def reaction(self):
        return self._reaction
        
    @property
    def energy(self):
        return self._energy

    @variable.setter
    def variable(self, var):
        if var is not None:
            cv.check_type('sensitivity variable', var, str)
            cv.check_value('sensitivity variable', var,
                           ('cross_section', 'multipole', 'curve_fit'))
        self._variable = var

    @nuclide.setter
    def nuclide(self, nuc):
        if nuc is not None:
            cv.check_type('sensitivity nuclide', nuc, str)
        self._nuclide = nuc
        
    @reaction.setter
    def reaction(self, rxn):
        if rxn is not None:
            cv.check_type('sensitivity reaction', rxn, str)
        self._reaction = rxn
        
    @energy.setter
    def energy(self, ene):
        self._energy = np.asarray(ene)

    def to_xml_element(self):
        """Return XML representation of the tally sensitivity

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing sensitivity data

        """

        element = ET.Element("sensitivity")
        element.set("id", str(self.id))

        element.set("variable", self.variable)
        element.set("nuclide", self.nuclide)
        
        if self.variable == 'cross_section':
            element.set("reaction", self.reaction)

            subelement = ET.SubElement(element, 'energy')
            subelement.text = ' '.join(str(e) for e in self.energy)

        return element
    
    @classmethod
    def from_xml_element(cls, elem):
        sens_id = int(elem.get("id"))
        variable = elem.get("variable")
        energy = [float(x) for x in get_text(elem, 'energy').split()]
        nuclide = elem.get("nuclide")
        reaction = elem.get("reaction") if variable == "cross_section" else None        
        return cls(sens_id, variable, nuclide, reaction, energy)