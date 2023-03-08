from collections import OrderedDict
from numbers import Real

import numpy as np
import h5py

from openmc.deplete.abc import (_SECONDS_PER_MINUTE, _SECONDS_PER_HOUR,
                                _SECONDS_PER_DAY, _SECONDS_PER_JULIAN_YEAR)
from openmc.checkvalue import check_type, check_value
from openmc import Materials, Material
from openmc.data import ELEMENT_SYMBOL

class MsrContinuous:
    """Class for defining continuous removals and feeds.

    Molten Salt Reactors (MSRs) benefit from continuous reprocessing,
    which removes fission products and feeds fresh fuel into the system. MSRs
    inspired the development of this class.

    An instance of this class can be passed directly to an instance of one of
    the :class:`openmc.deplete.Integrator` classes.

    .. versionadded:: 0.13.3
    Parameters
    ----------
    operator : openmc.TransportOperator
        Depletion operator
    model : openmc.Model
        OpenMC model containing materials and geometry. If using
        :class:`openmc.deplete.CoupledOperator`, the model must also contain
        a :class:`opnemc.Settings` object.

    Attributes
    ----------
    burnable_mats : list of str
        All burnable material IDs.
    removal_rates : OrderedDict of str and OrderedDict
        Container of removal rates, elements and destination material
    index_transfer : Set of pair of str
        Pair of strings needed to build final matrix (destination_material, mat)
    """

    def __init__(self, operator, model):

        self.operator = operator
        self.materials = model.materials
        self.burnable_mats = operator.burnable_mats

        #initialize removal rates container dict
        self.removal_rates = OrderedDict((mat, OrderedDict()) for mat in \
                                          self.burnable_mats)
        self.index_transfer = set()

    def _get_material_id(self, val):
        """Helper method for getting material id from Material obj or name.

        Parameters
        ----------
        val : openmc.Material or str or int representing material name/id

        Returns
        -------
        material_id : str

        """
        if isinstance(val, Material):
            check_value('Depeletable Material', str(val.id), self.burnable_mats)
            val = val.id

        elif isinstance(val, str):
            if val.isnumeric():
                check_value('Material ID', str(val), self.burnable_mats)
            else:
                check_value('Material name', val,
                        [mat.name for mat in self.materials if mat.depletable])
                val = [mat.id for mat in self.materials if mat.name == val][0]

        elif isinstance(val, int):
            check_value('Material ID', str(val), self.burnable_mats)

        return str(val)

    def get_removal_rate(self, material, element):
        """Return removal rate for given material and element.

        Parameters
        ----------
        material : openmc.Material or str or int
            Depletable material
        element : str
            Element to get removal rate value

        Returns
        -------
        removal_rate : float
            Removal rate value

        """
        material_id = self._get_material_id(material)
        check_value('element', element, ELEMENT_SYMBOL.values())
        return self.removal_rates[material_id][element][0]

    def get_destination_material(self, material, element):
        """Return destination (or transfer) material for given material and
        element, if defined.

        Parameters
        ----------
        material : openmc.Material or str or int
            Depletable material
        element : str
            Element that gets transferred to another material.

        Returns
        -------
        destination_material_id : str
            Depletable material ID to where the element gets transferred

        """
        material_id = self._get_material_id(material)
        check_value('element', element, ELEMENT_SYMBOL.values())
        if element in self.removal_rates[material_id]:
            return self.removal_rates[material_id][element][1]

    def get_elements(self, material):
        """Extract removing elements for a given material

        Parameters
        ----------
        material : openmc.Material or str or int
            Depletable material

        Returns
        -------
        elements : list
            List of elements where removal rates exist

        """
        material_id = self._get_material_id(material)
        if material_id in self.removal_rates.keys():
            return self.removal_rates[material_id].keys()

    def set_removal_rate(self, material, elements, removal_rate, removal_rate_units='1/s',
                         destination_material=None):
        """Set element removal rates in a depletable material.

        Parameters
        ----------
        material : openmc.Material or str or int
            Depletable material
        elements : list of str
            List of strings of elements that share removal rate
        removal_rate : float
            Removal rate
        destination_material : Openmc.Material or str or int, Optional
            Destination material to where nuclides get fed.
        removal_rate_units: {'1/s', '1/min', '1/h', '1/d', '1/a'}
            Units for values specified in the removal_rate argument. 's' means
            seconds, 'min' means minutes, 'h' means hours, 'a' means Julian years.

        """
        material_id = self._get_material_id(material)
        check_type('removal_rate', removal_rate, Real)

        if destination_material is not None:
            destination_material_id = self._get_material_id(destination_material)
            if len(self.burnable_mats) > 1:
                check_value('destination_material', str(destination_material_id),
                            self.burnable_mats)
            else:
                raise ValueError(f'Transfer to material {destination_material_id} '\
                        'is set, but there is only one depletable material')
        else:
            destination_material_id = None

        if removal_rate_units in ('1/s', '1/sec'):
            unit_conv = 1
        elif removal_rate_units in ('1/min', '1/minute'):
            unit_conv = _SECONDS_PER_MINUTE
        elif removal_rate_units in ('1/h', '1/hr', '1/hour'):
            unit_conv = _SECONDS_PER_HOUR
        elif removal_rate_units in ('1/d', '1/day'):
            unit_conv = _SECONDS_PER_DAY
        elif removal_rate_units in ('1/a', '1/year'):
            unit_conv = _SECONDS_PER_JULIAN_YEAR
        else:
            raise ValueError("Invalid removal rate unit '{}'".format(removal_rate_units))

        for element in elements:
            check_value('element', element, ELEMENT_SYMBOL.values())
            self.removal_rates[material_id][element] = removal_rate / unit_conv, destination_material_id
            if destination_material_id is not None:
                self.index_transfer.add((destination_material_id, material_id))
