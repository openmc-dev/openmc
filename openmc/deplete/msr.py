from abc import ABC, abstractmethod
from collections import OrderedDict
from copy import deepcopy
from warnings import warn
from numbers import Real

import numpy as np
import h5py

from openmc.checkvalue import check_type, check_value, check_less_than, \
check_iterable_type, check_length
from openmc import Materials, Material, Cell
from openmc.search import _SCALAR_BRACKETED_METHODS, search_for_keff
from openmc.data import atomic_mass, AVOGADRO, ELEMENT_SYMBOL
import openmc.lib

class MsrContinuous:
    """Class defining Molten salt reactor (msr) elements (e.g. fission products)
    continuous removal and transfer, based on removal rates and cycle time
    theory.

    An instance of this class can be passed directly to an instance of the
    integrator class, such as :class:`openmc.deplete.CECMIntegrator`.

    Parameters
    ----------
    operator : openmc.Operator
        OpenMC operator object
    model : openmc.Model
        OpenMC Model object

    Attributes
    ----------
    burn_mats : list of str
        All burnable material IDs.
    removal_rates : OrderedDict of str and OrderedDict
        Container of removal rates, elements and destination material
    index_transfer : Set of pair of str
        Pair of strings needed to build final depletion matrix (dest_mat, mat)
    """

    def __init__(self, operator, model):

        self.operator = operator
        self.materials = model.materials
        self.burn_mats = operator.burnable_mats

        #initialize removal rates container dict
        self.removal_rates = OrderedDict((mat, OrderedDict()) for mat in \
                                          self.burn_mats)
        self.index_transfer = set()

    def _get_mat_id(self, val):
        """Helper method for getting material id from Material obj or name.

        Parameters
        ----------
        val : Openmc,Material or str or int representing material name/id

        Returns
        ----------
        id : str
            Material id
        """
        if isinstance(val, Material):
            check_value('Material depletable', str(val.id), self.burn_mats)
            val = val.id

        elif isinstance(val, str):
            if val.isnumeric():
                check_value('Material id', str(val), self.burn_mats)
            else:
                check_value('Material name', val,
                        [mat.name for mat in self.materials if mat.depletable])
                val = [mat.id for mat in self.materials if mat.name == val][0]

        elif isinstance(val, int):
            check_value('Material id', str(val), self.burn_mats)

        return str(val)

    def get_removal_rate(self, mat, element):
        """Return removal rate for given material and element.

        Parameters
        ----------
        mat : Openmc,Material or str or int
            Depletable material
        element : str
            Element to get removal rate value

        Returns
        -------
        removal_rate : float
            Removal rate value
        """
        mat = self._get_mat_id(mat)
        check_value('Element', element, ELEMENT_SYMBOL.values())
        return self.removal_rates[mat][element][0]

    def get_destination_mat(self, mat, element):
        """Return destination (or transfer) material for given material and
        element, if defined.
        Parameters
        ----------
        mat : Openmc,Material or str or int
            Depletable material
        element : str
            Element that get transferred to another material.
        Returns:
        ----------
        destination_mat : str
            Depletable material id to where the aelement get transferred
        """
        mat = self._get_mat_id(mat)
        check_value('Element', element, ELEMENT_SYMBOL.values())
        if element in self.removal_rates[mat]:
            return self.removal_rates[mat][element][1]

    def get_elements(self, mat):
        """Extract removing elements for a given material
        Parameters
        ----------
        mat : Openmc,Material or str or int
            Depletable material
        Returns:
        ----------
        elements : list
            List of elements where a removal rates exist
        """
        mat = self._get_mat_id(mat)
        if mat in self.removal_rates.keys():
            return self.removal_rates[mat].keys()

    def set_removal_rate(self, mat, elements, removal_rate, units='1/s',
                         dest_mat=None):
        """Set removal rate to elements in a depletable material.
        Parameters
        ----------
        mat : Openmc,Material or str or int
            Depletable material
        elements : list[str]
            List of strings of elements that share removal rate
        removal_rate : float
            Removal rate value in [1/sec]
        dest_mat : Openmc,Material or str or int, Optional
            Destination (transfer) material if transfer or elements tracking
            is enabled.
        units: str, optional
            Removal rates units
            Default : '1/s'
        """
        mat = self._get_mat_id(mat)
        check_type('removal_rate', removal_rate, Real)

        if dest_mat is not None:
            dest_mat = self._get_mat_id(dest_mat)
            #prevent for setting tranfert to material if not set as depletable
            if len(self.burn_mats) > 1:
                check_value('transfert to material', str(dest_mat),
                            self.burn_mats)
            else:
                raise ValueError(f'Transfer to material {dest_mat} is set '\
                        'but there is only one depletable material')

        if units != '1/s':
            check_value('Units', units, ['1/h', '1/d'])
            if units == '1/h':
                unit_conv = 1/3600
            elif units == '1/d':
                unit_conv = 1/86400
        else:
            unit_conv = 1
        for element in elements:
            check_value('Element', element, ELEMENT_SYMBOL.values())
            self.removal_rates[mat][element] = removal_rate*unit_conv, dest_mat
            if dest_mat is not None:
                self.index_transfer.add((dest_mat, mat))
