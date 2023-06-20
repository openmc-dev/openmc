from numbers import Real

from openmc.checkvalue import check_type, check_value
from openmc import Material
from openmc.data import ELEMENT_SYMBOL


class TransferRates:
    """Class for defining continuous removals and feeds.

    Molten Salt Reactors (MSRs) benefit from continuous reprocessing,
    which removes fission products and feeds fresh fuel into the system. MSRs
    inspired the development of this class.

    An instance of this class can be passed directly to an instance of one of
    the :class:`openmc.deplete.Integrator` classes.

    .. versionadded:: 0.13.4

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
    local_mats : list of str
        All burnable material IDs being managed by a single process
    transfer_rates : dict of str to dict
        Container of transfer rates, elements and destination material
    index_transfer : Set of pair of str
        Pair of strings needed to build final matrix (destination_material, mat)
    """

    def __init__(self, operator, model):

        self.materials = model.materials
        self.burnable_mats = operator.burnable_mats
        self.local_mats = operator.local_mats

        #initialize transfer rates container dict
        self.transfer_rates = {mat: {} for mat in self.burnable_mats}
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

    def get_transfer_rate(self, material, element):
        """Return transfer rate for given material and element.

        Parameters
        ----------
        material : openmc.Material or str or int
            Depletable material
        element : str
            Element to get transfer rate value

        Returns
        -------
        transfer_rate : float
            Transfer rate value

        """
        material_id = self._get_material_id(material)
        check_value('element', element, ELEMENT_SYMBOL.values())
        return self.transfer_rates[material_id][element][0]

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
        if element in self.transfer_rates[material_id]:
            return self.transfer_rates[material_id][element][1]

    def get_elements(self, material):
        """Extract removing elements for a given material

        Parameters
        ----------
        material : openmc.Material or str or int
            Depletable material

        Returns
        -------
        elements : list
            List of elements where transfer rates exist

        """
        material_id = self._get_material_id(material)
        if material_id in self.transfer_rates:
            return self.transfer_rates[material_id].keys()

    def set_transfer_rate(self, material, elements, transfer_rate, transfer_rate_units='1/s',
                         destination_material=None):
        """Set element transfer rates in a depletable material.

        Parameters
        ----------
        material : openmc.Material or str or int
            Depletable material
        elements : list of str
            List of strings of elements that share transfer rate
        transfer_rate : float
            Rate at which elements are transferred. A positive or negative values
            set removal of feed rates, respectively.
        destination_material : openmc.Material or str or int, Optional
            Destination material to where nuclides get fed.
        transfer_rate_units : {'1/s', '1/min', '1/h', '1/d', '1/a'}
            Units for values specified in the transfer_rate argument. 's' means
            seconds, 'min' means minutes, 'h' means hours, 'a' means Julian years.

        """
        material_id = self._get_material_id(material)
        check_type('transfer_rate', transfer_rate, Real)

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

        if transfer_rate_units in ('1/s', '1/sec'):
            unit_conv = 1
        elif transfer_rate_units in ('1/min', '1/minute'):
            unit_conv = 60
        elif transfer_rate_units in ('1/h', '1/hr', '1/hour'):
            unit_conv = 60*60
        elif transfer_rate_units in ('1/d', '1/day'):
            unit_conv = 24*60*60
        elif transfer_rate_units in ('1/a', '1/year'):
            unit_conv = 365.25*24*60*60
        else:
            raise ValueError("Invalid transfer rate unit '{}'".format(transfer_rate_units))

        for element in elements:
            check_value('element', element, ELEMENT_SYMBOL.values())
            self.transfer_rates[material_id][element] = transfer_rate / unit_conv, destination_material_id
            if destination_material_id is not None:
                self.index_transfer.add((destination_material_id, material_id))
