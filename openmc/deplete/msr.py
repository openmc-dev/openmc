from collections import OrderedDict
from openmc import Materials, Material

class MsrContinuous:
    """Class defining Molten salt reactor (msr) elements (fission products)
    removal, based on removal rates and cycle time concepts.

    Parameters
    ----------
    local_mats : openmc.Material
        openmc.Material

    Attributes
    ----------
    local_mats : openmc.Material
        All local material defined in the model.
    index_burn : list of str
        List of depletable material id
    ordr_burn : OrderedDict of int and int
        OrderedDict of depletable material id and enuemerated indeces
    removal_rates : OrderedDict of str and OrderedDict
        OrderedDict of depletable material id and OrderedDict to fill
    """

    def __init__(self, local_mats):

        if not isinstance(local_mats, Materials):
            raise ValueError(f'{local_mats} is not a valid openmc Material')
        else:
            self.local_mats = local_mats

        self.index_burn = [mat.id for mat in self.local_mats if mat.depletable]
        self.index_msr = [(i, i) for i in range(self.n_burn)]

        self.ord_burn = self._order_burn_mats()
        self.removal_rates = self._initialize_removal_rates()

    @property
    def n_burn(self):
        return len(self.index_burn)

    def _order_burn_mats(self):
        """Order depletable material id
        Returns
        ----------
        OrderedDict of int and int
            OrderedDict of depletable material id and enuemrated indeces
        """
        return OrderedDict((int(id), i) for i, id in enumerate(self.index_burn))

    def _initialize_removal_rates(self):
        """Initialize removal rates container
        Returns
        ----------
        OrderedDict of str and OrderedDict
            OrderedDict of depletable material id and OrderedDict to fill
        """
        return OrderedDict((id, OrderedDict()) for id in self.index_burn)

    def _get_mat_index(self, mat):
        """Helper method for getting material index"""
        if isinstance(mat, Material):
            mat = str(mat.id)
        return self.index_burn[mat] if isinstance(mat, str) else mat

    def transfer_matrix_index(self):
        """Get transfer matrices indeces
        Returns
        ----------
        list of tuples :
            List of tuples pairs to order the transfer matrices
            when building the coupled matrix
        """
        transfer_index = OrderedDict()
        for id, val in self.removal_rates.items():
            if val:
                for elm, [tr, mat] in val.items():
                    if mat is not None:
                        j = self.ord_burn[id]
                        i = self.ord_burn[mat]
                        transfer_index[(i,j)] = None
        return list(transfer_index.keys())

    def get_removal_rate(self, mat, element):
        """Extract removal rates
        Parameters
        ----------
        mat : Openmc,Material or int
            Depletable material
        element : str
            Element to extract removal rate value
        Returns:
        ----------
        float :
            Removal rate value
        """
        mat = self._get_mat_index(mat)
        return self.removal_rates[mat][element][0]

    def get_destination_mat(self, mat, element):
        """Extract destination material
        Parameters
        ----------
        mat : Openmc,Material or int
            Depletable material
        element : str
            Element to extract removal rate value
        Returns:
        ----------
        str :
            Depletable material id to where elements get transferred
        """
        mat = self._get_mat_index(mat)
        return self.removal_rates[mat][element][1]

    def get_elements(self, mat):
        """Extract removing elements for a given material
        Parameters
        ----------
        mat : Openmc,Material or int
            Depletable material
        Returns:
        ----------
        list :
            list of elements
        """
        mat = self._get_mat_index(mat)
        elements=[]
        for k, v in self.removal_rates.items():
            if k == mat:
                for elm, _ in v.items():
                    elements.append(elm)
        return elements

    def set_removal_rate(self, mat, elements, removal_rate, dest_mat=None,
                         units='1/s'):
        """Set removal rates to depletable material
        Parameters
        ----------
        mat : Openmc,Material or int
            Depletable material to where add removal rates
        elements : list[str]
            List of strings of elements
        removal_rate : float
            Removal rate coefficient
        dest_mat : Openmc,Material or int, Optional
            Destination material if transfer or elements tracking is wanted
        units: str, optional
            Removal rates units (not set yet)
            Default : '1/s'
        """
        mat = self._get_mat_index(mat)
        if dest_mat is not None:
            dest_mat = self._get_mat_index(dest_mat)
        for element in elements:
            self.removal_rates[mat][element] = [removal_rate, dest_mat]
