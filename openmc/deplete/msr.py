from collections import OrderedDict

from openmc import Materials, Material

class MsrContinuous:
    """Class defining a msr continuous object.

    .. versionadded:: 0.13.2

    Parameters
    ----------
    local_mats : openmc.Material
        opoenmc.Material

    Attributes
    ----------
    local_mats : str or None
        Name of nuclide.
    index_burn : float or None
        Half life of nuclide in [s].
    ordr_burn : float
        Energy deposited from decay in [eV].
    removal_terms : int
        Number of decay pathways.
    """

    def __init__(self, local_mats):

        if not isinstance(local_mats, Materials):
            raise ValueError(f'{local_mats} is not a valid openmc Material')

        self.index_burn = [mat.id for mat in local_mats if mat.depletable]
        self.n_burn = len(self.index_burn)
        self.ord_burn = OrderedDict((int(id), i) for i, id in enumerate(self.index_burn))
        self.removal_terms = OrderedDict((id, OrderedDict()) for id in self.index_burn)
        self.index_msr = [(i, i) for i in range(self.n_burn)]

    def _get_mat_index(self, mat):
        """Helper method for getting material index"""
        if isinstance(mat, Material):
            mat = str(mat.id)
        return self.index_burn[mat] if isinstance(mat, str) else mat

    def get_transfer_index(self):
        transfer_index = OrderedDict()
        for id, val in self.removal_terms.items():
            if val:
                for elm, [tr, mat] in val.items():
                    if mat is not None:
                        j = self.ord_burn[id]
                        i = self.ord_burn[mat]
                        transfer_index[(i,j)] = None
        return list(transfer_index.keys())

    def get_removal_rate(self, mat, element):
        mat = self._get_mat_index(mat)
        return self.removal_terms[mat][element][0]

    def get_destination_mat(self, mat, element):
        mat = self._get_mat_index(mat)
        return self.removal_terms[mat][element][1]

    def get_elements(self, mat):
        mat = self._get_mat_index(mat)
        elements=[]
        for k, v in self.removal_terms.items():
            if k == mat:
                for elm, _ in v.items():
                    elements.append(elm)
        return elements

    def set_removal_rate(self, mat, elements, removal_rate, dest_mat=None, units='1/s'):
        mat = self._get_mat_index(mat)
        if dest_mat is not None:
            dest_mat = self._get_mat_index(dest_mat)
        for element in elements:
            self.removal_terms[mat][element] = [removal_rate, dest_mat]
