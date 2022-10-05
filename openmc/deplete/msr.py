from abc import ABC, abstractmethod
from collections import OrderedDict
from copy import deepcopy

import numpy as np

from openmc import Materials, Material
from openmc.search import search_for_keff
from openmc.data import atomic_mass, AVOGADRO
from openmc.lib import init_geom

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
    burn_mats : list of str
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

        self.burn_mats = [mat.id for mat in self.local_mats if mat.depletable]
        self.index_msr = [(i, i) for i in range(self.n_burn)]

        self.ord_burn = self._order_burn_mats()
        self.removal_rates = self._initialize_removal_rates()

    @property
    def n_burn(self):
        return len(self.burn_mats)

    def _get_mat_index(self, mat):
        """Helper method for getting material index"""
        if isinstance(mat, Material):
            mat = str(mat.id)
        return self.burn_mats[mat] if isinstance(mat, str) else mat

    def _order_burn_mats(self):
        """Order depletable material id
        Returns
        ----------
        OrderedDict of int and int
            OrderedDict of depletable material id and enuemrated indeces
        """
        return OrderedDict((int(id), i) for i, id in enumerate(self.burn_mats))

    def _initialize_removal_rates(self):
        """Initialize removal rates container
        Returns
        ----------
        OrderedDict of str and OrderedDict
            OrderedDict of depletable material id and OrderedDict to fill
        """
        return OrderedDict((id, OrderedDict()) for id in self.burn_mats)

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


class MsrBatchwise(ABC):
    """
    """
    def __init__(self, operator, model, range=None, bracketed_method='brentq',
                 tol=0.01, target=1.0, atom_density_limit=0.0,
                 refine_search=True, refuel=False, start_param=0.0):

        self.operator = operator
        self.burn_mats = operator.burnable_mats

        # Select nuclides with data that are also in the chain
        self.burn_nucs = [nuc.name for nuc in operator.chain.nuclides
                               if nuc.name in operator.nuclides_with_data]
        self.model = model
        self.geometry = model.geometry
        self.materials = model.materials

        if len(range) != 3:
            raise ValueError(f'{range} lenght must be 3')
        else:
            self.range = range

        self.bracketed_method = bracketed_method
        self.tol = tol
        self.target = target
        self.atom_density_limit = atom_density_limit
        self.refine_search = refine_search
        self.refuel = refuel
        self.start_param = start_param

    def _order_vector(self, x):
        """
        Order concentrations vector into a list of nuclides dictionaries.
        Parameters
        ------------
        x : atoms vector
        """
        ord_nucs = list()
        for idx, burn_id in enumerate(self.burn_mats):
            nucs = OrderedDict()
            for id_nuc, nuc in enumerate(self.operator.number.burnable_nuclides):
                nucs[nuc] = x[idx][id_nuc]
            ord_nucs.append(nucs)
        return ord_nucs

    @abstractmethod
    def msr_criticality_search(self, parametric_model):
        pass

    @abstractmethod
    def _build_parametric_model(self, param):
        pass

    def finalize(self):
        """
        """
        pass

class MsrBatchwiseGeom(MsrBatchwise):
    """
    """
    def __init__(self, operator, model, geom_id, type='surface', range=None,
                    bracketed_method='brentq', tol=0.01, target=1.0,
                    atom_density_limit=0.0, refine_search=True, refuel=False,
                    start_param=0.0):

        super().__init__(operator, model, range, bracketed_method, tol, target,
                         atom_density_limit, refine_search, refuel, start_param)

        self.geom_id = geom_id
        if type == 'surface':
            self.coeff = self._extract_geom_coeff()

    def _extract_geom_coeff(self):
        for surf in self.geometry.get_all_surfaces().items():
            if surf[1].id == self.geom_id:
                keys = list(surf[1].coefficients.keys())
                if len(keys) == 1:
                    coeff = getattr(surf[1], keys[0])
                else:
                    msg=(f'Surface coefficients {keys} are more than one')
                    raise Exception(msg)
        return coeff

    def _set_geom_coeff(self, res):
        """
        """
        for surf in self.geometry.get_all_surfaces().items():
            if surf[1].id == self.geom_id:
                keys = list(surf[1].coefficients.keys())
                if len(keys) == 1:
                    setattr(surf[1], keys[0], res)
                else:
                    msg = (f'Surface coefficients {keys} are more than one')
                    raise Exception(msg)

    def _normalize_nuclides(self, x):
        """
        Remove
        """
        materials = deepcopy(self.materials)
        nucs = super()._order_vector(x)

        for idx, burn_mat in enumerate(self.burn_mats):
            atoms_gram_per_mol = 0
            for nuc, val in nucs[idx].items():
                if val < self.atom_density_limit or nuc not in self.burn_nucs:
                    materials[int(burn_mat)-1].remove_nuclide(nuc)
                else:
                    materials[int(burn_mat)-1].remove_nuclide(nuc)
                    materials[int(burn_mat)-1].add_nuclide(nuc,val, 'ao')
                    atoms_gram_per_mol += val * atomic_mass(nuc)

            #ensure constant density is set and assign new volume
            density = materials[int(burn_mat)-1].get_mass_density()
            materials[int(burn_mat)-1].set_density('g/cm3', density)
            self.operator.number.volume[idx] = atoms_gram_per_mol /\
                                      AVOGADRO / density
        self.materials = materials

    def _build_parametric_model(self, param):
        """
        """
        model = deepcopy(self.model)
        self._set_geom_coeff(param)
        model.export_to_xml()
        return model

    def _finalize(self):
        init_geom()

    def msr_criticality_search(self, x):
        """
        """
        self._normalize_nuclides(x)

        lower_range = self.range[0]
        upper_range = self.range[1]

        if -1.0 < self.coeff < 1.0:
            self.tol /= 2
        else:
            self.tol /= abs(self.coeff)

        res = None
        while res == None:

            if self.refine_search and self.coeff >= abs(self.range[2])/2:
                check_brackets = True
            else:
                check_brackets = False

            search = search_for_keff(self._build_parametric_model,
                    bracket=[self.coeff + lower_range, self.coeff + upper_range],
                    tol=self.tol, bracketed_method= self.bracketed_method,
                    target=self.target, print_iterations=True)
                    #check_brackets=check_brackets)

            if len(search) == 3:
                _res, _, _ = search
                # Further check, in case upper limit gets hit
                if _res <= self.range[2]:
                    res = _res
                else:
                    if self.refuel:
                        msg = 'INFO: Hit upper limit {:.2f} cm Update geom' \
                        'coeff to {:.2f} and start refuel'.format(self.range[2],
                                                               self.init_param)
                        print(msg)
                        res = self.start_param
                        break
                    else:
                        msg = 'STOP: Hit upper limit and no further criteria' \
                               'defined'
                        raise Exception(msg)

            elif len(search) == 2:
                guesses, k = search

                if guesses[-1] > self.range[2]:
                    if self.refuel:
                        msg = 'INFO: Hit upper limit {:.2f} cm Update geom' \
                        'coeff to {:.2f} and start refuel'.format(self.range[2],
                                                               self.init_param)
                        print(msg)
                        res = self.start_param
                        break
                    else:
                        msg = 'STOP: Hit upper limit and no further criteria' \
                               'defined'
                        raise Exception(msg)

                if np.array(k).prod() < self.target:
                    print ('INFO: Function returned values BELOW target,' \
                           'adapting bracket range...')
                    if (self.target - np.array(k).prod()) <= 0.02:
                        lower_range = upper_range - 3
                    elif 0.02 < (self.target - np.array(k).prod()) < 0.03:
                        lower_range = upper_range - 1
                    else:
                        lower_range = upper_range - 0.5
                    upper_range += abs(self.range[1])/2

                else:
                    print ('INFO: Function returned values ABOVE target,' \
                           'adapting bracket range...')
                    upper_range = lower_range + 2
                    lower_range -= abs(self.range[0])

            else:
                raise ValueError(f'ERROR: search_for_keff output not valid')

        #probably not needed
        self._set_geom_coeff(res)
        self._finalize()
        msg = 'UPDATE: old coeff: {:.2f} cm --> ' \
              'new coeff: {:.2f} cm'.format(self.coeff, res)
        print(msg)
        diff = res - self.coeff
        return res, diff

class MsrBatchwiseMat(MsrBatchwise):
    """
    """
    def __init__(self, op, model, mat_id, refuel_vector, range=None,
                 bracketed_method='brentq', tol=0.01, target=1.0,
                 atom_density_limit=None):

        super().__init__(op, model, range, bracketed_method, tol, target,
                         atom_density_limit)
        if str(mat_id) not in self.burnable_mats:
                raise ValueError(f'Mat_id: {mat_id} is not a valid depletable material id')
        else:
            self.mat_id = mat_id
        self.refuel_vector = refuel_vector

    def initialize_geometry(self, geom_id=None, res=None):
        if res is not None:
            self.set_surface(geom_id, res)

    def _build_parametric_model(self, param):
        """
        """
        model = deepcopy(self.model)

        for idx, burn_id in enumerate(self.burnable_mats):
            for nuc, val in self.burn_nucs[idx].items():
                if nuc not in self.refuel_vector.keys():
                    if val < self.atom_density_limit or nuc not in self._burnable_nucs:
                        model.materials[int(burn_id)-1].remove_nuclide(nuc)
                    else:
                        model.materials[int(burn_id)-1].remove_nuclide(nuc)
                        model.materials[int(burn_id)-1].add_nuclide(nuc,val,'ao')
                else:
                    if model.materials[int(burn_id)-1].id == int(self.mat_id):
                        model.materials[int(burn_id)-1].remove_nuclide(nuc)
                        # convert grams into atoms
                        atoms = param/atomic_mass(nuc)*AVOGADRO*self.refuel_vector[nuc]
                        model.materials[int(burn_id)-1].add_nuclide(nuc,val+atoms,'ao')
                    else:
                        model.materials[int(burn_id)-1].remove_nuclide(nuc)
                        model.materials[int(burn_id)-1].add_nuclide(nuc,val,'ao')
            # ensure density is set
            density = model.materials[int(burn_id)-1].get_mass_density()
            model.materials[int(burn_id)-1].set_density('g/cm3',density)
        model.export_to_xml()
        return model

    def _update_x_vector_and_volumes(self, x, res):
        diff = dict()
        for idx, burn_id in enumerate(self.burnable_mats):
            #to calculate the new volume we need total atoms-gram per mol
            atoms_gram_per_mol = 0
            for id_nuc, nuc in enumerate(self.number.burnable_nuclides):
                if nuc in self.refuel_vector.keys() and int(burn_id) == int(self.mat_id):
                    # Convert res grams into atoms
                    res_atoms = res / atomic_mass(nuc) * AVOGADRO * self.refuel_vector[nuc]
                    diff[nuc] = x[idx][id_nuc] - res_atoms
                    x[idx][id_nuc] += res_atoms
                atoms_gram_per_mol += x[idx][id_nuc]*atomic_mass(nuc)
            # Calculate new volume and assign it in memory
            vol = atoms_gram_per_mol/AVOGADRO/self.materials[int(self.burnable_mats[idx])-1].get_mass_density()
            self.number.volume[idx] = vol
        return x, diff

    def perform_bw_search_for_keff(self, x):
        pass
