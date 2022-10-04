from collections import OrderedDict
from openmc import Materials, Material, search_for_keff
from copy import deepcopy
from .openmc_operator import OpenMCOperator

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

    def _get_mat_index(self, mat):
        """Helper method for getting material index"""
        if isinstance(mat, Material):
            mat = str(mat.id)
        return self.index_burn[mat] if isinstance(mat, str) else mat

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


class MsrBatchwise(OpenMCOperator):
    """
    """
    def __init__(self, model, type='surface', geom_id, range,
                    keff_search_bracket_method='brentq', keff_search_tolerance,
                    keff_search_target=1.0, atom_density_limit=None,
                    refine_search=True, refuel=False, start_param=0.0):

        self.model = model
        self.settings = model.settings
        self.geometry = self.model.geometry
        self.materials = self.model.materials
        self.geom_id = geom_id

        if type == 'surface':
            for surf in self.geometry.get_all_surfaces().items():
                if surf[1].id == self.geom_id:
                    keys = list(surf[1].coefficients.keys())
                    if len(keys) == 1:
                        self.geom = getattr(surf[1], keys[0])
                    else:
                        msg=(f'Surface coefficients {keys} are not one')
                        raise Exception(msg)
        #TODO implement a universe type
        # elif type == 'universe'

        if len(range) != 3:
            raise ValueError(f'{range} lenght must be 3')
        else:
            self.range = range

        self.keff_search_bracket_method = keff_search_bracket_method
        self.keff_search_tolerance = keff_search_tolerance
        self.keff_search_target = keff_search_target
        self.atom_density_limit = atom_density_limit
        self.refine_search = refine_search
        self.refuel = refuel
        self.start_param = start_param

        super().__init__()

    def _order_burn_nuclides(self, x):
        """
        x : atoms vector
        """
        burn_nucs = list()
        for idx,burn_id in enumerate(self.burnable_mats):
            nucs = OrderedDict()
            for id_nuc, nuc in enumerate(self.number.burnable_nuclides):
                nucs[nuc] = x[idx][id_nuc]
            burn_nucs.append(nucs)
        return burn_nucs

    def _build_parametric_model(self, param):
        """
        """
        model = deepcopy(self.model)
        for surf in model.geometry.get_all_surfaces().items():
            if surf[1].id == surf_id:
                keys = list(surf[1].coefficients.keys())
                if len(keys) == 1:
                    setattr(surf[1],keys[0], param)
                else:
                    raise Exception(f'Surface coefficients {keys} are not one')
        model.export_to_xml()
        return model

    def _reinitilize_materials(self, x):
        """
        """
        materials = deepcopy(self.model.materials)
        atoms_vec = self._order_burn_nuclides(x)

        for idx, burn_id in enumerate(self.burnable_mats):
            atoms_gram_per_mol = 0
            for nuc,val in atoms_vec[idx].items():
                if val < density_limit or nuc not in self._burnable_nucs:
                    materials[int(burn_id)-1].remove_nuclide(nuc)
                else:
                    materials[int(burn_id)-1].remove_nuclide(nuc)
                    materials[int(burn_id)-1].add_nuclide(nuc,val,'ao')
                    atoms_gram_per_mol += val*openmc.data.atomic_mass(nuc)

            #ensure constant density is set and assign new volume
            density = materials[int(burn_id)-1].get_mass_density()
            materials[int(burn_id)-1].set_density('g/cm3',density)
            self.number.volume[idx] = atoms_gram_per_mol/openmc.data.AVOGADRO/density
        self.materials = materials

    def _update_critical_param(self, res):
        """
        """
        # Set paramteric result in the geometry model (really needed?)
        for surf in self.geometry.get_all_surfaces().items():
            if surf[1].id == self.geom_id:
                keys = list(surf[1].coefficients.keys())
                if len(keys) == 1:
                    setattr(surf[1], keys[0], res)
                else:
                    msg = (f'Surface coefficients {keys} are not one')
                    raise Exception(msg)

class MsrBatchwiseGeom(MsrBatchwise):
    """
    """
    def __init__(self, model, type='surface', geom_id, range,
                    keff_search_bracket_method='brentq', keff_search_tolerance,
                    keff_search_target=1.0, atom_density_limit=None,
                    refine_search=True, refuel=False, start_param=0.0):

        super().__init__(model, type, geom_id, range, keff_search_bracket_method,
                         keff_search_tolerance, keff_search_target, atom_density_limit,
                         refine_search, refuel, start_param)

    def perform_bw_search_for_keff(self, x):
        """
        """
        self._reinitilize_materials(x)

        res = None
        lower_range = self.range[0]
        upper_range = self.range[1]
        if -1.0 < self.geom < 1.0:
            self.keff_search_tolerance /= 2
        else:
            self.keff_search_tolerance /= abs(self.geom)

        while res == None:
            if self.refine_search and self.geom >= abs(self.range[2])/2:
                check_brackets = True
            else:
                check_brackets = False

            search = search_for_keff(self._build_parametric_model,
                        bracket=[self.geom + lower_range,self.geom+upper_range],
                        tol=self.keff_search_tolerance,
                        bracketed_method= self.keff_search_bracket_method,
                        target=self.keff_search_target, print_iterations=True,
                        check_brackets = check_brackets)

            if len(search) == 3:
                _res, guesses, k = search
                # Further check, in case upper limit gets hit
                if _res <= self.range[2]:
                    res = _res
                else:
                    if self.refuel:
                        print(f'Upper limit reached, level brought down to {init_param}cm and \
                        start refueling...')
                        res = self.start_param
                        break
                    else:
                        raise Exception(f'Upper limit reached, stopping depletion')

            elif len(search) == 2:
                guesses, k = search

                if guesses[-1] > self.range[2]:
                    if self.refuel:
                        print(f'Upper limit reached, level brought down to {init_param}cm and \
                        start refueling...')
                        res = self.start_param
                        break
                    else:
                        raise Exception(f'Upper limit reached, stopping depletion')

                if np.array(k).prod() < self.keff_search_target:
                    print ("Invalid range, below target")
                    if (self.keff_search_target - np.array(k).prod()) <= 0.02:
                        lower_range = upper_range - 3
                    elif 0.02 < (self.keff_search_target - np.array(k).prod()) < 0.03:
                        lower_range = upper_range - 1
                    else:
                        lower_range = upper_range - 0.5
                    upper_range += abs(self.range[1])/2

                else:
                    print ("Invalid range, above target")
                    upper_range = lower_range + 2
                    lower_range -= abs(self.range[0])

            else:
                raise Exception(f'search_for_keff output not contemplated')

        #probably not needed
        self._update_critical_param(res)

        print(f'res: {res}, guess: {guess}')
        diff = res - guess

        return res

class MsrBatchwiseMat(MsrBatchwise):
    """
    """
    def __init__(self, model, refuel_vector, range, keff_search_bracket_method='brentq',
                 keff_search_tolerance, keff_search_target=1.0,
                atom_density_limit=None):

        super().__init__(model, range, keff_search_bracket_method,
                         keff_search_tolerance, keff_search_target, atom_density_limit)
                         
        self.refuel_vector = refuel_vector

    def perform_bw_search_for_keff(self, res, x):

        # In case of hitting top limit, refuel
        if res == sel.start_param:
            refuel = self.refuel
            x, diff_mat = keff_control_material(x, refuel['mat_id'], refuel['range'],
                        refuel['bracketed_method'], refuel['mat_comp'],
                        refuel['tol'], target, density_limit, surf_id, res)
        return x, diff
