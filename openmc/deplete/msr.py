from abc import ABC, abstractmethod
from collections import OrderedDict
from copy import deepcopy
from warnings import warn
from numbers import Real

import numpy as np

from openmc.checkvalue import check_type, check_value, check_less_than, \
check_iterable_type, check_length
from openmc import Materials, Material
from openmc.search import _SCALAR_BRACKETED_METHODS, search_for_keff
from openmc.data import atomic_mass, AVOGADRO, ELEMENT_SYMBOL
import openmc.lib

class MsrContinuous:
    """Class defining Molten salt reactor (msr) elements (fission products)
    removal, based on removal rates and cycle time concepts.

    An instance of this class can be passed directly to an instance of the
    integrator class, such as :class:`openmc.deplete.CECMIntegrator`.
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

    def __init__(self, operator, model, units='1/s'):

        self.operator = operator
        self.materials = model.materials
        self.burn_mats = operator.burnable_mats

        #initialize removal rates container
        self.removal_rates = OrderedDict((mat, OrderedDict()) for mat in \
                                          self.burn_mats)
        self.index_transfer = set()
        self.units = units

    @property
    def units(self):
        return self._units

    @units.setter
    def units(self, value):
        check_value('Units', value, ['1/s', '1/h', '1/d'])
        self._units = value

    def _get_mat_id(self, val):
        """Helper method for getting material id from Material obj or name"""

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
        """Extract removal rates
        Parameters
        ----------
        mat : Openmc,Material or int
            Depletable material
        element : str
            Element to extract removal rate value
        Returns:
        ----------
        removal_rate : float
            Removal rate value
        """
        mat = self._get_mat_id(mat)
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
        destination_mat : str
            Depletable material id to where elements get transferred
        """
        mat = self._get_mat_id(mat)
        if element in self.removal_rates[mat]:
            return self.removal_rates[mat][element][1]

    def get_elements(self, mat):
        """Extract removing elements for a given material
        Parameters
        ----------
        mat : Openmc,Material or int
            Depletable material
        Returns:
        ----------
        elements : list
            List of elements
        """
        mat = self._get_mat_id(mat)
        if mat in self.removal_rates.keys():
            return self.removal_rates[mat].keys()

    def set_removal_rate(self, mat, elements, removal_rate, dest_mat=None):
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
        mat = self._get_mat_id(mat)

        check_type('removal_rate', removal_rate, Real)

        if dest_mat is not None:
            #prevent for setting tranfert to material if not set as depletable
            if len(self.burn_mats) > 1:
                check_value('transfert to material', str(dest_mat),
                            self.burn_mats)
            else:
                raise ValueError(f'Transfer to material {dest_mat} is set '\
                        'but there is only one depletable material')

            dest_mat = self._get_mat_id(dest_mat)

        for element in elements:
            check_value('Element', element, ELEMENT_SYMBOL.values())
            self.removal_rates[mat][element] = removal_rate, dest_mat
            if dest_mat is not None:
                self.index_transfer.add((dest_mat, mat))


class MsrBatchwise(ABC):
    """Abstract Base Class for implementing msr batchwise classes.

    Users should instantiate
    :class:`openmc.deplete.msr.MsrBatchwiseGeom` or
    :class:`openmc.deplete.msr.MsrBatchwiseMat` rather than this class.

    Parameters
    ----------
    operator : openmc.deplete.Operator
        OpenMC operator object
    model : openmc.model.Model
        OpenMC model object
    bracket : list of float
        Bracketing interval to search for the solution as list of float.
        This is equivalent to the `bracket` parameter of the `search_for_keff`.
    bracketed_method : {'brentq', 'brenth', 'ridder', 'bisect'}, optional
        Solution method to use; only applies if
        `bracket` is set, otherwise the Secant method is used.
        Defaults to 'brentq'.
    tol : float
        Tolerance for search_for_keff method
        Default to 0.01
    target : Real, optional
        keff value to search for, defaults to 1.0.

    Attributes
    ----------
    operator : openmc.deplete.Operator
        OpenMC operator object
    burn_mats : list of str
        List of burnable materials ids
    burn_nucs : list of str
        List of nuclides with available data for burnup
    model : openmc.model.Model
        OpenMC model object
        bracket : list of float
        list of floats defining bracketed bracket for `search_for_keff`
    bracketed_method : string
        Bracketed method for search_for_keff
    tol : float
        Tolerance for search_for_keff method
    target : float
        Search_for_keff function target
    """
    def __init__(self, operator, model, bracket, bracket_limit,
                 bracketed_method='brentq', tol=0.01, target=1.0,
                 print_iterations=True, search_for_keff_output=False):

        self.operator = operator
        self.burn_mats = operator.burnable_mats
        self.model = model

        check_iterable_type('bracket', bracket, Real)
        check_length('bracket', bracket, 2)
        check_less_than('bracket values', bracket[0], bracket[1])
        self.bracket = bracket

        check_iterable_type('bracket_limit', bracket_limit, Real)
        check_length('bracket_limit', bracket_limit, 2)
        check_less_than('bracket limit values',
                         bracket_limit[0], bracket_limit[1])
        self.bracket_limit = bracket_limit

        self.bracketed_method = bracketed_method
        self.tol = tol
        self.target = target
        self.print_iterations = print_iterations
        self.search_for_keff_output = search_for_keff_output

    @property
    def bracketed_method(self):
        return self._bracketed_method

    @bracketed_method.setter
    def bracketed_method(self, value):
        check_value('bracketed_method', value, _SCALAR_BRACKETED_METHODS)
        if value != 'brentq':
            warn('brentq bracketed method is recomended')
        self._bracketed_method = value

    @property
    def tol(self):
        return self._tol

    @tol.setter
    def tol(self, value):
        check_type("tol", value, Real)
        self._tol = value

    @property
    def target(self):
        return self._target

    @target.setter
    def target(self, value):
        check_type("target", value, Real)
        self._target = value

    @abstractmethod
    def _model_builder(self, param):
        """
        Builds the parametric model to be passed to `search_for_keff`.
        Callable function which builds a model according to a passed
        parameter. This function must return an openmc.model.Model object.
        Parameters
        ------------
        param : parameter
            model function variable
        Returns
        ------------
        _model :  openmc.model.Model
            OpenMC parametric model
        """

    def _conc_dict(self, x):
        """
        Order concentrations vector into a list of nuclides ordered
        dictionaries for each depletable material for easier extraction.
        Parameters
        ------------
        x : list of numpy.ndarray
            Total atoms to be used in function.
        Returns
        ------------
        ord_x : list of OrderedDict
            List of OrderedDict of nuclides for each depletable material
        """
        conc = list()
        for i, mat in enumerate(self.burn_mats):
            nucs = OrderedDict()
            for j, nuc in enumerate(self.operator.number.burnable_nuclides):
                nucs[nuc] = x[i][j]
            conc.append(nucs)
        return conc

    def _msr_search_for_keff(self, x, val):
        """
        Perform the criticality search for a given parametric model.
        It can either be a geometrical or material `search_for_keff`.
        Parameters
        ------------
        x : list of numpy.ndarray
            Total atoms to be used in function for normalization
        Returns
        ------------
        res : float
            Estimated value of the variable parameter where keff is the
            targeted value
        """
        _bracket = deepcopy(self.bracket)
        _tol = self.tol
        # Normalize search_for_keff tolerance with guess value
        if not -1.0 < val < 1.0:
            _tol = self.tol / abs(val)

        # Run until a search_for_keff root is found or ouf ot limits
        res = None
        while res == None:

            search = search_for_keff(self._model_builder,
                    bracket = [_bracket[0]+val, _bracket[1]+val],
                    tol = _tol,
                    bracketed_method = self.bracketed_method,
                    target = self.target,
                    print_iterations = self.print_iterations,
                    run_args = {'output': self.search_for_keff_output})

            if len(search) == 3:
                _res, _, _ = search
                if self.bracket_limit[0] < _res <  self.bracket_limit[1]:
                    res = _res
                else:
                    # Set res with closest limit and continue
                    arg_min = abs(np.array(res.bracket_limit)-_res).argmin()
                    warn('WARNING: Search_for_keff returned root out of '\
                         'bracket limit. Set root to {:.2f} and continue.'
                         .format(self.bracket_limit[arg_min]))
                    res = self.bracket_limit[arg_min]

            elif len(search) == 2:
                guess, k = search

                if self.bracket_limit[0] < np.all(guess) < self.bracket_limit[1]:
                    #Simple method to iteratively adapt the bracket
                    print('INFO: Function returned values below or above ' \
                           'target. Adapt bracket...')

                    # if the bracket ends up being smaller than the std of the
                    # keff's closer value to target, no need to continue-
                    if np.all(k) <= max(k).s:
                        res = max(k).n

                    # Calculate gradient as ratio of delta bracket and delta k
                    grad = abs(np.diff(_bracket) / np.diff(k))[0].n

                    # Move the bracket closer to presumed keff root
                    if np.mean(k) < self.target:
                        sign = np.sign(guess[np.argmax(k)])
                        _bracket[np.argmin(k)] = _bracket[np.argmax(k)]
                        _bracket[np.argmax(k)] += grad * (self.target - \
                                                  max(k).n) * sign
                    else:
                        sign = np.sign(guess[np.argmin(k)])
                        _bracket[np.argmax(k)] = _bracket[np.argmin(k)]
                        _bracket[np.argmin(k)] += grad * (min(k).n - \
                                                  self.target) * sign

                else:
                    # Set res with closest limit and continue
                    arg_min = abs(np.array(res.bracket_limit)-guess).argmin()
                    warn('WARNING: Adaptive iterative bracket went off '\
                         'bracket limits. Set root to {:.2f} and continue.'
                         .format(self.bracket_limit[arg_min]))
                    res = self.bracket_limit[arg_min]

            else:
                raise ValueError(f'ERROR: Search_for_keff output is not valid')

        return x, res

class MsrBatchwiseGeom(MsrBatchwise):
    """ MsrBatchwise geoemtrical class

    Instances of this class can be used to define geometry-based  criticality
    actions during a transport-depletion calculation.

    An instance of this class can be passed directly to an instance of the
    integrator class, such as :class:`openmc.deplete.CECMIntegrator`.

    Parameters
    ----------
    operator : openmc.deplete.Operator
        OpenMC operator object
    model : openmc.model.Model
        OpenMC model object
    name_or_id : int
        Parametric goemetrical feature id from the geometry model
    obj_type : string, optional
        Parametric goemetrical feature obj_type. Right now only 'cell' is
        supported
        Default to 'cell'
    bracket : list of float
        List of floats defining bracketed bracket for search_for_keff
    bracketed_method : string
        Bracketed method for search_for_keff
        Default to 'brentq' --> more robuts
    tol : float
        Tolerance for search_for_keff method
        Default to 0.01
    target : float
        Search_for_keff function target
        Default to 1.0

    Attributes
    ----------
    name_or_id : int
        Parametric goemetrical feature id from the geometry model
    obj_type : string
        Parametric goemetrical feature obj_type. Right now only 'cell' is
        supported
    """
    def __init__(self, operator, model, name_or_id, axis, bracket,
                 bracket_limit, print_iterations=True, bracketed_method='brentq',
                 tol=0.01, target=1.0):

        super().__init__(operator, model, bracket, bracket_limit,
                         bracketed_method, tol, target, print_iterations)

        if isinstance(name_or_id, int):
            check_value('id', name_or_id, [cell.id for cell in \
                                self.model.geometry.get_all_cells().values()])
        elif isinstance(name_or_id, str):
            check_value('name', name_or_id, [cell.name for cell in \
                                self.model.geometry.get_all_cells().values()])
        else:
            ValueError(f'{name_or_id} is neither a str nor an int')
        self.name_or_id = name_or_id

        #index of vector axis
        if axis == 'x':
            self.axis = 0
        elif axis == 'y':
            self.axis = 1
        elif axis == 'z':
            self.axis = 2
        else:
            raise ValueError(f'{axis} value to recognized, should be x, y or z')

        # Initialize translation vector
        self.vector = np.zeros(3)

    def _get_cell_attrib(self):
        """
        Extract cell translation array from cell id or cell name and return
         The translation is only applied to cells filled with a universe
        Returns
        ------------
        coeff : float
            cell coefficient
        """
        for cell in openmc.lib.cells.values():
            if cell.id == self.name_or_id or cell.name == self.name_or_id:
                return cell.translation[self.axis]

    def _set_cell_attrib(self, val, attrib_name='translation'):
        """
        Set translation coeff to the cell.
        The translation is only applied to cells filled with a universe
        Parameters
        ------------
        var : float
            Surface coefficient to set
        geometry : openmc.model.geometry
            OpenMC geometry model
        """
        self.vector[self.axis] = val
        for cell in openmc.lib.cells.values():
            if cell.id == self.name_or_id or cell.name == self.name_or_id:
                setattr(cell, attrib_name, self.vector)

    def _update_materials(self, x):
        """
        Export a meterial xml to be red by the search_for_keff function
        based on available crosse section nuclides with data and
        `atom_density_limit` argument, if defined. This function also calculates
        the new total volume in [cc] from nuclides atom densities, assuming
        constant material density as set by the user in the initial material
        definition.
        The new volume is then passed to the `openmc.deplete.AtomNumber` class
        to renormalize the atoms vector in the model in memory before running
        the next transport iteration.
        Parameters
        ------------
        x : list of numpy.ndarray
            Total atoms to be used in function.
        """
        nucs = super()._conc_dict(x)

        for i, mat in enumerate(self.burn_mats):
            nuclides = []
            densities = []
            atoms_gram_per_mol = 0
            for nuc, val in nucs[i].items():
                if nuc in self.operator.nuclides_with_data:
                    if val > 0.0:
                        nuclides.append(nuc)
                        density = 1.0e-24 * val / self.operator.number.volume[i]
                        densities.append(density)
                atoms_gram_per_mol += val * atomic_mass(nuc)

            openmc.lib.materials[int(mat)].set_densities(nuclides, densities)

            #Assign new volume to AtomNumber
            mass_dens = [m.get_mass_density() for m in self.model.materials if
                    m.id == int(mat)][0]
            self.operator.number.volume[i] = atoms_gram_per_mol / AVOGADRO /\
                                                mass_dens

    def _model_builder(self, param):
        """
        Builds the parametric model to be passed to `msr_search_for_keff`
        Parameters
        ------------
        param :
            model function variable: geometrical coefficient
        Returns
        ------------
        _model :  openmc.model.Model
            OpenMC parametric model
        """
        self._set_cell_attrib(param)
        return self.model

    def msr_search_for_keff(self, x):
        """
        Perform the criticality search on the parametric geometrical model.
        Parameters
        ------------
        x : list of numpy.ndarray
            Total atoms to be used in function for normalization
        Returns
        ------------
        res : float
            geometrical coefficient root of search_for_keff function
        """
        self._update_materials(x)
        val = self._get_cell_attrib()
        x, res = super()._msr_search_for_keff(x, val)

        # set results value in the geometry model abd continue
        self._set_cell_attrib(res)
        print('UPDATE: old value: {:.2f} cm --> ' \
              'new value: {:.2f} cm'.format(val, res))
        return x

class MsrBatchwiseMat(MsrBatchwise):
    """ MsrBatchwise material class

    Instances of this class can be used to define material-based criticality
    actions during a transport-depletion calculation.

    An instance of this class can be passed directly to an instance of the
    integrator class, such as :class:`openmc.deplete.CECMIntegrator`.

    Parameters
    ----------
    operator : openmc.deplete.Operator
        OpenMC operator object
    model : openmc.model.Model
        OpenMC model object
    mat_id : int
        Material id to perform actions on
    refuel_vector : dict
        Refueling material nuclides composition in form of dict, where keys are
        the nuclides and values the fractions
    bracket : list of float
        List of floats defining bracketed bracket for search_for_keff
    bracketed_method : string
        Bracketed method for search_for_keff
        Default to 'brentq' --> more robuts
    tol : float
        Tolerance for search_for_keff method
        Default to 0.01
    target : float
        Search_for_keff function target
        Default to 1.0

    Attributes
    ----------
    mat_id : int
        Material id to perform actions on
    refuel_vector : dict
        Refueling material nuclides composition in form of dict, where keys are
        the nuclides str and values are the composition fractions.
    """
    def __init__(self, operator, model, mat_id, refuel_vector, bracket,
                 bracket_limit, print_iterations=True, bracketed_method='brentq',
                 tol=0.01, target=1.0):

        super().__init__(operator, model, bracket, bracket_limit,
                         bracketed_method, tol, target, print_iterations)

        check_value("mat id", str(mat_id), self.burn_mats)
        self.mat_id = mat_id

        check_type("refuel vector", refuel_vector, dict, str)
        for nuc in refuel_vector.keys():
            check_value("refuel vector nuclide", nuc,
                        self.operator.nuclides_with_data)
        self.refuel_vector = refuel_vector

    def _model_builder(self, param):
        """
        Builds the parametric model to be passed to `msr_search_for_keff`
        Parameters
        ------------
        param :
            Model function variable, mass of material to refuel in grams
        Returns
        ------------
        _model :  openmc.model.Model
            OpenMC parametric model
        """

        for i, mat in enumerate(self.burn_mats):
            nuclides = []
            densities = []
            for nuc, val in self.nucs[i].items():
                if nuc not in self.refuel_vector.keys():
                    if nuc in self.operator.nuclides_with_data and val > 0.0:
                        nuclides.append(nuc)
                        density = 1.0e-24 * val / self.operator.number.volume[i]
                        densities.append(density)

                else:
                    nuclides.append(nuc)
                    if openmc.lib.materials[int(mat)].id == self.mat_id:
                        # convert grams into atoms-b/cm3
                        density = 1.0e-24 * (val + (param / atomic_mass(nuc) * \
                                AVOGADRO * self.refuel_vector[nuc])) / \
                                self.operator.number.volume[i]
                    else:
                        density = 1.0e-24 * val / self.operator.number.volume[i]
                    densities.append(density)

            openmc.lib.materials[int(mat)].set_densities(nuclides, densities)
        return self.model

    def _update_x_vector_and_volumes(self, x, res):
        """
        Updates and returns total atoms with the root from the search_for_keff
        function. This function also calculates the new total volume
        in cc from nuclides atom densities, assuming constant material density
        as set by the user in the initial material definition.
        The new volume is then passed to the `openmc.deplete.AtomNumber` class
        to renormalize the atoms vector in the model in memory before running
        the next transport.
        Parameters
        ------------
        x : list of numpy.ndarray
            Total atoms to be used in function
        res : float
            Root of the search_for_keff function
        Returns
        ------------
        x : list of numpy.ndarray
            Updated total atoms to be used in function
        """

        for i, mat in enumerate(self.burn_mats):
            atoms_gram_per_mol = 0
            for j, nuc in enumerate(self.operator.number.burnable_nuclides):
                if nuc in self.refuel_vector.keys() and int(mat) == self.mat_id:
                    # Convert res grams into atoms
                    res_atoms = res / atomic_mass(nuc) * AVOGADRO * \
                                self.refuel_vector[nuc]
                    x[i][j] += res_atoms
                atoms_gram_per_mol += x[i][j] * atomic_mass(nuc)

            # Calculate new volume and assign it in memory
            mass_dens = [m.get_mass_density() for m in self.model.materials if
                    m.id == int(mat)][0]
            self.operator.number.volume[i] = atoms_gram_per_mol / AVOGADRO / \
                                             mass_dens

        return x

    def msr_search_for_keff(self, x):
        """
        Perform the criticality search on the parametric material model.
        Parameters
        ------------
        x : list of numpy.ndarray
            Total atoms to be used in function for normalization
        Returns
        ------------
        res : float
            Material mass root of search_for_keff function
        diff : float
            Difference from previous result
        """
        self.nucs = super()._conc_dict(x)
        x, res = super()._msr_search_for_keff(x, 0)

        print('UPDATE: material addition --> {:.2f} g --> '.format(res))
        x = self._update_x_vector_and_volumes(x, res)
        return  x
