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
        Returns:
        ----------
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
            #prevent for setting tranfert to material if not set as depletable
            if len(self.burn_mats) > 1:
                check_value('transfert to material', str(dest_mat),
                            self.burn_mats)
            else:
                raise ValueError(f'Transfer to material {dest_mat} is set '\
                        'but there is only one depletable material')
            dest_mat = self._get_mat_id(dest_mat)
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


class MsrBatchwise(ABC):
    """Abstract Base Class for implementing msr batchwise classes.

    Users should instantiate:
    :class:`openmc.deplete.msr.MsrBatchwiseGeom` or
    :class:`openmc.deplete.msr.MsrBatchwiseMat` rather than this class.

    Parameters
    ----------
    operator : openmc.deplete.Operator
        OpenMC operator object
    model : openmc.model.Model
        OpenMC model object
    bracket : list of float
        Bracketing range around the guess value to search for the solution as
        list of float.
        This is equivalent to the `bracket` parameter of the `search_for_keff`.
    bracket_limit : list of float
        Upper and lower limits for the search_for_keff. If search_for_keff root
        or guesses fall above the range, the closest limit will be taken and
        set as new result.
    bracketed_method : {'brentq', 'brenth', 'ridder', 'bisect'}, optional
        Solution method to use.
        This is equivalent to the `bracket_method` parameter of the
        `search_for_keff`.
        Defaults to 'brentq'.
    tol : float
        Tolerance for search_for_keff method.
        This is equivalent to the `tol` parameter of the `search_for_keff`.
        Default to 0.01
    target : Real, optional
        This is equivalent to the `target` parameter of the `search_for_keff`.
        Default to 1.0.
    print_iterations : Bool, Optional
        Wheter or not to print `search_for_keff` iterations.
        Default to True
    search_for_keff_output : Bool, Optional
        Wheter or not to print transport iterations during  `search_for_keff`.
        Default to False
    Attributes
    ----------
    operator : openmc.deplete.Operator
        OpenMC operator object
    burn_mats : list of str
        List of burnable materials ids
    model : openmc.model.Model
        OpenMC model object
    bracket : list of float
        Bracketing range around the guess value to search for the solution as
        list of float.
        This is equivalent to the `bracket` parameter of the `search_for_keff`.
    bracket_limit : list of float
        Upper and lower limits for the search_for_keff. If search_for_keff root
        or guesses fall above the range, the closest limit will be taken and
        set as new result.
    bracketed_method : {'brentq', 'brenth', 'ridder', 'bisect'}, optional
        Solution method to use.
        This is equivalent to the `bracket_method` parameter of the
        `search_for_keff`.
    tol : float
        Tolerance for search_for_keff method.
        This is equivalent to the `tol` parameter of the `search_for_keff`.
    target : Real, optional
        This is equivalent to the `target` parameter of the `search_for_keff`.
    print_iterations : Bool, Optional
        Wheter or not to print `search_for_keff` iterations.
    search_for_keff_output : Bool, Optional
        Wheter or not to print transport iterations during  `search_for_keff`.
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
        if bracket_limit[0] > bracket[0]:
            raise ValueError('Lower bracket limit {} is greater than lower '
                             'bracket {}'.format(bracket_limit[0], bracket[0]))
        if bracket_limit[1] < bracket[0]:
            raise ValueError('Upper bracket limit {} is less than upper '
                             'bracket {}'.format(bracket_limit[1], bracket[1]))
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
            warn('brentq bracketed method is recomended here')
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

    def _msr_search_for_keff(self, x, val):
        """
        Perform the criticality search for a given parametric model.
        It supports geometrical or material based `search_for_keff`.
        Parameters
        ------------
        x : list of numpy.ndarray
            Atoms concentration vector
        Returns
        ------------
        x : list of numpy.ndarray
            Modified atoms concentration vector
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

            # if len(search) is 3 search_for_keff was successful
            if len(search) == 3:
                _res, _, _ = search
                #Check if root is within bracket limits
                if self.bracket_limit[0] < _res <  self.bracket_limit[1]:
                    res = _res
                else:
                    # Set res with the closest limit and continue
                    arg_min = abs(np.array(self.bracket_limit)-_res).argmin()
                    warn('WARNING: Search_for_keff returned root out of '\
                         'bracket limit. Set root to {:.2f} and continue.'
                         .format(self.bracket_limit[arg_min]))
                    res = self.bracket_limit[arg_min]

            elif len(search) == 2:
                guess, k = search
                #Check if all guesses are within bracket limits
                if all(self.bracket_limit[0] < g < self.bracket_limit[1] \
                    for g in guess):
                    #Simple method to iteratively adapt the bracket
                    print('INFO: Function returned values below or above ' \
                           'target. Adapt bracket...')

                    # if the bracket ends up being smaller than the std of the
                    # keff's closer value to target, no need to continue-
                    if all(_k <= max(k).s for _k in k):
                        arg_min = abs(self.target-np.array(guess)).argmin()
                        res = guess[arg_min]

                    # Calculate gradient as ratio of delta bracket and delta k
                    grad = abs(np.diff(_bracket) / np.diff(k))[0].n
                    # Move the bracket closer to presumed keff root.
                    # 2 cases: both k are below or above target
                    if np.mean(k) < self.target:
                        # direction of moving bracket: +1 is up, -1 is down
                        if guess[np.argmax(k)] > guess[np.argmin(k)]:
                            dir = 1
                        else:
                            dir = -1
                        _bracket[np.argmin(k)] = _bracket[np.argmax(k)]
                        _bracket[np.argmax(k)] += grad * (self.target - \
                                                  max(k).n) * dir
                    else:
                        if guess[np.argmax(k)] > guess[np.argmin(k)]:
                            dir = -1
                        else:
                            dir = 1
                        _bracket[np.argmax(k)] = _bracket[np.argmin(k)]
                        _bracket[np.argmin(k)] += grad * (min(k).n - \
                                                  self.target) * dir

                else:
                    # Set res with closest limit and continue
                    arg_min = abs(np.array(self.bracket_limit)-guess).argmin()
                    warn('WARNING: Adaptive iterative bracket went off '\
                         'bracket limits. Set root to {:.2f} and continue.'
                         .format(self.bracket_limit[arg_min]))
                    res = self.bracket_limit[arg_min]

            else:
                raise ValueError(f'ERROR: Search_for_keff output is not valid')

        return x, res

    def _save_res(self, step_index, res):
        """
        Save results to msr_results.h5 file.
        Parameters
        ------------
        step_index : int
            depletion time step index
        res : float or dict
             Root of the search_for_keff function
        """
        kwargs = {'mode': "w" if step_index == 0 else "a"}
        with h5py.File('msr_results.h5', **kwargs) as h5:
            h5.create_dataset(str(step_index), data=res)


class MsrBatchwiseGeom(MsrBatchwise):
    """ MsrBatchwise geoemtrical class

    Instances of this class can be used to define geometrical based criticality
    actions during a transport-depletion calculation.
    Currently only translation is supported.
    In this case a geoemtrical cell translation coefficient will be used as
    parametric variable.
    The user should remember to fill the cell with a Universe.

    An instance of this class can be passed directly to an instance of the
    integrator class, such as :class:`openmc.deplete.CECMIntegrator`.

    Parameters
    ----------
    operator : openmc.deplete.Operator
        OpenMC operator object
    model : openmc.model.Model
        OpenMC model object
    cell_id_or_name : Openmc.Cell or int or str
        Identificative of parametric cell
    axis : int
        cell translation direction axis, where 0 is 'x', 1 is 'y' and 2 'z'.
    bracket : list of float
        Bracketing range around the guess value to search for the solution as
        list of float in cm.
        In this case the guess guess value is the translation coefficient result
        of the previous depletion step
    bracket_limit : list of float
        Upper and lower limits in cm. If search_for_keff root
        or guesses fall above the range, the closest limit will be taken and
        set as new result. In this case the limit should coincide with the
        cell geometrical boundary conditions.
    bracketed_method : {'brentq', 'brenth', 'ridder', 'bisect'}, optional
        Solution method to use.
        This is equivalent to the `bracket_method` parameter of the
        `search_for_keff`.
        Default to 'brentq'
    tol : float
        Tolerance for search_for_keff method.
        This is equivalent to the `tol` parameter of the `search_for_keff`.
        Default to 0.01
    target : Real, optional
        This is equivalent to the `target` parameter of the `search_for_keff`.
        Default to 1.0
    Attributes
    ----------
    cell_id : openmc.Cell or int or str
        Identificative of parametric cell
    axis : int
        cell translation direction axis, where 0 is 'x', 1 is 'y' and 2 'z'.
    vector : numpy.array
        translation vector
    """
    def __init__(self, operator, model, cell_id_or_name, axis, bracket,
                 bracket_limit, print_iterations=True, bracketed_method='brentq',
                 tol=0.01, target=1.0):

        super().__init__(operator, model, bracket, bracket_limit,
                         bracketed_method, tol, target, print_iterations)

        self.cell_id = self._get_cell_id(cell_id_or_name)

        #index of cell translation direction axis
        check_value('axis', axis, [0,1,2])
        self.axis = axis

        # Initialize translation vector
        self.vector = np.zeros(3)

    def _get_cell_id(self, val):
        """Helper method for getting cell id from Cell obj or name.
        Parameters
        ----------
        val : Openmc.Cell or str or int representing Cell
        Returns
        ----------
        id : str
            Cell id
        """
        if isinstance(val, Cell):
            check_value('Cell id', val.id, [cell.id for cell in \
                                self.model.geometry.get_all_cells().values()])
            val = val.id

        elif isinstance(val, str):
            if val.isnumeric():
                check_value('Cell id', val, [str(cell.id) for cell in \
                                self.model.geometry.get_all_cells().values()])
                val = int(val)
            else:
                check_value('Cell name', val, [cell.name for cell in \
                                self.model.geometry.get_all_cells().values()])

                val = [cell.id for cell in \
                    self.model.geometry.get_all_cells().values() \
                    if cell.name == val][0]

        elif isinstance(val, int):
            check_value('Cell id', val, [cell.id for cell in \
                                self.model.geometry.get_all_cells().values()])

        else:
            ValueError(f'Cell: {val} is not recognized')

        return val

    def _get_cell_attrib(self):
        """
        Get cell translation coefficient.
        The translation is only applied to cells filled with a universe
        Returns
        ------------
        coeff : float
            cell coefficient
        """
        for cell in openmc.lib.cells.values():
            if cell.id == self.cell_id:
                return cell.translation[self.axis]

    def _set_cell_attrib(self, val, attrib_name='translation'):
        """
        Set translation coeff to the cell in memeory.
        The translation is only applied to cells filled with a universe
        Parameters
        ------------
        var : float
            Surface coefficient to set
        geometry : openmc.model.geometry
            OpenMC geometry model
        attrib_name : str
            Currently only translation is implemented
            Default to 'translation'
        """
        self.vector[self.axis] = val
        for cell in openmc.lib.cells.values():
            if cell.id == self.cell_id or cell.name == self.cell_id:
                setattr(cell, attrib_name, self.vector)

    def _update_materials(self, x):
        """
        Updates nuclide densities concentration vector coming from
        Bateman quations solution and assign to the model in memory.
        By doing so it also calculates the total number of atoms-grams per mol.
        This quantity divided by the total mass density gives a volume
        in cc. This can be inteded as the new material volume if we were to fix
        the total material mass density. If not set, the volume will remain
        constant and the density will update at the next depletion step.
        Parameters
        ------------
        x : list of numpy.ndarray
            Total atom concentrations
        """
        self.operator.number.set_density(x)
        for i, mat in enumerate(self.burn_mats):
            nuclides = []
            densities = []
            density = 0
            for nuc in self.operator.number.nuclides:
                # get nuclide density [atoms/b-cm]
                val = 1.0e-24 * self.operator.number.get_atom_density(mat, nuc)
                if nuc in self.operator.nuclides_with_data:
                    if val > 0.0:
                        nuclides.append(nuc)
                        densities.append(val)
                # density in [atoms-g/b-cm-mol]
                density +=  val * atomic_mass(nuc)

            openmc.lib.materials[int(mat)].set_densities(nuclides, densities)

            # Mass density in [g/cc] that will be kept constant
            mass_dens = [m.get_mass_density() for m in self.model.materials if
                    m.id == int(mat)][0]
            #In the internal version we assign new volume to AtomNumber
            self.operator.number.volume[i] *= 1.0e24 * density / AVOGADRO /\
                                                mass_dens

    def _model_builder(self, param):
        """
        Builds the parametric model that is passed to the `msr_search_for_keff`
        function by setting the parametric variable to the geoemetry cell.
        Parameters
        ------------
        param : model parametricl variable
            cell translation coefficient
        Returns
        ------------
        _model :  openmc.model.Model
            OpenMC parametric model
        """
        self._set_cell_attrib(param)
        return self.model

    def msr_search_for_keff(self, x, step_index):
        """
        Perform the criticality search on the parametric geometrical model.
        Will set the root of the `search_for_keff` function to the cell
        attribute.
        Parameters
        ------------
        x : list of numpy.ndarray
            Total atoms concentrations
        Returns
        ------------
        x : list of numpy.ndarray
            Updated total atoms concentrations
        """
        val = self._get_cell_attrib()
        check_type('Cell coeff', val, Real)

        if step_index > 0:
            self._update_materials(x)
            x, res = super()._msr_search_for_keff(x, val)

            # set results value in the geometry model abd continue
            self._set_cell_attrib(res)
            print('UPDATE: old value: {:.2f} cm --> ' \
                  'new value: {:.2f} cm'.format(val, res))

            #Store results
            super()._save_res(step_index, res)
        else:
            super()._save_res(step_index, val)
        return x

class MsrBatchwiseMat(MsrBatchwise):
    """ MsrBatchwise material class

    Instances of this class can be used to define material based criticality
    actions during a transport-depletion calculation.

    An instance of this class can be passed directly to an instance of the
    integrator class, such as :class:`openmc.deplete.CECMIntegrator`.

    Parameters
    ----------
    operator : openmc.deplete.Operator
        OpenMC operator object
    model : openmc.model.Model
        OpenMC model object
    mat_id_or_name : openmc.Material or int or str
        Identificative of parametric material
    refuel_vector : dict
        Refueling material nuclides composition in form of dict, where keys are
        nuclides and values fractions.
        E.g., refuel_vector = {'U235':0.3,'U238':0.7}
    bracket : list of float
        Bracketing range quantity of material to add in grams.
    bracket_limit : list of float
        Upper and lower limits of material to add in grams.
    bracketed_method : {'brentq', 'brenth', 'ridder', 'bisect'}, optional
        Solution method to use.
        This is equivalent to the `bracket_method` parameter of the
        `search_for_keff`.
        Default to 'brentq'
    tol : float
        Tolerance for search_for_keff method.
        This is equivalent to the `tol` parameter of the `search_for_keff`.
        Default to 0.01
    target : Real, optional
        This is equivalent to the `target` parameter of the `search_for_keff`.
        Default to 1.0
    Attributes
    ----------
    mat_id : openmc.Material or int or str
        Identificative of parametric material
    refuel_vector : dict
        Refueling material nuclides composition in form of dict, where keys are
        the nuclides str and values are the composition fractions.
    """
    def __init__(self, operator, model, mat_id_or_name, refuel_vector, bracket,
                 bracket_limit, print_iterations=True, bracketed_method='brentq',
                 tol=0.01, target=1.0):

        super().__init__(operator, model, bracket, bracket_limit,
                         bracketed_method, tol, target, print_iterations)

        self.mat_id = self._get_mat_id(mat_id_or_name)

        check_type("refuel vector", refuel_vector, dict, str)
        for nuc in refuel_vector.keys():
            check_value("check nuclide exists", nuc,
                        self.operator.nuclides_with_data)
        if sum(refuel_vector.values()) != 1.0:
            raise ValueError('Refuel vector fractions {} do not sum up to 1.0'
                             .format(refuel_vector.values()))
        self.refuel_vector = refuel_vector

    def _check_nuclides(self, nucs):
        """Checks if refueling nuclides exist in the parametric material.
        Parameters
        ----------
        nucs : list of str
            Refueling nuclides
        """
        for nuc in nucs:
            check_value("check nuclide to refuel exists in mat", nuc,
                [mat.nuclides for mat_id, mat in openmc.lib.materials.items() \
                if mat_id == self.mat_id][0])

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
            check_value('Material id', str(val.id), self.burn_mats)
            val = val.id

        elif isinstance(val, str):
            if val.isnumeric():
                check_value('Material id', val, self.burn_mats)
                val = int(val)
            else:
                check_value('Material name', val,
                   [mat.name for mat in self.model.materials if mat.depletable])
                val = [mat.id for mat in self.model.materials \
                        if mat.name == val][0]

        elif isinstance(val, int):
            check_value('Material id', str(val), self.burn_mats)

        return val

    def _model_builder(self, param):
        """
        Builds the parametric model that is passed to the `msr_search_for_keff`
        function by updating the material densities and setting the parametric
        variable to the material nuclides to add.
        Parameters
        ------------
        param :
            Model function variable, quantity of material to refuel in grams
        Returns
        ------------
        _model :  openmc.model.Model
            OpenMC parametric model
        """

        for i, mat in enumerate(self.burn_mats):
            nuclides = []
            densities = []
            for nuc in self.operator.number.nuclides:
                if nuc not in self.refuel_vector.keys():
                    if nuc in self.operator.nuclides_with_data:
                        # get atoms density [atoms/b-cm]
                        val = 1.0e-24 * \
                              self.operator.number.get_atom_density(mat, nuc)
                        if val > 0.0:
                            nuclides.append(nuc)
                            densities.append(val)

                else:
                    nuclides.append(nuc)
                    val = 1.0e-24 * self.operator.number.get_atom_density(mat,
                                                                         nuc)
                    if int(mat) == self.mat_id:
                        # convert params [grams] into [atoms/cc]
                        param *= 1.0e-24 / atomic_mass(nuc) * AVOGADRO * \
                               self.refuel_vector[nuc] / \
                               self.operator.number.volume[i]
                    densities.append(val + param)

            openmc.lib.materials[int(mat)].set_densities(nuclides, densities)
        return self.model

    def _update_x_vector_and_volumes(self, x, res):
        """
        Updates and returns the total atoms concentrations vector with the root
        from the `search_for_keff`. It also calculates the new total volume
        in cc from the nuclides atom densities, assuming constant material mass
        density. The volume is then passed to the `openmc.deplete.AtomNumber`
        class to renormalize the atoms vector in the model in memory before
        running the next transport.
        Parameters
        ------------
        x : list of numpy.ndarray
            Total atoms concentrations
        res : float
            Root of the search_for_keff function
        Returns
        ------------
        x : list of numpy.ndarray
            Updated total atoms concentrations
        """

        for i, mat in enumerate(self.burn_mats):
            density = 0
            for j, nuc in enumerate(self.operator.number.burnable_nuclides):
                if nuc in self.refuel_vector.keys() and int(mat) == self.mat_id:
                    # Convert res grams into atoms
                    res_atoms = res / atomic_mass(nuc) * AVOGADRO * \
                                self.refuel_vector[nuc]
                    x[i][j] += res_atoms
                # Density in [atoms-g/mol]
                density += x[i][j] * atomic_mass(nuc)

            # Mass density in [g/cc] that will be kept constant
            mass_dens = [m.get_mass_density() for m in self.model.materials if
                    m.id == int(mat)][0]
            #In the internal version we assign new volume to AtomNumber
            self.operator.number.volume[i] = atoms_gram_per_mol / AVOGADRO / \
                                             mass_dens

        return x

    def msr_search_for_keff(self, x, step_index):
        """
        Perform the criticality search on the parametric material model.
        Will set the root of the `search_for_keff` function to the atoms
        concentrations vector.
        Parameters
        ------------
        x : list of numpy.ndarray
            Total atoms concentrations
        Returns
        ------------
        x : list of numpy.ndarray
            Updated total atoms concentrations
        """
        if step_index > 0:
            self.operator.number.set_density(x)
            self._check_nuclides(self.refuel_vector.keys())
            x, res = super()._msr_search_for_keff(x, 0)

            print('UPDATE: material addition --> {:.2f} g --> '.format(res))
            x = self._update_x_vector_and_volumes(x, res)

            #Store results
            super()._save_res(step_index, res)
        else:
            super()._save_res(step_index, 0)
        return  x
