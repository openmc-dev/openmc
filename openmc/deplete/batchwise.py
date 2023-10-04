from abc import ABC, abstractmethod
from copy import deepcopy
from numbers import Real
from warnings import warn
import numpy as np
import re

import openmc.lib
from openmc.mpi import comm
from openmc.search import _SCALAR_BRACKETED_METHODS, search_for_keff
from openmc import Material, Cell
from openmc.data import atomic_mass, AVOGADRO
from openmc.checkvalue import (check_type, check_value, check_less_than,
    check_iterable_type, check_length)
from ._density_funcs import oxidation_state, flithu, flithtru

class Batchwise(ABC):
    """Abstract class defining a generalized batch wise scheme.

    Batchwise schemes, such as control rod adjustment or material refueling to
    control reactivity and maintain keff constant and equal to one.

    A batch wise scheme can be added here to an integrator instance,
    such as  :class:`openmc.deplete.CECMIntegrator`, to parametrize one system
    variable with the aim of satisfy certain design criteria, such as keeping
    keff equal to one, while running transport-depletion caculations.

    Specific classes for running batch wise depletion calculations are
    implemented as derived class of Batchwise
    .
    .. versionadded:: 0.13.4

    Parameters
    ----------
    operator : openmc.deplete.Operator
        OpenMC operator object
    model : openmc.model.Model
        OpenMC model object
    bracket : list of float
        Bracketing interval to search for the solution, always relative to the
        solution at previous step.
    bracket_limit : list of float
        Absolute bracketing interval lower and upper; if during the adaptive
        algorithm the search_for_keff solution lies off these limits the closest
        limit will be set as new result.
    density_treatment : str
        Whether or not to keep constant volume or density after a depletion step
        before the next one.
        Default to 'constant-volume'
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
        Whether or not to print `search_for_keff` iterations.
        Default to True
    search_for_keff_output : Bool, Optional
        Whether or not to print transport iterations during  `search_for_keff`.
        Default to False
    Attributes
    ----------
    burn_mats : list of str
        List of burnable materials ids
    local_mats : list of str
        All burnable material IDs being managed by a single process
        List of burnable materials ids
    """
    def __init__(self, operator, model, bracket, bracket_limit,
                 density_treatment = 'constant-volume', bracketed_method='brentq',
                 redox_vec=None, tol=0.01, target=1.0, print_iterations=True,
                 search_for_keff_output=True):

        self.operator = operator
        self.burn_mats = operator.burnable_mats
        self.local_mats = operator.local_mats
        self.model = model
        self.geometry = model.geometry

        check_value('density_treatment', density_treatment,
                                        ('constant-density','constant-volume'))
        self.density_treatment = density_treatment
        # initialize material density functions container. density_treatment
        # is ignored if density functions are set
        self.density_functions = {}

        check_iterable_type('bracket', bracket, Real)
        check_length('bracket', bracket, 2)
        check_less_than('bracket values', bracket[0], bracket[1])
        self.bracket = bracket

        check_iterable_type('bracket_limit', bracket_limit, Real)
        check_length('bracket_limit', bracket_limit, 2)
        check_less_than('bracket limit values',
                         bracket_limit[0], bracket_limit[1])
        self.bracket_limit = bracket_limit

        check_value('bracketed_method', bracketed_method,
                   _SCALAR_BRACKETED_METHODS)
        self.bracketed_method = bracketed_method

        if redox_vec is not None:
            check_type("redox vector", redox_vec, dict, str)
            for nuc in redox_vec:
                check_value("check nuclide exists", nuc, self.operator.nuclides_with_data)

            if round(sum(redox_vec.values()), 2) != 1.0:
                # Normalize material elements vector
                sum_values = sum(redox_vec.values())
                for nuc in redox_vec:
                    redox_vec[nuc] /= sum_values
        self.redox_vec = redox_vec

        check_type('tol', tol, Real)
        self.tol = tol

        check_type('target', target, Real)
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
            warn('WARNING: brentq bracketed method is recommended')
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

    @classmethod
    def from_params(cls, obj, attr, operator, model, **kwargs):
        return cls(obj, attr, operator, model, **kwargs)

    @abstractmethod
    def _model_builder(self, param):
        """
        Builds the parametric model to be passed to the
        :meth:`openmc.search.search_for_keff` method.
        Callable function which builds a model according to a passed
        parameter. This function must return an openmc.model.Model object.
        Parameters
        ----------
        param : parameter
            model function variable
        Returns
        -------
        _model :  openmc.model.Model
            OpenMC parametric model
        """

    def _search_for_keff(self, val):
        """
        Perform the criticality search for a given parametric model.
        If the solution lies off the initial bracket, this method iteratively
        adapt it until :meth:`openmc.search.search_for_keff` return a valid
        solution.
        The adapting bracket algorithm takes the ratio between the bracket and
        the returned keffs values to scale the bracket until a valid solution
        is found.
        A bracket limit pose the upper and lower boundaries to the adapting
        bracket. If one limit is hit, the algorithm will stop and the closest
        limit value will be used.
        Parameters
        ----------
        val : float
            Previous result value
        Returns
        -------
        root : float
            Estimated value of the variable parameter where keff is the
            targeted value
        """
        #make sure we don't modify original bracket and tol values
        bracket = deepcopy(self.bracket)

        #search_for_keff tolerance should vary according to the first guess value
        if abs(val) > 1.0:
            tol = self.tol / abs(val)
        else:
            tol = self.tol

        # Run until a search_for_keff root is found or out of limits
        root = None

        while root == None:
            search = search_for_keff(self._model_builder,
                            bracket = np.array(bracket) + val,
                            tol = tol,
                            bracketed_method = self.bracketed_method,
                            target = self.target,
                            print_iterations = self.print_iterations,
                            run_args = {'output': self.search_for_keff_output},
                            run_in_memory = True)

            # if len(search) is 3 search_for_keff was successful
            if len(search) == 3:
                res,_,_ = search

                #Check if root is within bracket limits
                if self.bracket_limit[0] < res <  self.bracket_limit[1]:
                    root = res

                else:
                    # Set res with the closest limit and continue
                    arg_min = abs(np.array(self.bracket_limit) - res).argmin()
                    warn("WARNING: Search_for_keff returned root out of "
                         "bracket limit. Set root to {:.2f} and continue."
                         .format(self.bracket_limit[arg_min]))
                    root = self.bracket_limit[arg_min]

            elif len(search) == 2:
                guesses, keffs = search

                #Check if all guesses are within bracket limits
                if all(self.bracket_limit[0] < guess < self.bracket_limit[1] \
                    for guess in guesses):
                    #Simple method to iteratively adapt the bracket
                    print("Search_for_keff returned values below or above "
                          "target. Trying to iteratively adapt bracket...")

                    # if the bracket ends up being smaller than the std of the
                    # keff's closer value to target, no need to continue-
                    if all(keff <= max(keffs).s for keff in keffs):
                        arg_min = abs(self.target - np.array(guesses)).argmin()
                        root = guesses[arg_min]

                    # Calculate gradient as ratio of delta bracket and delta keffs
                    grad = abs(np.diff(bracket) / np.diff(keffs))[0].n
                    # Move the bracket closer to presumed keff root.

                    # Two cases: both keffs are below or above target
                    if np.mean(keffs) < self.target:
                        # direction of moving bracket: +1 is up, -1 is down
                        if guesses[np.argmax(keffs)] > guesses[np.argmin(keffs)]:
                            dir = 1
                        else:
                            dir = -1
                        bracket[np.argmin(keffs)] = bracket[np.argmax(keffs)]
                        bracket[np.argmax(keffs)] += grad * (self.target - \
                                                  max(keffs).n) * dir
                    else:
                        if guesses[np.argmax(keffs)] > guesses[np.argmin(keffs)]:
                            dir = -1
                        else:
                            dir = 1
                        bracket[np.argmax(keffs)] = bracket[np.argmin(keffs)]
                        bracket[np.argmin(keffs)] += grad * (min(keffs).n - \
                                                  self.target) * dir

                else:
                    # Set res with closest limit and continue
                    arg_min = abs(np.array(self.bracket_limit) - guesses).argmin()
                    warn("WARNING: Adaptive iterative bracket went off "
                         "bracket limits. Set root to {:.2f} and continue."
                         .format(self.bracket_limit[arg_min]))
                    root = self.bracket_limit[arg_min]

            else:
                raise ValueError('ERROR: Search_for_keff output is not valid')

        return root

    def _get_material(self, val):
        """Helper method for getting openmc material from Material instance or
        material name or id.
        Parameters
        ----------
        val : Openmc.Material or str or int representing material name or id
        Returns
        -------
        val : openmc.Material
            Openmc Material
        """
        if isinstance(val, Material):
            check_value('Material', str(val.id), self.burn_mats)

        elif isinstance(val, str):
            if val.isnumeric():
                check_value('Material id', val, self.burn_mats)
                val = [mat for mat in self.model.materials \
                        if mat.id == int(val)][0]

            else:
                check_value('Material name', val,
                    [mat.name for mat in self.model.materials if mat.depletable])
                val = [mat for mat in self.model.materials \
                        if mat.name == val][0]

        elif isinstance(val, int):
            check_value('Material id', str(val), self.burn_mats)
            val = [mat for mat in self.model.materials \
                    if mat.id == val][0]

        else:
            ValueError(f'Material: {val} is not supported')

        return val

    def _update_volumes(self):
        """
        Update volumes stored in AtomNumber.
        After a depletion step, both material volume and density change, due to
        changes in nuclides composition.
        At present we lack an implementation to calculate density and volume
        changes due to the different molecules speciation. Therefore, OpenMC
        assumes by default that the depletable volume does not change and only
        updates the nuclide densities and consequently the total material density.
        This method, vice-versa, assumes that the total material density does not
        change and update the material volumes instead.
        """
        number_i = self.operator.number
        for mat_idx, mat_id in enumerate(self.local_mats):
            # Total number of atoms-gram per mol
            agpm = 0
            for nuc in number_i.nuclides:
                agpm +=  number_i[mat_id, nuc] * atomic_mass(nuc)
            # Get mass dens from beginning, intended to be held constant
            density = openmc.lib.materials[int(mat_id)].get_density('g/cm3')
            number_i.volume[mat_idx] = agpm / AVOGADRO / density

    def set_density_function(self, mat, density_func):
        #check mat exists and is depletable
        mat = self._get_material(mat)
        # check density_func exists
        check_value('density_func', density_func, ('flithu', 'flithtru'))
        self.density_functions[str(mat.id)] = density_func

    def _update_densities(self):
        for mat in self.local_mats:
            if mat in self.density_functions:
                mol_comp = self._get_molar_composition(mat)
                if self.density_functions[mat] == 'flithu':
                    density = flithu(**mol_comp)
                elif self.density_functions[mat] == 'flithtru':
                    density = flithtru(**mol_comp)
                if not np.isnan(density):
                    openmc.lib.materials[int(mat)].set_density(density,'g/cm3')
                else:
                    print(f'Mat: {mat}, with comp: {mol_comp}, returned nan density')

    def _update_materials(self, x):
        """
        Update number density and material compositions in OpenMC on all processes.
        If density_treatment is set to 'constant-density'
        :meth:`openmc.deplete.batchwise._update_volumes` is called to update
        material volumes in AtomNumber, keeping the material total density
        constant, before re-normalizing the atom densities and assigning them
        to the model in memory.
        Parameters
        ----------
        x : list of numpy.ndarray
            Total atom concentrations
        """
        self.operator.number.set_density(x)

        #if self.redox_vec is not None:
        #    x = self._balance_redox(x)
        if self.density_functions:
            self._update_densities()
            self._update_volumes()
        elif self.density_treatment == 'constant-density':
            self._update_volumes()
        #else
            # constant-volume

        for rank in range(comm.size):
            number_i = comm.bcast(self.operator.number, root=rank)

            for mat in number_i.materials:
                nuclides = []
                densities = []

                for nuc in number_i.nuclides:
                    # get atom density in atoms/b-cm
                    val = 1.0e-24 * number_i.get_atom_density(mat, nuc)
                    if nuc in self.operator.nuclides_with_data:
                        if val > 0.0:
                            nuclides.append(nuc)
                            densities.append(val)
                # Update densities on C API side
                openmc.lib.materials[int(mat)].set_densities(nuclides, densities)

    def _update_x_and_set_volumes(self, x, volumes):
        """
        Update x vector with new volumes, before assign them to AtomNumber.
        Parameters
        ----------
        x : list of numpy.ndarray
            Total atoms concentrations
        volumes : dict
            Updated volumes, where key is material id and value material
            volume, in cm3
        Returns
        -------
        x : list of numpy.ndarray
            Updated total atoms concentrations
        """
        number_i = self.operator.number

        for mat_idx, mat in enumerate(self.local_mats):

            if mat in volumes:
                res_vol = volumes[mat]
                # Normalize burnable nuclides in x vector without cross section data
                for nuc_idx, nuc in enumerate(number_i.burnable_nuclides):
                    if nuc not in self.operator.nuclides_with_data:
                        # normalize with new volume
                        x[mat_idx][nuc_idx] *= number_i.get_mat_volume(mat) / \
                                               res_vol

                # update all concentration data with the new updated volumes
                for nuc, dens in zip(openmc.lib.materials[int(mat)].nuclides,
                                     openmc.lib.materials[int(mat)].densities):
                    nuc_idx = number_i.index_nuc[nuc]
                    n_of_atoms = dens / 1.0e-24 * res_vol

                    if nuc in number_i.burnable_nuclides:
                        # convert [#atoms/b-cm] into [#atoms]
                        x[mat_idx][nuc_idx] = n_of_atoms
                    # when the nuclide is not in depletion chain update the AtomNumber
                    else:
                        #Atom density needs to be in [#atoms/cm3]
                        number_i[mat, nuc] = n_of_atoms

                #Now we can update the new volume in AtomNumber
                number_i.volume[mat_idx] = res_vol
        return x

    def _get_molar_composition(self, mat_id):
        """
        Parameters
        ----------
        mats : list of openmc.materials
            materials where to calculate molar composition on
        Return
        ------
         mat_molar_comp: dict of dict of floats
            dictionary where keys are material id and values dict of molar
            compositions, where keys are element and values molar fractions
        """
        number_i = self.operator.number
        mol_comp = {}
        for mat in self.local_mats:
            if mat == mat_id:
                #list of TRU
                tru = ['Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr']
                # atoms per fluoride element
                atoms_per_elm = {elm:0 for elm in ['Li', 'Th', 'U', 'TRU']}
                for nuc in number_i.nuclides:
                    elm = re.split(r'\d+', nuc)[0]
                    if elm in ['Li', 'Th', 'U']:
                        # alwyas consider 1 anion of Fluorine
                        ions = oxidation_state[elm] + 1
                        atoms_per_elm[elm] += number_i[mat, nuc] * ions
                    elif elm in tru:
                        ions = oxidation_state[elm] + 1
                        atoms_per_elm['TRU'] += number_i[mat, nuc] * ions

                mol_comp = {k:100*v/sum(atoms_per_elm.values()) for k,v in \
                                                atoms_per_elm.items()}

        return mol_comp


    def _get_redox(self, mat_id):
        """
        Parameters
        ----------
        mats : list of openmc.materials
            materials where to calculate excess (definciency) of fluorine atoms
        Return
        ------
        excess_f_dict : dict of float
            dictionary where keys are material id and values excess (definciency)
            of fluorine atoms
        """
        # initialize excess fluorine for each material in list
        redox = 0
        number_i = self.operator.number

        for mat_idx, mat in enumerate(self.local_mats):
            if mat == mat_id:
                for nuc in number_i.nuclides:
                    elm = re.split(r'\d+', nuc)[0]
                    if elm in oxidation_state:
                    #excess (definciency) or fluorine is given by the number of atoms
                    # of each element multiplied by its oxidation state, assuming
                    # every element whats to bind with fluorine.
                        redox += number_i[mat, nuc] * oxidation_state[elm]
        return redox

    def _balance_redox(self, x):
        number_i = self.operator.number
        for mat_idx,mat in enumerate(self.local_mats):
            # number of fluorine atoms to balance
            redox = self._get_redox(mat)
            # excess/deficiency of F, add or remove mat_vec nuclides to balance it
            # if redox is negative there is an excess of fluorine and we need to
            # add mat_vector, if positive we need to remove it.
            for nuc,fraction in self.redox_vec.items():
                elm = re.split(r'\d+', nuc)[0]
                number_i[mat, nuc] -= redox * fraction / oxidation_state[elm]
            #Also update x vector
            nuc_idx = number_i.index_nuc[nuc]
            x[mat_idx][nuc_idx] = number_i[mat, nuc]
        return x

class BatchwiseCell(Batchwise):
    """Abstract class holding batch wise cell-based functions.

    Specific classes for running batch wise depletion calculations are
    implemented as derived class of BatchwiseCell.

    .. versionadded:: 0.13.4

    Parameters
    ----------
    cell : openmc.Cell or int or str
        OpenMC Cell identifier to where apply batch wise scheme
    operator : openmc.deplete.Operator
        OpenMC operator object
    model : openmc.model.Model
        OpenMC model object
    bracket : list of float
        Bracketing interval to search for the solution, always relative to the
        solution at previous step.
    bracket_limit : list of float
        Absolute bracketing interval lower and upper; if during the adaptive
        algorithm the search_for_keff solution lies off these limits the closest
        limit will be set as new result.
    density_treatment : str
        Whether or not to keep constant volume or density after a depletion step
        before the next one.
        Default to 'constant-volume'
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
        Whether or not to print `search_for_keff` iterations.
        Default to True
    search_for_keff_output : Bool, Optional
        Whether or not to print transport iterations during  `search_for_keff`.
        Default to False
    Attributes
    ----------
    cell : openmc.Cell or int or str
        OpenMC Cell identifier to where apply batch wise scheme
    cell_materials : list of openmc.Material
        Depletable materials that fill the Cell Universe. Only valid for
        translation or rotation attributes
    axis : int {0,1,2}
        Directional axis for geometrical parametrization, where 0, 1 and 2 stand
        for 'x', 'y' and 'z', respectively.
    """
    def __init__(self, cell, operator, model, bracket, bracket_limit, axis=None,
                 density_treatment='constant-volume', bracketed_method='brentq',
                 redox_vec=None, tol=0.01, target=1.0, print_iterations=True,
                 search_for_keff_output=True):

        super().__init__(operator, model, bracket, bracket_limit, density_treatment,
                         bracketed_method, redox_vec, tol, target, print_iterations,
                         search_for_keff_output)

        self.cell = self._get_cell(cell)

        if axis is not None:
            #index of cell directional axis
            check_value('axis', axis, (0,1,2))
        self.axis = axis

        # list of material fill the attribute cell, if depletables
        self.cell_materials = None

    @abstractmethod
    def _get_cell_attrib(self):
        """
        Get cell attribute coefficient.
        Returns
        -------
        coeff : float
            cell coefficient
        """

    @abstractmethod
    def _set_cell_attrib(self, val):
        """
        Set cell attribute to the cell instance.
        Parameters
        ----------
        val : float
            cell coefficient to set
        """

    def _get_cell(self, val):
        """Helper method for getting cell from cell instance or cell name or id.
        Parameters
        ----------
        val : Openmc.Cell or str or int representing Cell
        Returns
        -------
        val : openmc.Cell
            Openmc Cell
        """
        if isinstance(val, Cell):
            check_value('Cell exists', val, [cell for cell in \
                                self.geometry.get_all_cells().values()])

        elif isinstance(val, str):
            if val.isnumeric():
                check_value('Cell id exists', int(val), [cell.id for cell in \
                                self.geometry.get_all_cells().values()])
                val = [cell for cell in \
                       self.geometry.get_all_cells().values() \
                       if cell.id == int(val)][0]
            else:
                check_value('Cell name exists', val, [cell.name for cell in \
                                self.geometry.get_all_cells().values()])

                val = [cell for cell in \
                       self.geometry.get_all_cells().values() \
                       if cell.name == val][0]

        elif isinstance(val, int):
            check_value('Cell id exists', val, [cell.id for cell in \
                                self.geometry.get_all_cells().values()])
            val = [cell for cell in \
                       self.geometry.get_all_cells().values() \
                       if cell.id == val][0]
        else:
            ValueError(f'Cell: {val} is not supported')

        return val

    def _model_builder(self, param):
        """
        Builds the parametric model to be passed to the
        :meth:`openmc.search.search_for_keff` method.
        Callable function which builds a model according to a passed
        parameter. This function must return an openmc.model.Model object.
        Parameters
        ----------
        param : model parametric variable
            cell translation coefficient or cell rotation coefficient
        Returns
        -------
        self.model :  openmc.model.Model
            Openmc parametric model
        """
        self._set_cell_attrib(param)
        # At this stage it not important to reassign the new volume, as the
        # nuclide densities remain constant. However, if at least one of the cell
        # materials is set as depletable, the volume change needs to be accounted for
        # as an increase or reduction of number of atoms, i.e. vector x, before
        # solving the depletion equations
        return self.model

    def search_for_keff(self, x, step_index):
        """
        Perform the criticality search on the parametric cell coefficient and
        update materials accordingly.
        The :meth:`openmc.search.search_for_keff` solution is then set as the
        new cell attribute.
        Parameters
        ----------
        x : list of numpy.ndarray
            Total atoms concentrations
        Returns
        -------
        x : list of numpy.ndarray
            Updated total atoms concentrations
        """
        # Get cell attribute from previous iteration
        val = self._get_cell_attrib()
        check_type('Cell coeff', val, Real)

        # Update all material densities from concentration vectors
        #before performing the search_for_keff. This is needed before running
        # the  transport equations in the search_for_keff algorithm
        super()._update_materials(x)

        # Calculate new cell attribute
        root = super()._search_for_keff(val)

        # set results value as attribute in the geometry
        self._set_cell_attrib(root)

        # if at least one of the cell materials is depletable, calculate new
        # volume and update x and number accordingly
        # new volume
        if self.cell_materials:
            volumes = self._calculate_volumes()
            x = super()._update_x_and_set_volumes(x, volumes)

        return x, root


class BatchwiseCellGeometrical(BatchwiseCell):
    """
    Batch wise cell-based with geometrical-attribute class.

    A user doesn't need to call this class directly.
    Instead an instance of this class is automatically created by calling
    :meth:`openmc.deplete.Integrator.add_batchwise` method from an integrator
    class, such as  :class:`openmc.deplete.CECMIntegrator`.

    .. versionadded:: 0.13.4

    Parameters
    ----------
    cell : openmc.Cell or int or str
        OpenMC Cell identifier to where apply batch wise scheme
    attrib_name : str
        Cell attribute type
    operator : openmc.deplete.Operator
        OpenMC operator object
    model : openmc.model.Model
        OpenMC model object
    axis : int {0,1,2}
        Directional axis for geometrical parametrization, where 0, 1 and 2 stand
        for 'x', 'y' and 'z', respectively.
    bracket : list of float
        Bracketing interval to search for the solution, always relative to the
        solution at previous step.
    bracket_limit : list of float
        Absolute bracketing interval lower and upper; if during the adaptive
        algorithm the search_for_keff solution lies off these limits the closest
        limit will be set as new result.
    density_treatment : str
        Whether or not to keep constant volume or density after a depletion step
        before the next one.
        Default to 'constant-volume'
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
        Whether or not to print `search_for_keff` iterations.
        Default to True
    search_for_keff_output : Bool, Optional
        Whether or not to print transport iterations during  `search_for_keff`.
        Default to False
    Attributes
    ----------
    cell : openmc.Cell or int or str
        OpenMC Cell identifier to where apply batch wise scheme
    attrib_name : str {'translation', 'rotation'}
        Cell attribute type
    axis : int {0,1,2}
        Directional axis for geometrical parametrization, where 0, 1 and 2 stand
        for 'x', 'y' and 'z', respectively.
    samples : int
        Number of samples used to generate volume estimates for stochastic
        volume calculations.
    """
    def __init__(self, cell, attrib_name, operator, model, bracket,
                 bracket_limit, axis, density_treatment='constant-volume',
                 bracketed_method='brentq', redox_vec=None, tol=0.01, target=1.0, samples=1000000,
                 print_iterations=True, search_for_keff_output=True):

        super().__init__(cell, operator, model, bracket, bracket_limit, axis,
                         density_treatment, bracketed_method, redox_vec, tol, target,
                         print_iterations, search_for_keff_output)

        check_value('attrib_name', attrib_name,
                    ('rotation', 'translation'))
        self.attrib_name = attrib_name

        # check if cell is filled with 2 cells
        check_length('fill materials', self.cell.fill.cells, 2)

        # Initialize vector
        self.vector = np.zeros(3)

        self.cell_materials = [cell.fill for cell in \
                        self.cell.fill.cells.values() if cell.fill.depletable]

        check_type('samples', samples, int)
        self.samples = samples

        if self.cell_materials:
            self._initialize_volume_calc()

    def _get_cell_attrib(self):
        """
        Get cell attribute coefficient.
        Returns
        -------
        coeff : float
            cell coefficient
        """
        for cell in openmc.lib.cells.values():
            if cell.id == self.cell.id:
                if self.attrib_name == 'translation':
                    return cell.translation[self.axis]
                elif self.attrib_name == 'rotation':
                    return cell.rotation[self.axis]

    def _set_cell_attrib(self, val):
        """
        Set cell attribute to the cell instance.
        Attributes are only applied to a cell filled with a universe containing
        two cells itself.
        Parameters
        ----------
        val : float
            Cell coefficient to set, in cm for translation and deg for rotation
        """
        self.vector[self.axis] = val
        for cell in openmc.lib.cells.values():
            if cell.id == self.cell.id:
                setattr(cell, self.attrib_name, self.vector)

    def _initialize_volume_calc(self):
        """
        Set volume calculation model settings of depletable materials filling
        the parametric Cell.
        """
        ll, ur = self.geometry.bounding_box
        mat_vol = openmc.VolumeCalculation(self.cell_materials, self.samples,
                                           ll, ur)
        self.model.settings.volume_calculations = mat_vol

    def _calculate_volumes(self):
        """
        Perform stochastic volume calculation
        Returns
        -------
        volumes : dict
            Dictionary of calculated volumes, where key is mat id and value
            material volume, in cm3
        """
        openmc.lib.calculate_volumes()
        volumes = {}
        if comm.rank == 0:
            res = openmc.VolumeCalculation.from_hdf5('volume_1.h5')
            for mat_idx, mat in enumerate(self.local_mats):
                if int(mat) in [mat.id for mat in self.cell_materials]:
                    volumes[mat] = res.volumes[int(mat)].n
        volumes = comm.bcast(volumes)
        return volumes

class BatchwiseCellTemperature(BatchwiseCell):
    """
    Batch wise cell-based with temperature-attribute class.

    A user doesn't need to call this class directly.
    Instead an instance of this class is automatically created by calling
    :meth:`openmc.deplete.Integrator.add_batchwise` method from an integrator
    class, such as  :class:`openmc.deplete.CECMIntegrator`.

    .. versionadded:: 0.13.4

    Parameters
    ----------
    cell : openmc.Cell or int or str
        OpenMC Cell identifier to where apply batch wise scheme
    operator : openmc.deplete.Operator
        OpenMC operator object
    model : openmc.model.Model
        OpenMC model object
    attrib_name : str
        Cell attribute type
    bracket : list of float
        Bracketing interval to search for the solution, always relative to the
        solution at previous step.
    bracket_limit : list of float
        Absolute bracketing interval lower and upper; if during the adaptive
        algorithm the search_for_keff solution lies off these limits the closest
        limit will be set as new result.
    density_treatment : str
        Whether or not to keep constant volume or density after a depletion step
        before the next one.
        Default to 'constant-volume'
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
        Whether or not to print `search_for_keff` iterations.
        Default to True
    search_for_keff_output : Bool, Optional
        Whether or not to print transport iterations during  `search_for_keff`.
        Default to False
    Attributes
    ----------
    cell : openmc.Cell or int or str
        OpenMC Cell identifier to where apply batch wise scheme
    """
    def __init__(self, cell, attrib_name, operator, model, bracket, bracket_limit,
                 axis=None, density_treatment='constant-volume',
                 bracketed_method='brentq', redox_vec=None, tol=0.01, target=1.0,
                 print_iterations=True, search_for_keff_output=True):

        super().__init__(cell, operator, model, bracket, bracket_limit, axis,
                         density_treatment, bracketed_method, redox_vec, tol, target,
                         print_iterations, search_for_keff_output)

        # Not needed but used for consistency with other classes
        check_value('attrib_name', attrib_name, 'temperature')
        self.attrib_name = attrib_name
        # check if initial temperature has been set to right cell material
        if isinstance(self.cell.fill, openmc.Universe):
            cells = [cell for cell in self.cell.fill.cells.values()  \
                         if cell.fill.temperature]
            check_length('Only one cell with temperature',cells,1)
            self.cell = cells[0]

        check_type('temperature cell real', self.cell.fill.temperature, Real)

    def _get_cell_attrib(self):
        """
        Get cell temperature.
        Returns
        -------
        coeff : float
            cell temperature, in Kelvin
        """
        for cell in openmc.lib.cells.values():
            if cell.id == self.cell.id:
                return cell.get_temperature()

    def _set_cell_attrib(self, val):
        """
        Set temperature value to the cell instance.
        Parameters
        ----------
        val : float
            Cell temperature to set, in Kelvin
        """
        for cell in openmc.lib.cells.values():
            if cell.id == self.cell.id:
                cell.set_temperature(val)

class BatchwiseMaterial(Batchwise):
    """Abstract class holding batch wise material-based functions.

    Specific classes for running batch wise depletion calculations are
    implemented as derived class of BatchwiseMaterial.

    .. versionadded:: 0.13.4

    Parameters
    ----------
    material : openmc.Material or int or str
        OpenMC Material identifier to where apply batch wise scheme
    operator : openmc.deplete.Operator
        OpenMC operator object
    model : openmc.model.Model
        OpenMC model object
    mat_vector : dict
        Dictionary of material composition to parameterize, where a pair key value
        represents a nuclide and its weight fraction, respectively.
    bracket : list of float
        Bracketing interval to search for the solution, always relative to the
        solution at previous step.
    bracket_limit : list of float
        Absolute bracketing interval lower and upper; if during the adaptive
        algorithm the search_for_keff solution lies off these limits the closest
        limit will be set as new result.
    density_treatment : str
        Whether or not to keep constant volume or density after a depletion step
        before the next one.
        Default to 'constant-volume'
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
        Whether or not to print `search_for_keff` iterations.
        Default to True
    search_for_keff_output : Bool, Optional
        Whether or not to print transport iterations during  `search_for_keff`.
        Default to False
    Attributes
    ----------
    material : openmc.Material or int or str
        OpenMC Material identifier to where apply batch wise scheme
    mat_vector : dict
        Dictionary of material composition to parameterize, where a pair key value
        represents a nuclide and its weight fraction, respectively.
    """

    def __init__(self, material, operator, model, mat_vector, bracket,
                 bracket_limit, density_treatment='constant-volume',
                 bracketed_method='brentq', redox_vec=None, tol=0.01, target=1.0,
                 print_iterations=True, search_for_keff_output=True):

        super().__init__(operator, model, bracket, bracket_limit,
                         density_treatment, bracketed_method, redox_vec, tol, target,
                         print_iterations, search_for_keff_output)

        self.material = super()._get_material(material)

        check_type("material vector", mat_vector, dict, str)
        for nuc in mat_vector.keys():
            check_value("check nuclide exists", nuc, self.operator.nuclides_with_data)

        if round(sum(mat_vector.values()), 2) != 1.0:
            # Normalize material elements vector
            sum_values = sum(mat_vector.values())
            for elm in mat_vector:
                mat_vector[elm] /= sum_values
        self.mat_vector = mat_vector

    @abstractmethod
    def _model_builder(self, param):
        """
        Builds the parametric model to be passed to the
        :meth:`openmc.search.search_for_keff` method.
        Callable function which builds a model according to a passed
        parameter. This function must return an openmc.model.Model object.
        Parameters
        ----------
        param :
            Model material function variable
        Returns
        -------
        _model :  openmc.model.Model
            Openmc parametric model
        """

    def search_for_keff(self, x, step_index):
        """
        Perform the criticality search on parametric material variable.
        The :meth:`openmc.search.search_for_keff` solution is then used to
        calculate the new material volume and update the atoms concentrations.
        Parameters
        ----------
        x : list of numpy.ndarray
            Total atoms concentrations
        Returns
        -------
        x : list of numpy.ndarray
            Updated total atoms concentrations
        """
        # Update AtomNumber with new conc vectors x. Materials are also updated
        # even though they are also re-calculated when running the search_for_kef
        super()._update_materials(x)

        # Solve search_for_keff and find new value
        root = super()._search_for_keff(0)

        #Update concentration vector and volumes with new value
        volumes = self._calculate_volumes(root)
        x = super()._update_x_and_set_volumes(x, volumes)

        return  x, root

class BatchwiseMaterialRefuel(BatchwiseMaterial):
    """
    Batch wise material-based class for refuelling (addition or removal) scheme.

    A user doesn't need to call this class directly.
    Instead an instance of this class is automatically created by calling
    :meth:`openmc.deplete.Integrator.add_batchwise` method from an integrator
    class, such as  :class:`openmc.deplete.CECMIntegrator`.

    .. versionadded:: 0.13.4

    Parameters
    ----------
    material : openmc.Material or int or str
        OpenMC Material identifier to where apply batch wise scheme
    operator : openmc.deplete.Operator
        OpenMC operator object
    attrib_name : str
        Material attribute
    model : openmc.model.Model
        OpenMC model object
    mat_vector : dict
        Dictionary of material composition to parameterize, where a pair key value
        represents a nuclide and its weight fraction, respectively.
    bracket : list of float
        Bracketing interval to search for the solution, always relative to the
        solution at previous step.
    bracket_limit : list of float
        Absolute bracketing interval lower and upper; if during the adaptive
        algorithm the search_for_keff solution lies off these limits the closest
        limit will be set as new result.
    density_treatment : str
        Whether or not to keep constant volume or density after a depletion step
        before the next one.
        Default to 'constant-volume'
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
        Whether or not to print `search_for_keff` iterations.
        Default to True
    search_for_keff_output : Bool, Optional
        Whether or not to print transport iterations during  `search_for_keff`.
        Default to False
    Attributes
    ----------
    material : openmc.Material or int or str
        OpenMC Material identifier to where apply batch wise scheme
    mat_vector : dict
        Dictionary of material composition to parameterize, where a pair key value
        represents a nuclide and its weight fraction, respectively.
    """

    def __init__(self, material, attrib_name, operator, model, mat_vector, bracket,
                 bracket_limit, density_treatment='constant-volume',
                 bracketed_method='brentq', redox_vec=None, tol=0.01, target=1.0,
                 print_iterations=True, search_for_keff_output=True):

        super().__init__(material, operator, model, mat_vector, bracket,
                         bracket_limit, density_treatment, bracketed_method, redox_vec, tol,
                         target, print_iterations, search_for_keff_output)

        # Not needed but used for consistency with other classes
        check_value('attrib_name', attrib_name, 'refuel')
        self.attrib_name = attrib_name

    def _model_builder(self, param):
        """
        Callable function which builds a model according to a passed
        parameter. This function must return an openmc.model.Model object.
        Builds the parametric model to be passed to the
        :meth:`openmc.search.search_for_keff` method.
        The parametrization can either be at constant volume or constant
        density, according to user input.
        Default is constant volume.
        Parameters
        ----------
        param : float
            Model function variable, total grams of material to add or remove
        Returns
        -------
        model :  openmc.model.Model
            Openmc parametric model
        """
        for rank in range(comm.size):
            number_i = comm.bcast(self.operator.number, root=rank)

            for mat in number_i.materials:
                nuclides = []
                densities = []

                if int(mat) == self.material.id:

                    if self.density_treatment == 'constant-density':
                        vol = number_i.get_mat_volume(mat) + \
                                    (param / self.material.get_mass_density())

                    elif self.density_treatment == 'constant-volume':
                        vol = number_i.get_mat_volume(mat)

                    for nuc in number_i.index_nuc:
                        # check only nuclides with cross sections data
                        if nuc in self.operator.nuclides_with_data:
                            if nuc in self.mat_vector:
                                # units [#atoms/cm-b]
                                val = 1.0e-24 * (number_i.get_atom_density(mat,
                                    nuc) + param / atomic_mass(nuc) * \
                                    AVOGADRO * self.mat_vector[nuc] / vol)

                            else:
                                # get normalized atoms density in [atoms/b-cm]
                                val = 1.0e-24 * number_i[mat, nuc] / vol

                            if val > 0.0:
                                nuclides.append(nuc)
                                densities.append(val)

                else:
                    # for all other materials, still check atom density limits
                    for nuc in number_i.nuclides:
                        if nuc in self.operator.nuclides_with_data:
                            # get normalized atoms density in [atoms/b-cm]
                            val = 1.0e-24 * number_i.get_atom_density(mat, nuc)

                            if val > 0.0:
                                nuclides.append(nuc)
                                densities.append(val)

                #set nuclides and densities to the in-memory model
                openmc.lib.materials[int(mat)].set_densities(nuclides, densities)

        # always need to return a model
        return self.model

    def _calculate_volumes(self, res):
        """
        Uses :meth:`openmc.batchwise._search_for_keff` solution as grams of
        material to add or remove to calculate new material volume.
        Parameters
        ----------
        res : float
            Solution in grams of material, coming from
            :meth:`openmc.batchwise._search_for_keff`
        Returns
        -------
        volumes : dict
            Dictionary of calculated volume, where key mat id and value
            material volume, in cm3
        """
        number_i = self.operator.number
        volumes = {}

        for mat in self.local_mats:
            if int(mat) == self.material.id:
                if self.density_treatment == 'constant-density':
                    volumes[mat] =  number_i.get_mat_volume(mat) + \
                              (res / self.material.get_mass_density())
                elif self.density_treatment == 'constant-volume':
                    volumes[mat] = number_i.get_mat_volume(mat)
        return volumes
