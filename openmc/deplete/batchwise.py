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
from openmc.data import atomic_mass, AVOGADRO, ELEMENT_SYMBOL
from openmc.checkvalue import (check_type, check_value, check_less_than,
    check_iterable_type, check_length)
from ._density_funcs import flithu, flithtru

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
        # List of materials to add as single action in time
        self.materials_to_add = {}

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
    def get_root(self):
        """
        """

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

    @abstractmethod
    def update_from_restart(self, step_index, x, root):
        """
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

                #Check if all guesses are within bracket limits + bracket dimensions
                if all(self.bracket_limit[0]-self.bracket[0] <= guess <= self.bracket_limit[1]+self.bracket[1] \
                    for guess in guesses):
                    #Simple method to iteratively adapt the bracket
                    print("Search_for_keff returned values below or above "
                          "target. Trying to iteratively adapt bracket...")

                    # if the difference between keffs is smaller than 1 pcm,
                    # the grad calculation will be overshot, so let's set the root
                    # to the closest guess value
                    if abs(np.diff(keffs))*1e5 < 1:
                        arg_min = abs(self.target - np.array(keffs)).argmin()
                        print("Difference between keff values is "
                              "smaller than 1 pcm. "
                              "Set root to guess value with "
                              "closest keff to target and continue...")
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

                    #check if new bracket lies completely outside of limit and in that case
                    # set root to closest limit and continue
                    msg = ("WARNING: Adaptive iterative bracket {} went off "
                           "bracket limits. Set root to {:.2f} and continue."
                           )
                    if all(np.array(bracket)+val <= self.bracket_limit[0]):
                        warn(msg.format(bracket, self.bracket_limit[0]))
                        root = self.bracket_limit[0]
                    if all(np.array(bracket)+val >= self.bracket_limit[1]):
                        warn(msg.format(bracket, self.bracket_limit[1]))
                        root = self.bracket_limit[1]

                    #check if one bracket is outside of limit
                    if bracket[1] + val > self.bracket_limit[1]:
                        bracket[1] = self.bracket_limit[1] - val
                    if bracket[0] + val < self.bracket_limit[0]:
                        bracket[0] = self.bracket_limit[0] - val

                else:
                    # Set res with closest limit and continue
                    arg_min = abs(np.array(self.bracket_limit) - guesses).argmin()
                    warn("WARNING: Search_for_keff returned values off "
                        "bracket limits. Set root to {:.2f} and continue."
                         .format(self.bracket_limit[arg_min]))
                    root = self.bracket_limit[arg_min]

            else:
                raise ValueError('ERROR: Search_for_keff output is not valid')

        return root

    def _get_materials(self, vals):
        """Helper method for getting openmc material from Material instance or
        material name or id.
        Parameters
        ----------
        val : Openmc.Material or str or int representing material name or id
        Returns
        -------
        mat_val : Dict of openmc.Material and their ids
        """
        mat_val = {}
        for val in vals:
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

            mat_val[str(val.id)] = val
        return mat_val

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

    def set_density_function(self, mats, density_func, oxidation_states):
        #check mat exists and is depletable
        mats = self._get_materials(mats)
        # check density_func exists
        check_value('density_func', density_func, ('flithu', 'flithtru'))
        for mat_id in mats:
            self.density_functions[mat_id] = density_func
        for elm in oxidation_states:
            if elm not in ELEMENT_SYMBOL.values():
                raise ValueError(f'{elm} is not a valid element.')
        self.oxidation_states = oxidation_states

    def _update_densities(self):
        for mat_id in self.local_mats:
            if mat_id in self.density_functions:
                mol_comp = self._get_molar_composition(mat_id)
                if self.density_functions[mat_id] == 'flithu':
                    density = flithu(**mol_comp)
                elif self.density_functions[mat_id] == 'flithtru':
                    density = flithtru(**mol_comp)
                if not np.isnan(density):
                    openmc.lib.materials[int(mat_id)].set_density(density,'g/cm3')
                else:
                    print(f'Mat ID: {mat_id}, with comp: {mol_comp}, returned nan density')

    def _update_materials(self, x, step_index):
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
        number_i = self.operator.number
        for mat_idx, mat in enumerate(self.local_mats):
            # Total number of atoms-gram per mol
            agpm = 0
            for nuc in number_i.nuclides:
                agpm +=  number_i[mat, nuc] * atomic_mass(nuc)
            # Get mass dens from beginning, intended to be held constant
            density = openmc.lib.materials[int(mat)].get_density('g/cm3')
            number_i.volume[mat_idx] = agpm / AVOGADRO / density

    def _update_materials(self, x):
        """
        Update number density and material compositions in OpenMC on all processes.
        If density_treatment is set to 'constant-density', "_update_volumes" will
        update material volumes, keeping the material total density constant, in
        AtomNumber to renormalize the atom densities before assign them to the
        model in memory.
        Parameters
        ----------
        x : list of numpy.ndarray
            Total atom concentrations
        """
        self.operator.number.set_density(x)

        if self.materials_to_add:
            for (mat, t), values in self.materials_to_add.items():
                if t == step_index:
                    x = self._add_material(x, mat, values)

        if self.density_functions:
            self._update_densities()
            self._update_volumes()
        elif self.density_treatment == 'constant-density':
            self._update_volumes()
        #else
            # constant-volume

        elif self.density_treatment == 'constant-density':
            self._update_volumes()
        #else
            # constant-volume

        for rank in range(comm.size):
            number_i = comm.bcast(self.operator.number, root=rank)

            for mat_id in number_i.materials:
                nuclides = []
                densities = []

                for nuc in number_i.nuclides:
                    # get atom density in atoms/b-cm
                    val = 1.0e-24 * number_i.get_atom_density(mat_id, nuc)
                    if nuc in self.operator.nuclides_with_data:
                        if val > 0.0:
                            nuclides.append(nuc)
                            densities.append(val)
                # Update densities on C API side
<<<<<<< HEAD
                openmc.lib.materials[int(mat_id)].set_densities(nuclides, densities)

        return x
=======
                print(mat,nuclides,densities)
                openmc.lib.materials[int(mat)].set_densities(nuclides, densities)
>>>>>>> 0c3eee56a (run batchwise redox)

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

        for mat_idx, mat_id in enumerate(self.local_mats):

            if mat_id in volumes:
                res_vol = volumes[mat_id]
                # Normalize burnable nuclides in x vector without cross section data
                for nuc_idx, nuc in enumerate(number_i.burnable_nuclides):
                    if nuc not in self.operator.nuclides_with_data:
                        # normalize with new volume
                        x[mat_idx][nuc_idx] *= number_i.get_mat_volume(mat_id) / \
                                               res_vol

                # update all concentration data with the new updated volumes
                for nuc, dens in zip(openmc.lib.materials[int(mat_id)].nuclides,
                                     openmc.lib.materials[int(mat_id)].densities):
                    nuc_idx = number_i.index_nuc[nuc]
                    n_of_atoms = dens / 1.0e-24 * res_vol

                    if nuc in number_i.burnable_nuclides:
                        # convert [#atoms/b-cm] into [#atoms]
                        x[mat_idx][nuc_idx] = n_of_atoms
                    # when the nuclide is not in depletion chain update the AtomNumber
                    else:
                        #Atom density needs to be in [#atoms/cm3]
                        number_i[mat_id, nuc] = n_of_atoms

                #Now we can update the new volume in AtomNumber
                number_i.volume[mat_idx] = res_vol
        return x

    def _get_molar_composition(self, mat_mol_id):
        """
        Parameters
        ----------
        mat_mol_id : Openmc.material Id
            Material  where to calculate molar composition on
        Return
        ------
         mat_molar_comp: dict of dict of floats
            dictionary where keys are material id and values dict of molar
            compositions, where keys are element and values molar fractions
        """
        number_i = self.operator.number
        mol_comp = {}
        for mat_id in self.local_mats:
            if mat_id == mat_mol_id:
                #list of TRU
                tru = ['Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr']
                # atoms per fluoride element
                atoms_per_elm = {elm:0 for elm in ['Li', 'Th', 'U', 'TRU']}
                for nuc in number_i.nuclides:
                    elm = re.split(r'\d+', nuc)[0]
                    if elm in ['Li', 'Th', 'U']:
                        # alwyas consider 1 anion of Fluorine
                        ions = self.oxidation_states[elm] + 1
                        atoms_per_elm[elm] += number_i[mat_id, nuc] * ions
                    elif elm in tru:
                        ions = self.oxidation_states[elm] + 1
                        atoms_per_elm['TRU'] += number_i[mat_id, nuc] * ions

                mol_comp = {k:100*v/sum(atoms_per_elm.values()) for k,v in \
                                                atoms_per_elm.items()}

        return mol_comp

    def _get_redox(self, mat_rx_id):
        """
        Parameters
        ----------
        mat_rx_id : Openmc material id
            material where to calculate redox
        Return
        ------
        excess_f_dict : dict of float
            dictionary where keys are material id and values excess (definciency)
            of fluorine atoms
        """
        # initialize excess fluorine for each material in list
        redox = 0
        number_i = self.operator.number

        for mat_id in self.local_mats:
            if mat_id == mat_rx_id:
                for nuc in number_i.nuclides:
                    elm = re.split(r'\d+', nuc)[0]
                    if elm in self.oxidation_states:
                    #excess (definciency) or fluorine is given by the number of atoms
                    # of each element multiplied by its oxidation state, assuming
                    # every element whats to bind with fluorine.
                        redox += number_i[mat_id, nuc] * self.oxidation_states[elm]
        return redox

    def _balance_redox(self, x):
        number_i = self.operator.number

        for mat_idx, mat_id in enumerate(self.local_mats):
            if mat_id != '9':
                # number of fluorine atoms to balance
                redox = self._get_redox(mat_id)
                # excess/deficiency of F, add or remove mat_vec nuclides to balance it
                # if redox is negative there is an excess of fluorine and we need to
                # add mat_vector, if positive we need to remove it.
                for nuc,fraction in self.redox_vec.items():
                    elm = re.split(r'\d+', nuc)[0]
                    number_i[mat_id, nuc] -= redox * fraction / self.oxidation_states[elm]
                #Now updates x vector
                nuc_idx = number_i.index_nuc[nuc]
                x[mat_idx][nuc_idx] = number_i[mat_id, nuc]
        return x

    def add_material(self, mat, value, mat_vector, timestep, units='grams'):
        """
        Add material entries to self.materials_to_add dictionary for later use,
        before converting grams into atoms.
        Parameters
        ----------
        mat : openmc.Material or str or int
            Material identifier
        value : float
            Material to add in units of 'units'
        mat_vector : dict of float
            Nuclide composition of material to add where keys are the nuclides
            and values the fractions.
        units : str
            Units of material value.
            Default : grams
        """
        # convert grams in atoms
        self.materials_to_add[(mat,timestep)] = dict()
        for nuc,frac in mat_vector.items():
            self.materials_to_add[(mat,timestep)][nuc] = \
                    frac * value / atomic_mass(nuc) * AVOGADRO

    def _add_material(self, x, mat, values):
        """
        Private method to add materials as number of atoms to number and x.
        Parameters
        ----------
        x : list of numpy.ndarray
            Total atoms concentrations
        mat : openmc.Material or str or int
            Material identifier
        values : dict of floats
            Numebr of atoms per nuclide to add, where keys are nuclides and
            values number of atoms.
        Return
        ------
        x : list of numpy.ndarray
            Total atoms concentrations
        """
        number_i = self.operator.number
        #get material id
        mats = self._get_materials([mat])

        for mat_idx, mat_id in enumerate(self.local_mats):
            if mat_id in mats:
                for nuc, val in values.items():
                    number_i[mat_id, nuc] += val
                    #Now updates x vector
                    if nuc in number_i.burnable_nuclides:
                        nuc_idx = number_i.index_nuc[nuc]
                        x[mat_idx][nuc_idx] += val
        return x

class BatchwisePure(Batchwise):
    def __init__(self, operator, model, density_treatment='constant-volume',
                 redox_vec=None):
        super().__init__(operator, model, [0,1], [-100,100])
        self.cell_materials = None

    def get_root(self):
        return np.nan

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
        x = super()._update_materials(x, step_index)
        # if at least one of the cell materials is depletable, calculate new
        # volume and update x and number accordingly.
        # Alternatively, if material has been edded x needs to be updated
        # TODO: improve this
        if self.cell_materials:
            volumes = self._calculate_volumes()
            x = super()._update_x_and_set_volumes(x, volumes)

        return x, np.nan

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

    def update_from_restart(self, step_index, x, root):
        """
        """

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

    def get_root(self):
        return self._get_cell_attrib()

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

class BatchwiseCell(Batchwise):
    """Abstract class holding batchwise cell-based functions.

    Specific classes for running batchwise depletion calculations are
    implemented as derived class of BatchwiseCell.

    .. versionadded:: 0.13.4

    Parameters
    ----------
    cell : openmc.Cell or int or str
        OpenMC Cell identifier to where apply batchwise scheme
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
        Wether ot not to keep contant volume or density after a depletion step
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
        Wheter or not to print `search_for_keff` iterations.
        Default to True
    search_for_keff_output : Bool, Optional
        Wheter or not to print transport iterations during  `search_for_keff`.
        Default to False
    Attributes
    ----------
    cell : openmc.Cell or int or str
        OpenMC Cell identifier to where apply batchwise scheme
    cell_materials : list of openmc.Material
        Depletable materials that fill the Cell Universe. Only valid for
        translation or rotation atttributes
    """
    def __init__(self, cell, operator, model, bracket, bracket_limit,
                 density_treatment='constant-volume', bracketed_method='brentq',
                 tol=0.01, target=1.0, print_iterations=True,
                 search_for_keff_output=True):

        super().__init__(operator, model, bracket, bracket_limit, density_treatment,
                         bracketed_method, tol, target, print_iterations,
                         search_for_keff_output)

        self.cell = self._get_cell(cell)

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
        x = super()._update_materials(x, step_index)

        # Calculate new cell attribute
        root = super()._search_for_keff(val)

        # set results value as attribute in the geometry
        self._set_cell_attrib(root)

        # if at least one of the cell materials is depletable, calculate new
        # volume and update x and number accordingly.
        # Alternatively, if material has been edded x needs to be updated
        # TODO: improve this
        if self.cell_materials:
            volumes = self._calculate_volumes()
            x = super()._update_x_and_set_volumes(x, volumes)

        return x, root

    def update_from_restart(self, step_index, x, root):
        self._set_cell_attrib(root)
        super()._update_materials(x, step_index)

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
        if not isinstance(self.cell.fill, openmc.universe.DAGMCUniverse):

            check_length('fill materials', self.cell.fill.cells, 2)
            self.cell_materials = [cell.fill for cell in \
                        self.cell.fill.cells.values() if cell.fill.depletable]
        else:
            self.cell_materials = None

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
        for cell in openmc.lib.cells.values():
            if cell.id == self.cell.id:
                if self.attrib_name == "translation":
                    vector = cell.translation
                elif self.attrib_name == "rotation":
                    vector = cell.rotation
                vector[self.axis] = val
                setattr(cell, self.attrib_name, vector)

    def _initialize_volume_calc(self):
        """
        Set volume calculation model settings of depletable materials filling
        the parametric Cell.
        """
        if not self.model.settings.volume_calculations:
            mat_vol = openmc.VolumeCalculation(self.cell_materials, self.samples,
                                           ll.values(), ur.values())
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
    materias : List of openmc.Material or int or str
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
    materials : List of openmc.Material or int or str
        List OpenMC Material identifier to where apply batch wise scheme
    mat_vector : dict
        Dictionary of material composition to parameterize, where a pair key value
        represents a nuclide and its weight fraction, respectively.
    """

    def __init__(self, materials, operator, model, mat_vector, bracket,
                 bracket_limit, density_treatment='constant-volume',
                 bracketed_method='brentq', redox_vec=None, tol=0.01, target=1.0,
                 print_iterations=True, search_for_keff_output=True):

        super().__init__(operator, model, bracket, bracket_limit,
                         density_treatment, bracketed_method, redox_vec, tol, target,
                         print_iterations, search_for_keff_output)

        self.materials = super()._get_materials(materials)

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

    def get_root(self):
        for mat in openmc.lib.materials.values():
            if mat.id in [m.id for m in self.materials]:
                return mat.get_density()

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
        x = super()._update_materials(x, step_index)

        # Solve search_for_keff and find new value
        root = super()._search_for_keff(0)

        #Update concentration vector and volumes with new value
        volumes = self._calculate_volumes(root)
        x = super()._update_x_and_set_volumes(x, volumes)

        return  x, root

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

        for mat_id in self.local_mats:
            if mat_id in self.materials:
                if self.density_treatment == 'constant-density':
                    volumes[mat_id] =  number_i.get_mat_volume(mat_id) + \
                              (res / self.materials[mat_id].get_mass_density())
                elif self.density_treatment == 'constant-volume':
                    volumes[mat_id] = number_i.get_mat_volume(mat_id)
        return volumes

    def update_from_restart(self, step_index, x, root):
        super()._update_materials(x, step_index)

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
    materials : List of openmc.Material or int or str
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
    materials : List of openmc.Material or int or str
        OpenMC Material identifier to where apply batch wise scheme
    mat_vector : dict
        Dictionary of material composition to parameterize, where a pair key value
        represents a nuclide and its weight fraction, respectively.
    """

    def __init__(self, materials, attrib_name, operator, model, mat_vector, bracket,
                 bracket_limit, density_treatment='constant-volume',
                 bracketed_method='brentq', redox_vec=None, tol=0.01, target=1.0,
                 print_iterations=True, search_for_keff_output=True):

        super().__init__(materials, operator, model, mat_vector, bracket,
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

            for mat_id in number_i.materials:
                nuclides = []
                densities = []

                if mat_id in self.materials:

                    if self.density_treatment == 'constant-density':
                        vol = number_i.get_mat_volume(mat_id) + \
                                    (param / self.materials[mat_id].get_mass_density())

                    elif self.density_treatment == 'constant-volume':
                        vol = number_i.get_mat_volume(mat_id)

                    for nuc in number_i.index_nuc:
                        # check only nuclides with cross sections data
                        if nuc in self.operator.nuclides_with_data:
                            if nuc in self.mat_vector:
                                # units [#atoms/cm-b]
                                val = 1.0e-24 * (number_i.get_atom_density(mat_id,
                                    nuc) + param / atomic_mass(nuc) * \
                                    AVOGADRO * self.mat_vector[nuc] / vol)

                            else:
                                # get normalized atoms density in [atoms/b-cm]
                                val = 1.0e-24 * number_i[mat_id, nuc] / vol

                            if val > 0.0:
                                nuclides.append(nuc)
                                densities.append(val)

                else:
                    # for all other materials, still check atom density limits
                    for nuc in number_i.nuclides:
                        if nuc in self.operator.nuclides_with_data:
                            # get normalized atoms density in [atoms/b-cm]
                            val = 1.0e-24 * number_i.get_atom_density(mat_id, nuc)

                            if val > 0.0:
                                nuclides.append(nuc)
                                densities.append(val)

                #set nuclides and densities to the in-memory model
                openmc.lib.materials[int(mat_id)].set_densities(nuclides, densities)

        # always need to return a model
        return self.model

class BatchwiseMaterialDilute(BatchwiseMaterial):
    """
    Batch wise material-based class for refuelling (addition or removal) scheme.

    A user doesn't need to call this class directly.
    Instead an instance of this class is automatically created by calling
    :meth:`openmc.deplete.Integrator.add_batchwise` method from an integrator
    class, such as  :class:`openmc.deplete.CECMIntegrator`.

    .. versionadded:: 0.13.4

    Parameters
    ----------
    materials : List of openmc.Material or int or str
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
    materials : List of openmc.Material or int or str
        OpenMC Material identifier to where apply batch wise scheme
    mat_vector : dict
        Dictionary of material composition to parameterize, where a pair key value
        represents a nuclide and its weight fraction, respectively.
    """

    def __init__(self, materials, attrib_name, operator, model, mat_vector, bracket,
                 bracket_limit, density_treatment='constant-volume',
                 bracketed_method='brentq', redox_vec=None, tol=0.01, target=1.0,
                 print_iterations=True, search_for_keff_output=True):

        super().__init__(materials, operator, model, mat_vector, bracket,
                         bracket_limit, density_treatment, bracketed_method, redox_vec, tol,
                         target, print_iterations, search_for_keff_output)

        # Not needed but used for consistency with other classes
        check_value('attrib_name', attrib_name, 'dilute')
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

            for mat_id in number_i.materials:
                nuclides = []
                densities = []

                if mat_id in self.materials:

                    if self.density_treatment == 'constant-density':
                        vol = number_i.get_mat_volume(mat_id) + \
                                    (param / self.materials[mat_id].get_mass_density())

                    elif self.density_treatment == 'constant-volume':
                        vol = number_i.get_mat_volume(mat_id)

                    # Sum all atoms present and convert in [#atoms/b-cm]
                    mat_ind = number_i.index_mat[mat_id]
                    tot_atoms = 1.0e-24 * sum(number_i.number[mat_ind]) / vol

                    for nuc in number_i.index_nuc:
                        # check only nuclides with cross sections data
                        if nuc in self.operator.nuclides_with_data:
                            # [#atoms/b-cm]
                            val = 1.0e-24 * number_i.get_atom_density(mat_id,nuc)
                            # Build parametric function, where param is the
                            # dilute fraction to replace.
                            # it assumes all nuclides in material vector have
                            # cross sections data
                            if nuc in self.mat_vector:
                                # units [#atoms/cm-b]
                                val = (1-param) * val + param * \
                                      self.mat_vector[nuc] * tot_atoms

                            else:
                                val *= (1-param)

                            if val > 0.0:
                                nuclides.append(nuc)
                                densities.append(val)

                else:
                    # for all other materials, still check atom density limits
                    for nuc in number_i.nuclides:
                        if nuc in self.operator.nuclides_with_data:
                            # get normalized atoms density in [atoms/b-cm]
                            val = 1.0e-24 * number_i.get_atom_density(mat_id, nuc)

                            if val > 0.0:
                                nuclides.append(nuc)
                                densities.append(val)

                #set nuclides and densities to the in-memory model
                openmc.lib.materials[int(mat_id)].set_densities(nuclides, densities)

        # always need to return a model
        return self.model

class BatchwiseMaterialAdd(BatchwiseMaterial):
    """
    Batch wise material-based class point-wise material addition or removal.
    Not based on criticality control but on user-defined values.

    A user doesn't need to call this class directly.
    Instead an instance of this class is automatically created by calling
    :meth:`openmc.deplete.Integrator.add_batchwise` method from an integrator
    class, such as  :class:`openmc.deplete.CECMIntegrator`.

    .. versionadded:: 0.13.4

    Parameters
    ----------
    materials : List of openmc.Material or int or str
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
    materials : List of openmc.Material or int or str
        OpenMC Material identifier to where apply batch wise scheme
    mat_vector : dict
        Dictionary of material composition to parameterize, where a pair key value
        represents a nuclide and its weight fraction, respectively.
    """

    def __init__(self, materials, attrib_name, operator, model, mat_vector, bracket,
                 bracket_limit, density, volume, density_treatment='constant-volume',
                 bracketed_method='brentq', redox_vec=None, tol=0.01, target=1.0,
                 print_iterations=True, search_for_keff_output=True):

        super().__init__(materials, operator, model, mat_vector, bracket,
                         bracket_limit, density_treatment, bracketed_method,
                         redox_vec, tol, target, print_iterations,
                         search_for_keff_output)

        check_value('attrib_name', attrib_name, 'addition')
        self.attrib_name = attrib_name

        self.density = density
        self.volume_to_add = volume

    def _model_builder(self, param):
        """
        Builds the parametric model that is passed to the `msr_search_for_keff`
        function by updating the material densities and setting the parametric
        variable to the material nuclides to add. Here we fix the total number
        of atoms per material and try to conserve this quantity. Both the
        material volume and density are let free to vary.
        Parameters
        ----------
        param :
            Model function variable, fraction of total atoms to dilute
        Returns
        -------
        _model :  openmc.model.Model
            OpenMC parametric model
        """

    def _update_densities(self, mat_id, atoms, vol):
        """
        Calculate new atom densities dividing new total atoms
        with new volume and assign them to openmc.lib.materials.
        Parameters
        ----------
        x : list of numpy.ndarray
            Total atoms concentrations
        mat_id : str
            openmc.materials id
        atoms : dict
            dictionary of update atoms concentrations
        vol : float
            new update volume
        """
        for rank in range(comm.size):
            number_i = comm.bcast(self.operator.number, root=rank)
            if str(mat_id) in number_i.materials:
                mat_idx = number_i.index_mat[str(mat_id)]
                nuclides = []
                densities = []
                for nuc in number_i.nuclides:
                    # Only modify nuclides with cross section data available
                    if nuc in self.operator.nuclides_with_data:
                        # number of atoms before volume change
                        val = number_i.get_atom_density(str(mat_id) ,nuc) * \
                                                number_i.volume[mat_idx]
                        if nuc in atoms:
                            # add atoms
                            val += atoms[nuc]
                        # obtain density [atoms/b-cm], dividing by new volume
                        if val > 0.0:
                            # divide by new volume to obtain density
                            val *= 1.0e-24 / vol
                            nuclides.append(nuc)
                            densities.append(val)
                openmc.lib.materials[mat_id].set_densities(nuclides, densities)

    def _update_x(self, x, mat_id, vol):
        """
        Update x vector with newly update atoms densities and volume.
        Parameters
        ----------
        x : list of numpy.ndarray
            Total atoms concentrations
        mat_id : str
            openmc.materials id
        vol : float
            new update volume
        Returns
        -------
        x : list of numpy.ndarray
            Total atoms concentrations
        """
        number_i = self.operator.number
        if str(mat_id) in self.local_mats:
            mat_idx = self.local_mats.index(str(mat_id))
            #update x vector
            for nuc, dens in zip(openmc.lib.materials[mat_id].nuclides,
                                 openmc.lib.materials[mat_id].densities):
                if nuc in number_i.burnable_nuclides:
                    nuc_idx = number_i.burnable_nuclides.index(nuc)
                    # convert [#atoms/b-cm] into [#atoms]
                    x[mat_idx][nuc_idx] = dens / 1.0e-24 * vol
                else:
                    #Divide by 1.0e-24, atom density here must be in [#atoms/cm3]
                    number_i.set_atom_density(mat_idx, nuc, dens / 1.0e-24)

        return x

    def _dilute_x(self, x, mat_id, vol):
        """
        Dilute atom vector x within new volume
        Parameters
        ----------
        x : list of numpy.ndarray
            Total atoms concentrations
        mat_id : str
            openmc.materials id
        vol : float
            new update volume
        Returns
        -------
        x : list of numpy.ndarray
            Total atoms concentrations
        """
        number_i = self.operator.number
        if str(mat_id) in self.local_mats:
            mat_idx = self.local_mats.index(str(mat_id))
            for nuc in number_i.burnable_nuclides:
                nuc_idx = number_i.burnable_nuclides.index(nuc)
                # update number of atoms in new volume
                x[mat_idx][nuc_idx] *= vol / number_i.volume[mat_idx]
        return x

    def _find_ofc_material_id(self, mat_name):
        """
        Giving a material name, search if exist out-of-core region materials
        and return material id.
        Parameters
        ----------
        mat_name : str
            openmc material name
        Returns
        -------
            mat_ids : list
            list of ofc materials
        """
        # split by "-" and "_" characters
        return [mat.id for mat in self.model.materials if \
            all(x in re.split(', |-|-|\+', mat.name) for x in [mat_name, 'ofc'])]

    def _calc_volumes(self):
        """
        Calculate new volumes after geometrical cell change
        Returns
        -------
        openmc.Volume : dict
            Dictionary of calculate volumes
        """
        openmc.lib.calculate_volumes()
        comm.barrier()
        return openmc.VolumeCalculation.from_hdf5('volume_1.h5')

    def _distribute_volumes(self, materials):
        """
        Following a volumes update, the new volume difference must be added
        or removed from the material compositions.
        Prepare dictionary with new volume and new total atoms for each
        material to update.
        Parameters
        ----------
        materials : list
            list of openmc.materials to update material compositions
        Returns
        -------
        mats_to_update : dict of dict
            Dictionary of materials to update where keys are material ids and
            values dictionaries, where keys are 'volume' and 'atoms' and values
            calculated new volume and atoms to update.
        """
        # firstly, calculate the new volumes
        res = self._calc_volumes()
        # Gather all AtomNumbers and x in rank 0
        number = comm.gather(self.operator.number)
        mats_to_update = {}
        # Perform all on rank 0
        if comm.rank == 0:
            for mat in materials:
                # get rank index of AtomNumber that contains mat.id
                rank_index = [i for i,n in enumerate(number) \
                                if str(mat.id) in n.materials][0]
                mat_index = number[rank_index].index_mat[str(mat.id)]
                vol_new = res.volumes[mat.id].n
                #calculate volume difference
                vol_diff = vol_new - number[rank_index].volume[mat_index]

                # add updated volume to materials in user defined list
                if str(mat.id) in self.materials:
                    #convert user defined mass density in atomic density [#atoms/cm3]
                    atom_dens = {nuc: self.density / \
                             atomic_mass(nuc) * AVOGADRO * frac \
                             for nuc, frac in self.mat_vector.items()}

                    #if the material has an ofc region, the volume_to_add user
                    # defined value is added directly to the ofc region, subtracted
                    # by the volume change in core, which is added to the in-core
                    #materials
                    if self._find_ofc_material_id(mat.name):
                        mats_to_update[mat.id] = {'vol':vol_new}
                        mats_to_update[mat.id]['atoms'] = \
                                {nuc: i*vol_diff for nuc,i in atom_dens.items()}
                        #update out of core densities using the remaining volume
                        # to self.volume_to_add
                        mat_ofc_id = self._find_ofc_material_id(mat.name)[0]
                        # mat_ofc_id might not be in the same AtomNumber chunk
                        # as mat_id
                        rank_index = [i for i,n in enumerate(number) \
                                if str(mat_ofc_id) in n.materials][0]
                        mat_ofc_index = number[rank_index].index_mat[str(mat_ofc_id)]
                        vol_ofc_new = number[rank_index].volume[mat_ofc_index] + \
                                        (self.volume_to_add-vol_diff)
                        mats_to_update[mat_ofc_id] = {'vol':vol_ofc_new}
                        mats_to_update[mat_ofc_id]['atoms'] = {nuc: i * \
                            (self.volume_to_add-vol_diff) for nuc,i in atom_dens.items()}
                    else:
                        vol_new = number[rank_index].volume[mat_index] + \
                                    self.volume_to_add
                        mats_to_update[mat.id] = {'vol':vol_new}
                        mats_to_update[mat.id]['atoms'] = {nuc: i * \
                            self.volume_to_add for nuc,i in atom_dens.items()}
                else:
                    #if the material is not in the user defined material list
                    #, check anyway if it has an ofc region. In this case the
                    #volume difference is added or removed there
                    if self._find_ofc_material_id(mat.name):
                        mat_ofc_id = self._find_ofc_material_id(mat.name)[0]
                        rank_index = [i for i,n in enumerate(number) \
                                        if str(mat_ofc_id) in n.materials][0]
                        mat_ofc_index = number[rank_index].index_mat[str(mat_ofc_id)]
                        vol_ofc_new = number[rank_index].volume[mat_ofc_index] - \
                                                                    vol_diff
                        mats_to_update[mat.id] = {'vol':vol_new}
                        mats_to_update[mat_ofc_id] = {'vol':vol_ofc_new}
                    else:
                        # the geometrical volume changes but the total volume remains
                        # unchanged
                        continue

        return mats_to_update

    def add_materials(self, x, materials):
        """
        Materials update after a geometrical change.
        Updates x composition vector, material densities and material volumes
        Parameters
        ----------
        x : list of numpy.ndarray
            Total atoms concentrations
        materials : list
            list of openmc.materials to update
        Returns
        -------
        x : list of numpy.ndarray
            Total atoms concentrations
        """
        mats_to_update = self._distribute_volumes(materials)

        for rank in range(comm.size):
            # update in all comm.ranks
            mats_to_update = comm.bcast(mats_to_update, root=rank)
            number_i = comm.bcast(self.operator.number, root=rank)
            for mat_id, vals in mats_to_update.items():
                if 'atoms' in vals:
                    self._update_densities(mat_id, vals['atoms'], vals['vol'])
                    x = self._update_x(x, mat_id, vals['vol'])
                else:
                    x = self._dilute_x(x, mat_id, vals['vol'])

                if str(mat_id) in self.operator.number.materials:
                    mat_index = self.operator.number.index_mat[str(mat_id)]
                    self.operator.number.volume[mat_index] = vals['vol']

        return x

class BatchwiseSchemeStd():
    """
    Batchwise wrapper class, it wraps BatchwiseGeom and BatchwiseMat instances,
    with some user defined logic.

    This class should probably not be defined here, but we can keep it now for
    convenience

    The loop logic of this wrapper class is the following:

    1. Run BatchwiseGeom and return geometrical coefficient
    2. check if step index equals user definde dilute interval
    3.1 if not, update geometry
    3.2 if yes, set geometrical coefficient to user-defined restart level and
    run BatchwiseMatDilute

    In this case if the bracket upper geometrical limit is hitted,
    simply stop the simulation.

    An instance of this class can be passed directly to an instance of the
    integrator class, such as :class:`openmc.deplete.CECMIntegrator`.

    Parameters
    ----------
    bw_geom : BatchwiseGeom
        openmc.deplete.batchwise.BatchwiseGeom object
    bw_mat : BatchwiseMat
        openmc.deplete.batchwise.BatchwiseMat object
    dilute_interval : int
        Frequency of dilution in number of timesteps
    first_dilute : int or None
        Timestep index for first dilution, to be used during restart simulation
        Default to None
    """

    def __init__(self, bw_list, n_timesteps, dilute_interval, restart_level,
                 first_dilute=None, interrupt=False):

        if isinstance(bw_list, list):
            for bw in bw_list:
                if isinstance(bw, BatchwiseCell):
                    self.bw_geom = bw
                elif isinstance(bw, BatchwiseMaterial):
                    self.bw_mat = bw
                else:
                    raise ValueError(f'{bw} is not a valid instance of'
                                      ' Batchwise class')
        else:
            raise ValueError(f'{bw_list} is not a list')

        self.bw_list = bw_list
        self.n_timesteps = n_timesteps

        if not isinstance(restart_level, (float, int)):
            raise ValueError(f'{restart_level} is of type {type(restart_level)},'
                             ' while it should be int or float')
        else:
            self.restart_level = restart_level

        #TODO check these values
        self.first_dilute = first_dilute
        self.step_interval = dilute_interval

        # if first dilute is set, the dilute interval needs to be updated
        if self.first_dilute is not None:
            self.dilute_interval = dilute_interval + self.first_dilute
        else:
            self.dilute_interval = dilute_interval
        self.interrupt = interrupt

    def get_root(self):
        return self.bw_geom.get_root()

    def set_density_function(self, mats, density_func, oxidation_states):
        for bw in self.bw_list:
            bw.set_density_function(mats, density_func, oxidation_states)

    def add_material(self, mat, value, mat_vector, timestep, units):
        self.bw_geom.add_material(mat, value, mat_vector, timestep, units)

    def update_from_restart(self, step_index, x, root):
        """
        This is for a restart. TODO update abc class
        Parameters
        ----------
        x : list of numpy.ndarray
            Total atom concentrations
        """
        self.bw_geom.update_from_restart(step_index, x, root)

    def search_for_keff(self, x, step_index):
        """
        Perform the criticality search on the parametric material model.
        Will set the root of the `search_for_keff` function to the atoms
        concentrations vector.
        Parameters
        ----------
        x : list of numpy.ndarray
            Total atoms concentrations
        Returns
        -------
        x : list of numpy.ndarray
            Updated total atoms concentrations
        """
        #Check if index lies in dilution timesteps
        if step_index in np.arange(self.dilute_interval, self.n_timesteps,
                                                self.dilute_interval) and \
                        self.bw_geom._get_cell_attrib() <= self.restart_level:

            # restart level and perform dilution
            self.bw_geom._set_cell_attrib(self.restart_level)
            x, _ = self.bw_mat.search_for_keff(x, step_index)
            root = self.restart_level
            #update dulution interval
            #if step_index == self.dilute_interval:
            #    self.dilute_interval += self.step_interval


        else:
            x, root = self.bw_geom.search_for_keff(x, step_index)
            # in this case if upper limit gets hit, stop directly
            if root >= self.bw_geom.bracket_limit[1]:
                from pathlib import Path
                print(f'Reached maximum of {self.bw_geom.bracket_limit[1]} cm'
                       ' exit..')
                Path('sim.done').touch()

                for rank in range(comm.size):
                    comm.Abort()
        return x, root

class BatchwiseSchemeRefuel():
    """
    Batchwise wrapper class, it wraps BatchwiseGeom and BatchwiseMat instances,
    with some user defined logic.

    This class should probably not be defined here, but we can keep it now for
    convenience

    The loop logic of this wrapper class is the following:

    1. Run BatchwiseGeom and return geometrical coefficient
    2. check if geometrical coefficient hit upper limit
    3.1 if not, continue
    3.2 if yes, refuel and reset geometrical coefficient

    An instance of this class can be passed directly to an instance of the
    integrator class, such as :class:`openmc.deplete.CECMIntegrator`.

    Parameters
    ----------
    bw_geom : BatchwiseGeom
        openmc.deplete.batchwise.BatchwiseGeom object
    bw_mat : BatchwiseMat
        openmc.deplete.batchwise.BatchwiseMat object
    restart_level : int
        Geometrical coefficient after reset
    """

    def __init__(self, bw_list, restart_level):

        if isinstance(bw_list, list):
            for bw in bw_list:
                if isinstance(bw, BatchwiseCell):
                    self.bw_geom = bw
                elif isinstance(bw, BatchwiseMaterial):
                    self.bw_mat = bw
                else:
                    raise ValueError(f'{bw} is not a valid instance of'
                                      ' Batchwise class')
        else:
            raise ValueError(f'{bw_list} is not a list')

        self.bw_list = bw_list

        if not isinstance(restart_level, (float, int)):
            raise ValueError(f'{restart_level} is of type {type(restart_level)},'
                             ' while it should be int or float')
        else:
            self.restart_level = restart_level

    def get_root(self):
        return self.bw_geom.get_root()

    def set_density_function(self, mats, density_func, oxidation_states):
        for bw in self.bw_list:
            bw.set_density_function(mats, density_func, oxidation_states)

    def add_material(self, mat, value, mat_vector, timestep, units):
        self.bw_geom.add_material(mat, value, mat_vector, timestep, units)

    def update_from_restart(self, step_index, x, root):
        """
        This is for a restart. TODO update abc class
        Parameters
        ----------
        x : list of numpy.ndarray
            Total atom concentrations
        """
        self.bw_geom.update_from_restart(step_index, x, root)

    def search_for_keff(self, x, step_index):
        """
        Perform the criticality search on the parametric material model.
        Will set the root of the `search_for_keff` function to the atoms
        concentrations vector.
        Parameters
        ----------
        x : list of numpy.ndarray
            Total atoms concentrations
        Returns
        -------
        x : list of numpy.ndarray
            Updated total atoms concentrations
        """
        #Start by doing a geometrical search§
        x, root = self.bw_geom.search_for_keff(x, step_index)

        #check if upper geometrical limit gets hit
        if root >= self.bw_geom.bracket_limit[1]:
            # Reset geometry and refuel
            self.bw_geom._set_cell_attrib(self.restart_level)
            x,_ = self.bw_mat.search_for_keff(x, step_index)
            root = self.restart_level

        return x, root

class BatchwiseSchemeFlex():
    """
    Batchwise wrapper class, it wraps BatchwiseGeom and BatchwiseMat instances,
    with some user defined logic.

    This class should probably not be defined here, but we can keep it now for
    convenience

    The loop logic of this wrapper class is the following:

    1. Run BatchwiseGeom and return geometrical coefficient
    2. check if nuclide density of user defined nuclide is within limit
    3.1 if not, continue
    3.2 if yes, make a flex rotation and update all the materials before
    proceeding.
    4 When all flex rotations have been made, conclude the simulation.

    In this case if the bracket upper geometrical limit is hitted,
    simply stop the simulation.

    An instance of this class can be passed directly to an instance of the
    integrator class, such as :class:`openmc.deplete.CECMIntegrator`.

    Parameters
    ----------
    bw_list : list of Batchwise classes
        list of openmc.deplete.batchwise objects
    nuclide : str
        Nuclide to monitor nuclide density
    limit : float
        Nuclide density of nuclide to monitor
    flex_fraction : float
        Initial flex fraction in core
    span : float
        Flex rotation angle
    Attributes
    ----------
    rotation : int
        Initial flex rotation
        Default to 1.
    """
    def __init__(self, bw_list, n_timesteps, nuclide, limit, flex_fraction, span,
                    dilute_interval=None, restart_level = 0):

        if isinstance(bw_list, list):
            for bw in bw_list:
                if isinstance(bw, BatchwiseCell):
                    if bw.attrib_name == 'translation':
                        self.bw_geom_trans = bw
                    elif bw.attrib_name == 'rotation':
                        self.bw_geom_rot = bw
                elif isinstance(bw, BatchwiseMaterial):
                    if bw.attrib_name == 'addition':
                        self.bw_mat_add = bw
                    elif bw.attrib_name == 'dilute':
                        self.bw_mat_dil = bw
                else:
                    raise ValueError(f'{bw} is not a valid instance of'
                                      ' Batchwise class')
        else:
            raise ValueError(f'{bw_list} is not a list')

        self.bw_list = bw_list
        self.n_timesteps = n_timesteps
        check_value("nuclide exists", nuclide, \
                    self.bw_geom_trans.operator.nuclides_with_data)
        self.nuclide = nuclide
        self.limit = limit
        self.span = span
        self.flex_fraction = flex_fraction
        # start roation counter at 1
        self.rotation = 1

        #container for water level
        self.levels = []
        #Linear fit water level 1st coeff (line slope, or 1st derivative )
        self.a = 0

        # Dilution interval after all sections are converted
        self.dilute_interval = dilute_interval
        self.restart_level = restart_level
        #
        self.converged_flex = False

    def get_root(self):
        return self.bw_geom_trans.get_root()

    def set_density_function(self, mats, density_func, oxidation_states):
        for bw in self.bw_list:
            bw.set_density_function(mats, density_func, oxidation_states)

    def _nuclide_density(self, mat_name, nuc):
        """
        Return nuclide of specific nuclide in a material
        Parameters
        ----------
        mat_name : str
            Name of material to look for nuclide
        nuclide : str
            Name of nuclide to look for
        Returns
        -------
        atoms density : float
            Atom density in [atoms/b-cm]

        """
        mat_id = [id for id,mat in openmc.lib.materials.items() \
                    if mat.name == mat_name][0]

        if nuc not in openmc.lib.materials[mat_id].nuclides:
            return 0
        else:
            nuc_index = openmc.lib.materials[mat_id].nuclides.index(nuc)

        return openmc.lib.materials[mat_id].densities[nuc_index]

    def search_for_keff(self, x, step_index):
        """
        Perform the criticality search on the parametric material model.
        Will set the root of the `search_for_keff` function to the atoms
        concentrations vector.
        Parameters
        ----------
        x : list of numpy.ndarray
            Total atoms concentrations
        Returns
        -------
        x : list of numpy.ndarray
            Updated total atoms concentrations
        """

        nuc_density = self._nuclide_density('flex', self.nuclide)
        print('{} density: {:.7f} [atoms/b-cm]'.format(self.nuclide, nuc_density))
        #Compute liner fit of water level variation (needs at least 20 points)
        if len(self.levels) >= 8:
            self.a = np.polyfit(np.arange(0, len(self.levels), 1), \
                                np.array(self.levels), 1)[0]

        #If water level too close to the upper border don't apply rotation (just temporary ) or
        # water level first derivative is negative:
        if (nuc_density >= self.limit and \
            self.bw_geom_trans._get_cell_attrib() < self.bw_geom_trans.bracket_limit[1] - 15) or \
            self.a < 0:

            if self.rotation < 1/self.flex_fraction:
                print(f'Flex rotation nr: {self.rotation}')
                print(f'Slope coeff: {self.a}')
                # each time make a rotation anti-clockwise, i.e increase the
                #volume of flex
                self.bw_geom_rot._set_cell_attrib(self.span * \
                                            self.flex_fraction * self.rotation)
                # For every rotation, the volume and the materials must be updated
                x = self.bw_mat_add.add_materials(x, self.bw_geom_rot.cell_materials)
                self.rotation += 1
                # reset water level container
                self.levels = []
                self.a = 0
                if self.rotation == 1/self.flex_fraction:
                    print('All flex sections converted')
                    # set start dilution time step index from next one
                    self.start_dilute_index = step_index + 1
                    self.converged_flex = True
            else:
                if self.dilute_interval is not None:
                    if step_index in np.arange(self.start_dilute_index,
                                    self.n_timesteps,self.dilute_interval) and \
                    self.bw_geom_trans._get_cell_attrib() <= self.restart_level:
                        # restart level and perform dilution
                        self.bw_geom_trans._set_cell_attrib(self.restart_level)
                        x,_ = self.bw_mat_dil.search_for_keff(x, step_index)
                    else:
                        x, root = self.bw_geom_trans.search_for_keff(x, step_index)
                        # in this case if upper limit gets hit, stop directly
                        if root >= self.bw_geom_trans.bracket_limit[1]:
                            from pathlib import Path
                            print(f""" Reached maximum of
                                    {self.bw_geom_trans.bracket_limit[1]} cm
                                    exit..""")
                            Path('sim.done').touch()
                            for rank in range(comm.size):
                                comm.Abort()
                else:
                    print('Stop simulation')
                    from pathlib import Path
                    Path('sim.done').touch()
                    for rank in range(comm.size):
                        comm.Abort()

        x, root = self.bw_geom_trans.search_for_keff(x, step_index)
        # update water level container after geom keff search
        self.levels.append(root)
        # in this case if upper limit gets hit, stop directly
        if root >= self.bw_geom_trans.bracket_limit[1]:
            from pathlib import Path
            print(f'Reached maximum of {self.bw_geom_trans.bracket_limit[1]} cm'
                   ' exit..')
            Path('sim.done').touch()

            for rank in range(comm.size):
                comm.Abort()

        return x, root
