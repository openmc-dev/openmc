from abc import ABC, abstractmethod
from copy import deepcopy
from numbers import Real
from warnings import warn
import numpy as np

import openmc.lib
from openmc.mpi import comm
from openmc.search import _SCALAR_BRACKETED_METHODS, search_for_keff
from openmc import Material, Cell
from openmc.data import atomic_mass, AVOGADRO
from openmc.checkvalue import (check_type, check_value, check_less_than,
    check_iterable_type, check_length)

class Batchwise(ABC):
    """Abstract class defining a generalized batchwise scheme.

    In between transport and depletion steps batchwise operations can be added
    with the aim of maintain a system critical. These operations act on OpenMC
    Cell or Material, parametrizing one or more coefficients while running a
    search_for_kef algorithm.

    Specific classes for running batchwise depletion calculations are
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
    burn_mats : list of str
        List of burnable materials ids
    local_mats : list of str
        All burnable material IDs being managed by a single process
        List of burnable materials ids
    """
    def __init__(self, operator, model, bracket, bracket_limit,
                 density_treatment = 'constant-volume', bracketed_method='brentq',
                 tol=0.01, target=1.0, print_iterations=True,
                 search_for_keff_output=True):

        self.operator = operator
        self.burn_mats = operator.burnable_mats
        self.local_mats = operator.local_mats
        self.model = model
        self.geometry = model.geometry

        check_value('density_treatment', density_treatment,
                                        ('constant-density','constant-volume'))
        self.density_treatment = density_treatment

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

    @abstractmethod
    def _model_builder(self, param):
        """
        Builds the parametric model to be passed to the `search_for_keff`
        algorithm.
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
        It supports cell or material based `search_for_keff`.
        If the solution lies off the initial bracket, the method will iteratively
        adapt the bracket to be able to run  `search_for_keff` effectively.
        The ratio between the bracket and the returned keffs values will be
        the scaling factor to iteratively adapt the brackt unitl a valid solution
        is found.
        If the adapted bracket is moved too far off, i.e. above or below user
        inputted bracket limits interval, the algorithm will stop and the closest
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

        # Run until a search_for_keff root is found or ouf ot limits
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

    def _update_volumes(self):
        """
        Update number volume.
        After a depletion step, both material volume and density change, due to
        decay, transmutation reactions and transfer rates, if set.
        At present we lack an implementation to calculate density and volume
        changes due to the different molecules speciation. Therefore, OpenMC
        assumes by default that the depletable volume does not change and only
        updates the nuclide densities and consequently the total material density.
        This method, vice-versa, assumes that the total material density does not
        change and updates the materials volume in AtomNumber dataset.

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

        if self.density_treatment == 'constant-density':
            self._update_volumes()

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
        Update x with volumes, before assign them to AtomNumber materials.
        Parameters
        ----------
        x : list of numpy.ndarray
            Total atoms concentrations
        volumes : dict

        Returns
        -------
        x : list of numpy.ndarray
            Updated total atoms concentrations
        """
        number_i = self.operator.number

        for mat_idx, mat in enumerate(self.local_mats):

            if mat in volumes:
                res_vol = volumes[mat]
                # Normalize burbable nuclides in x vector without cross section data
                for nuc_idx, nuc in enumerate(number_i.burnable_nuclides):
                    if nuc not in self.operator.nuclides_with_data:
                        # normalzie with new volume
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
        Builds the parametric model that is passed to the `search_for_keff`
        function by setting the parametric variable to the cell.
        Parameters
        ----------
        param : model parametricl variable
            for examlple: cell translation coefficient
        Returns
        -------
        self.model :  openmc.model.Model
            OpenMC parametric model
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
        Perform the criticality search parametrizing the cell coefficient.
        The `search_for_keff` solution is then set as the new cell attribute
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
    Batchwise cell-based with geometrical-attribute class.

    A user doesn't need to call this class directly.
    Instead an instance of this class is automatically created by calling
    :meth:`openmc.deplete.Integrator.add_batchwise` method from an integrator
    class, such as  :class:`openmc.deplete.CECMIntegrator`.

    .. versionadded:: 0.13.4

    Parameters
    ----------
    cell : openmc.Cell or int or str
        OpenMC Cell identifier to where apply batchwise scheme
    attrib_name : str {'translation', 'rotation'}
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
    attrib_name : str {'translation', 'rotation'}
        Cell attribute type
    axis : int {0,1,2}
        Directional axis for geometrical parametrization, where 0, 1 and 2 stand
        for 'x', 'y' and 'z', respectively.
    """
    def __init__(self, cell, attrib_name, operator, model, axis, bracket,
                 bracket_limit, density_treatment='constant-volume',
                 bracketed_method='brentq', tol=0.01, target=1.0,
                 print_iterations=True, search_for_keff_output=True):

        super().__init__(cell, operator, model, bracket, bracket_limit,
                         density_treatment, bracketed_method, tol, target,
                         print_iterations, search_for_keff_output)

        check_value('attrib_name', attrib_name,
                    ('rotation', 'translation'))
        self.attrib_name = attrib_name

        # check if cell is filled with 2 cells
        check_length('fill materials', self.cell.fill.cells, 2)
        #index of cell directionnal axis
        check_value('axis', axis, (0,1,2))
        self.axis = axis

        # Initialize vector
        self.vector = np.zeros(3)

        self.cell_materials = [cell.fill for cell in \
                        self.cell.fill.cells.values() if cell.fill.depletable]

        if self.cell_materials:
            self._initialize_volume_calc()

    @classmethod
    def from_params(cls, obj, attr, operator, model, **kwargs):
        return cls(obj, attr, operator, model, **kwargs)

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
        Attributes are only applied to cells filled with a universe
        Parameters
        ----------
        val : float
            Cell coefficient to set in cm for translation and def for rotation
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
        mat_vol = openmc.VolumeCalculation(self.cell_materials, 100000, ll, ur)
        self.model.settings.volume_calculations = mat_vol

    def _calculate_volumes(self):
        """
        Perform stochastic volume calculation
        Returns
        -------
        volumes : dict
            Dictionary of calculate volumes, where key is the mat id
        """
        openmc.lib.calculate_volumes()
        volumes = {}
        res = openmc.VolumeCalculation.from_hdf5('volume_1.h5')
        for mat_idx, mat in enumerate(self.local_mats):
            if int(mat) in [mat.id for mat in self.cell_materials]:
                volumes[mat] = res.volumes[int(mat)].n
        return volumes

class BatchwiseCellTemperature(BatchwiseCell):
    """
    Batchwise cell-based with temperature-attribute class.

    A user doesn't need to call this class directly.
    Instead an instance of this class is automatically created by calling
    :meth:`openmc.deplete.Integrator.add_batchwise` method from an integrator
    class, such as  :class:`openmc.deplete.CECMIntegrator`.

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
    """
    def __init__(self, cell, operator, model, bracket, bracket_limit,
                 density_treatment='constant-volume', bracketed_method='brentq',
                 tol=0.01, target=1.0, print_iterations=True,
                 search_for_keff_output=True):

        super().__init__(cell, operator, model, bracket, bracket_limit,
                         density_treatment, bracketed_method, tol, target,
                         print_iterations, search_for_keff_output)

        # check if initial temperature has been set to right cell material
        if isinstance(self.cell.fill, openmc.Universe):
            cells = [cell for cell in self.cell.fill.cells.values()  \
                         if cell.fill.temperature]
            check_length('Only one cell with temperature',cells,1)
            self.cell = cells[0]

        check_type('temperature cell real', self.cell.fill.temperature, Real)

    @classmethod
    def from_params(cls, obj, attr, operator, model, **kwargs):
        return cls(obj, operator, model, **kwargs)

    def _get_cell_attrib(self):
        """
        Get cell attribute coefficient.
        Returns
        -------
        coeff : float
            cell temperature in Kelvin
        """
        for cell in openmc.lib.cells.values():
            if cell.id == self.cell.id:
                return cell.get_temperature()

    def _set_cell_attrib(self, val):
        """
        Set cell attribute to the cell instance.
        Parameters
        ----------
        val : float
            Cell temperature to set in Kelvin
        """
        for cell in openmc.lib.cells.values():
            if cell.id == self.cell.id:
                cell.set_temperature(val)

class BatchwiseMaterial(Batchwise):
    """Abstract class holding batchwise material-based functions.

    Specific classes for running batchwise depletion calculations are
    implemented as derived class of BatchwiseMaterial.

    .. versionadded:: 0.13.4

    Parameters
    ----------
    material : openmc.Material or int or str
        OpenMC Material identifier to where apply batchwise scheme
    operator : openmc.deplete.Operator
        OpenMC operator object
    model : openmc.model.Model
        OpenMC model object
    mat_vector : dict
        Dictionary of material composition to paremetrize, where a pair key, value
        represents a nuclide and its weight fraction.
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
    material : openmc.Material or int or str
        OpenMC Material identifier to where apply batchwise scheme
    mat_vector : dict
        Dictionary of material composition to paremetrize, where a pair key, value
        represents a nuclide and its weight fraction.
    """

    def __init__(self, material, operator, model, mat_vector, bracket,
                 bracket_limit, density_treatment='constant-volume',
                 bracketed_method='brentq', tol=0.01, target=1.0,
                 print_iterations=True, search_for_keff_output=True):

        super().__init__(operator, model, bracket, bracket_limit,
                         density_treatment, bracketed_method, tol, target,
                         print_iterations, search_for_keff_output)

        self.material = self._get_material(material)

        check_type("material vector", mat_vector, dict, str)
        for nuc in mat_vector.keys():
            check_value("check nuclide exists", nuc, self.operator.nuclides_with_data)

        if round(sum(mat_vector.values()), 2) != 1.0:
            # Normalize material elements vector
            sum_values = sum(mat_vector.values())
            for elm in mat_vector:
                mat_vector[elm] /= sum_values
        self.mat_vector = mat_vector

    @classmethod
    def from_params(cls, obj, attr, operator, model, **kwargs):
        return cls(obj, operator, model, **kwargs)

    def _get_material(self, val):
        """Helper method for getting openmc material from Material instance or
        material name or id.
        Parameters
        ----------
        val : Openmc.Material or str or int representing material name/id
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

    @abstractmethod
    def _model_builder(self, param):
        """
        Builds the parametric model that is passed to the `msr_search_for_keff`
        function by updating the material densities and setting the parametric
        variable as function of the nuclides vector. Since this is a paramteric
        material addition (or removal), we can parametrize the volume as well.
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
    Batchwise material-based class for refueling (addition or removal) scheme.

    A user doesn't need to call this class directly.
    Instead an instance of this class is automatically created by calling
    :meth:`openmc.deplete.Integrator.add_batchwise` method from an integrator
    class, such as  :class:`openmc.deplete.CECMIntegrator`.

    .. versionadded:: 0.13.4

    Parameters
    ----------
    material : openmc.Material or int or str
        OpenMC Material identifier to where apply batchwise scheme
    operator : openmc.deplete.Operator
        OpenMC operator object
    model : openmc.model.Model
        OpenMC model object
    mat_vector : dict
        Dictionary of material composition to paremetrize, where a pair key, value
        represents a nuclide and its weight fraction.
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
    material : openmc.Material or int or str
        OpenMC Material identifier to where apply batchwise scheme
    mat_vector : dict
        Dictionary of material composition to paremetrize, where a pair key, value
        represents a nuclide and its weight fraction.
    """

    def __init__(self, material, operator, model, mat_vector, bracket,
                 bracket_limit, density_treatment='constant-volume',
                 bracketed_method='brentq', tol=0.01, target=1.0,
                 print_iterations=True, search_for_keff_output=True):

        super().__init__(material, operator, model, mat_vector, bracket,
                         bracket_limit, density_treatment, bracketed_method, tol,
                         target, print_iterations, search_for_keff_output)

    @classmethod
    def from_params(cls, obj, attr, operator, model, **kwargs):
        return cls(obj, operator, model, **kwargs)

    def _model_builder(self, param):
        """
        Builds the parametric model that is passed to the `msr_search_for_keff`
        function by updating the material densities and setting the parametric
        variable to the material nuclides to add. We can either fix the total
        volume or the material density according to user input.
        Parameters
        ----------
        param :
            Model function variable, total grams of material to add or remove
        Returns
        -------
        model :  openmc.model.Model
            OpenMC parametric model
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

        # alwyas need to return a model
        return self.model

    def _calculate_volumes(self, res):
        """
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
