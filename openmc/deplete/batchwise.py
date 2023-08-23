from abc import ABC, abstractmethod
from collections.abc import Iterable
from openmc.search import _SCALAR_BRACKETED_METHODS, search_for_keff
from openmc import Materials, Material, Cell
from openmc.data import atomic_mass, AVOGADRO, ELEMENT_SYMBOL

import numpy as np
import os
import h5py

class Batchwise(ABC):

    def __init__(self, operator, model, bracket, bracket_limit,
                 bracketed_method='brentq', tol=0.01, target=1.0,
                 print_iterations=True, search_for_keff_output=True):

        self.operator = operator
        self.burn_mats = operator.burnable_mats
        self.local_mats = operator.local_mats
        self.model = model
        self.geometry = model.geometry

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

    def _get_cell_id(self, val):
        """Helper method for getting cell id from cell instance or cell name.
        Parameters
        ----------
        val : Openmc.Cell or str or int representing Cell
        Returns
        -------
        id : str
            Cell id
        """
        if isinstance(val, Cell):
            check_value('Cell id', val.id, [cell.id for cell in \
                                self.geometry.get_all_cells().values()])
            val = val.id

        elif isinstance(val, str):
            if val.isnumeric():
                check_value('Cell id', val, [str(cell.id) for cell in \
                                self.geometry.get_all_cells().values()])
                val = int(val)
            else:
                check_value('Cell name', val, [cell.name for cell in \
                                self.geometry.get_all_cells().values()])

                val = [cell.id for cell in \
                    self.geometry.get_all_cells().values() \
                    if cell.name == val][0]

        elif isinstance(val, int):
            check_value('Cell id', val, [cell.id for cell in \
                                self.geometry.get_all_cells().values()])

        else:
            ValueError(f'Cell: {val} is not recognized')

        return val

    def _search_for_keff(self, val):
        """
        Perform the criticality search for a given parametric model.
        It supports geometrical or material based `search_for_keff`.
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

        while res == None:
            search = search_for_keff(self._model_builder,
                    bracket = [bracket[0] + val, bracket[1] + val],
                    tol = tol,
                    bracketed_method = self.bracketed_method,
                    target = self.target,
                    print_iterations = self.print_iterations,
                    run_args = {'output': self.search_for_keff_output})

            # if len(search) is 3 search_for_keff was successful
            if len(search) == 3:
                res,_,_ = search

                #Check if root is within bracket limits
                if self.bracket_limit[0] < res <  self.bracket_limit[1]:
                    root = res

                else:
                    # Set res with the closest limit and continue
                    arg_min = abs(np.array(self.bracket_limit) - res).argmin()
                    warn('WARNING: Search_for_keff returned root out of '\
                         'bracket limit. Set root to {:.2f} and continue.'
                         .format(self.bracket_limit[arg_min]))
                    root = self.bracket_limit[arg_min]

            elif len(search) == 2:
                guesses, keffs = search

                #Check if all guesses are within bracket limits
                if all(self.bracket_limit[0] < guess < self.bracket_limit[1] \
                    for guess in guesses):
                    #Simple method to iteratively adapt the bracket
                    print('INFO: Function returned values below or above ' \
                           'target. Adapt bracket...')

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
                        if guess[np.argmax(keffs)] > guess[np.argmin(keffs)]:
                            dir = 1
                        else:
                            dir = -1
                        bracket[np.argmin(keffs)] = bracket[np.argmax(keffs)]
                        bracket[np.argmax(keffs)] += grad * (self.target - \
                                                  max(keffs).n) * dir
                    else:
                        if guess[np.argmax(keffs)] > guess[np.argmin(keffs)]:
                            dir = -1
                        else:
                            dir = 1
                        bracket[np.argmax(keffs)] = bracket[np.argmin(keffs)]
                        bracket[np.argmin(keffs)] += grad * (min(keffs).n - \
                                                  self.target) * dir

                else:
                    # Set res with closest limit and continue
                    arg_min = abs(np.array(self.bracket_limit) - guesses).argmin()
                    warn('WARNING: Adaptive iterative bracket went off '\
                         'bracket limits. Set root to {:.2f} and continue.'
                         .format(self.bracket_limit[arg_min]))
                    root = self.bracket_limit[arg_min]

            else:
                raise ValueError(f'ERROR: Search_for_keff output is not valid')

        return root

    def _save_res(self, type, step_index, root):
        """
        Save results to msr_results.h5 file.
        Parameters
        ----------
        type : str
            String to characterize geometry and material results
        step_index : int
            depletion time step index
        root : float or dict
             Root of the search_for_keff function
        """
        filename = 'msr_results.h5'
        kwargs = {'mode': "a" if os.path.isfile(filename) else "w"}

        if comm.rank == 0:
            with h5py.File(filename, **kwargs) as h5:
                name = '_'.join([type, str(step_index)])
                if name in list(h5.keys()):
                    last = sorted([int(re.split('_',i)[1]) for i in h5.keys()])[-1]
                    step_index = last + 1
                h5.create_dataset('_'.join([type, str(step_index)]), data=root)

    def _update_volumes_after_depletion(self, x):
        """
        After a depletion step, both materials volume and density change, due to
        decay, transmutation reactions and transfer rates, if set.
        At present we lack an implementation to calculate density and volume
        changes due to the different molecules speciation. Therefore, the assumption
        is to consider the density constant and let the material volume
        vary with the change in nuclides concentrations.
        The method uses the nuclides concentrations coming from the previous Bateman
        solution and calculates a new volume, keeping the mass density of the material
        constant. It will then assign the volumes to the AtomNumber instance.

        Parameters
        ----------
        x : list of numpy.ndarray
            Total atom concentrations
        """
        self.operator.number.set_density(x)

        for rank in range(comm.size):
            number_i = comm.bcast(self.operator.number, root=rank)

            for i, mat in enumerate(number_i.materials):
                # Total nuclides density
                density = 0
                for nuc in number_i.nuclides:
                    # total number of atoms
                    val = number_i[mat, nuc]
                    # obtain nuclide density in atoms-g/mol
                    density +=  val * atomic_mass(nuc)
                # Get mass dens from beginning, intended to be held constant
                rho = openmc.lib.materials[int(mat)].get_density('g/cm3')
                number_i.volume[i] = density / AVOGADRO / rho

class BatchwiseCell(Batchwise):

    def __init__(self, operator, model, cell, attrib_name, axis, bracket,
                 bracket_limit, bracketed_method='brentq', tol=0.01, target=1.0,
                 print_iterations=True, search_for_keff_output=True):

        super().__init__(operator, model, bracket, bracket_limit,
                         bracketed_method, tol, target, print_iterations,
                         search_for_keff_output)

        self.cell_id = super()._get_cell_id(cell_id_or_name)
        check_value('attrib_name', attrib_name,
                    ('rotation', 'translation'))
        self.attrib_name = attrib_name

        #index of cell directionnal axis
        check_value('axis', axis, (0,1,2))
        self.axis = axis

        # Initialize vector
        self.vector = np.zeros(3)

        # materials that fill the attribute cell, if depletables
        self.cell_material_ids = [cell.fill.id for cell in \
            self.geometry.get_all_cells()[self.cell_id].fill.cells.values() \
            if cell.fill.depletable]

    def _get_cell_attrib(self):
        """
        Get cell attribute coefficient.
        Returns
        -------
        coeff : float
            cell coefficient
        """
        for cell in openmc.lib.cells.values():
            if cell.id == self.cell_id:
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
        var : float
            Surface coefficient to set
        geometry : openmc.model.geometry
            OpenMC geometry model
        attrib_name : str
            Currently only translation is implemented
        """
        self.vector[self.axis] = val
        for cell in openmc.lib.cells.values():
            if cell.id == self.cell_id or cell.name == self.cell_id:
                setattr(cell, self.attrib_name, self.vector)

    def _update_materials(self, x):
        """
        Assign concentration vectors from Bateman solution at previous
        timestep to the in-memory model materials, after having recalculated the
        material volume.

        Parameters
        ----------
        x : list of numpy.ndarray
            Total atom concentrations
        """
        super()._update_volumes_after_depletion(x)

        for rank in range(comm.size):
            number_i = comm.bcast(self.operator.number, root=rank)

            for mat in number_i.materials:
                nuclides = []
                densities = []

                for nuc in number_i.nuclides:
                    # get atom density in atoms/b-cm
                    val = 1.0e-24 * number_i.get_atom_density(mat, nuc)
                    if nuc in self.operator.nuclides_with_data:
                        if val > self.atom_density_limit:
                            nuclides.append(nuc)
                            densities.append(val)

                #set nuclide densities to model in memory (C-API)
                openmc.lib.materials[int(mat)].set_densities(nuclides, densities)

    def _update_volumes(self):
        openmc.lib.calculate_volumes()
        res = openmc.VolumeCalculation.from_hdf5('volume_1.h5')

        number_i = self.operator.number
        for mat_idx, mat_id in enumerate(self.local_mats):
            if mat_id in self.cell_material_ids:
                number_i.volume[mat_idx] = res.volumes[int(mat_id)].n

    def _update_x(self):
        number_i = self.operator.number
        for mat_idx, mat_id in enumerate(self.local_mats):
            if mat_id in self.cell_material_ids:
                for nuc_idx, nuc in enumerate(number_i.burnable_nuclides):
                    x[mat_idx][nuc_idx] = number_i.volume[mat_idx] * \
                            number_i.get_atom_density(mat_idx, nuc)
        return x

    def _model_builder(self, param):
        """
        Builds the parametric model that is passed to the `msr_search_for_keff`
        function by setting the parametric variable to the geoemetrical cell.
        Parameters
        ----------
        param : model parametricl variable
            for examlple: cell translation coefficient
        Returns
        -------
        _model :  openmc.model.Model
            OpenMC parametric model
        """
        self._set_cell_attrib(param)
        #Calulate new volume and update if materials filling the cell are
        # depletable
        if self.cell_material_ids:
            self._update_volumes()
        return self.model

    def search_for_keff(self, x, step_index):
        """
        Perform the criticality search on the parametric geometrical model.
        Will set the root of the `search_for_keff` function to the cell
        attribute.
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

        # Update volume and concentration vectors before performing the search_for_keff
        self._update_materials(x)

        # Calculate new cell attribute
        root = super().search_for_keff(val)

        # set results value as attribute in the geometry
        self._set_cell_attrib(root)
        print('UPDATE: old value: {:.2f} cm --> ' \
              'new value: {:.2f} cm'.format(val, root))

        # x needs to be updated after the search with new volumes if materials cell
        # depletable
        if self.cell_material_ids:
            self._update_x()

        #Store results
        super()._save_res('geometry', step_index, root)
        return x
