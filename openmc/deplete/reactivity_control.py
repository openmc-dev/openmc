from abc import ABC, abstractmethod
from copy import deepcopy
from numbers import Real
from warnings import warn
import numpy as np
from math import isclose

import openmc.lib
from openmc.mpi import comm
from openmc.search import _SCALAR_BRACKETED_METHODS, search_for_keff
from openmc import Material, Cell
from openmc.data import atomic_mass, AVOGADRO
from openmc.checkvalue import (
    check_type,
    check_value,
    check_less_than,
    check_greater_than,
    check_iterable_type,
    check_length,
)


class ReactivityController(ABC):
    """Abstract class defining a generalized reactivity control.

    Reactivity control schemes, such as control rod adjustment or material
    refuelling to control reactivity and maintain keff constant and equal to a
    desired value, usually one.

    A reactivity control scheme can be added here to an integrator instance,
    such as  :class:`openmc.deplete.CECMIntegrator`, to parametrize one system
    variable with the aim of satisfy certain design criteria, such as keeping
    keff equal to one, while running transport-depletion calculations.

    This abstract class sets the requirements
    for the reactivity control set. Users should instantiate
    :class:`openmc.deplete.reactivity_control.CellReactivityController`
    or
    :class:`openmc.deplete.reactivity_control.MaterialReactivityController`
    rather than this class.
    .. versionadded:: 0.14.1

    Parameters
    ----------
    operator : openmc.deplete.Operator
        OpenMC operator object
    bracket : list of float
        Initial bracketing interval to search for the solution, relative to the
        solution at previous step.
    bracket_limit : list of float
        Absolute bracketing interval lower and upper; if during the adaptive
        algorithm the search_for_keff solution lies off these limits the closest
        limit will be set as new result.
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
        Print a status message each iteration
        Default to True
    search_for_keff_output : Bool, Optional
        Print full transport logs during iterations (e.g., inactive generations, active
        generations). Default to False

    Attributes
    ----------
    burn_mats : list of str
        List of burnable materials ids
    local_mats : list of str
        All burnable material IDs being managed by a single process
        List of burnable materials ids

    """

    def __init__(
        self,
        operator,
        bracket,
        bracket_limit,
        bracketed_method="brentq",
        tol=0.01,
        target=1.0,
        print_iterations=True,
        search_for_keff_output=True,
    ):

        self.operator = operator
        self.burn_mats = operator.burnable_mats
        self.local_mats = operator.local_mats
        self.model = operator.model
        self.geometry = operator.model.geometry
        self.bracket = bracket

        check_iterable_type("bracket_limit", bracket_limit, Real)
        check_length("bracket_limit", bracket_limit, 2)
        check_less_than("bracket limit values", bracket_limit[0], bracket_limit[1])

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
        check_value("bracketed_method", value, _SCALAR_BRACKETED_METHODS)
        if value != "brentq":
            warn("""brentq bracketed method is recommended to ensure
                    convergence of the method""")
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
    def from_params(cls, obj, operator, **kwargs):
        return cls(obj, operator, **kwargs)

    @abstractmethod
    def _model_builder(self, param):
        """Build the parametric model to be solved.

        Callable function which builds a model according to a passed
        parameter. This function must return an openmc.model.Model object.

        Parameters
        ----------
        param : parameter
            model function variable

        Returns
        -------
        openmc.model.Model
            OpenMC parametric model
        """

    @abstractmethod
    def search_for_keff(self, x, step_index):
        """Perform the criticality search on parametric material variable.

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
        root : float
             Search_for_keff returned root value
        """

    @abstractmethod
    def _adjust_volumes(self, root=None):
        """Adjust volumes after criticality search when cell or material are
        modified.

        Parameters
        ----------
        root : float, Optional
            :meth:`openmc.search.search_for_keff` volume root for
            :class:`openmc.deplete.reactivity_control.CellReactivityController`
            instances, in cm3
        Returns
        -------
        volumes : dict
            Dictionary of calculated volume, where key mat id and value
            material volume, in cm3
        """

    def _search_for_keff(self, val):
        """Perform the criticality search for a given parametric model.

        If the solution lies off the initial bracket, this method iteratively
        adapts it until :meth:`openmc.search.search_for_keff` return a valid
        solution.
        It calculates the ratio between guessed and corresponding keffs values
        as the proportional term, so that the bracket can be moved towards the
        target.
        A bracket limit poses the upper and lower boundaries to the adapting
        bracket. If one limit is hit, the algoritm will stop and the closest
        limit value will be set.

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
        # make sure we don't modify original bracket
        bracket = deepcopy(self.bracket)

        # Run until a search_for_keff root is found or out of limits
        root = None
        while root == None:
            search = search_for_keff(
                self._model_builder,
                bracket=np.array(bracket) + val,
                tol=self.tol,
                bracketed_method=self.bracketed_method,
                target=self.target,
                print_iterations=self.print_iterations,
                run_args={"output": self.search_for_keff_output},
                run_in_memory=True,
            )

            # if len(search) is 3 search_for_keff was successful
            if len(search) == 3:
                res, _, _ = search

                # Check if root is within bracket limits
                if self.bracket_limit[0] < res < self.bracket_limit[1]:
                    root = res

                else:
                    # Set res with the closest limit and continue
                    arg_min = abs(np.array(self.bracket_limit) - res).argmin()
                    warn(
                        "Search_for_keff returned root out of "
                        "bracket limit. Set root to {:.2f} and continue.".format(
                            self.bracket_limit[arg_min]
                        )
                    )
                    root = self.bracket_limit[arg_min]

            elif len(search) == 2:
                guesses, keffs = search

                # Check if all guesses are within bracket limits
                if all(
                    self.bracket_limit[0] <= guess <= self.bracket_limit[1]
                    for guess in guesses
                ):
                    # Simple method to iteratively adapt the bracket
                    warn(
                        "Search_for_keff returned values below or above "
                        "target. Trying to iteratively adapt bracket..."
                    )

                    # if the difference between keffs is smaller than 1 pcm,
                    # the slope calculation will be overshoot, so let's set the root
                    # to the closest guess value
                    if abs(np.diff(keffs)) < 1.0e-5:
                        arg_min = abs(self.target - np.array(keffs)).argmin()
                        warn(
                            "Difference between keff values is "
                            "smaller than 1 pcm. "
                            "Set root to guess value with "
                            "closest keff to target and continue..."
                        )
                        root = guesses[arg_min]

                    # Calculate slope as ratio of delta bracket and delta keffs
                    slope = abs(np.diff(bracket) / np.diff(keffs))[0].n
                    # Move the bracket closer to presumed keff root by using the
                    # slope

                    # Two cases: both keffs are below or above target
                    if np.mean(keffs) < self.target:
                        # direction of moving bracket: +1 is up, -1 is down
                        if guesses[np.argmax(keffs)] > guesses[np.argmin(keffs)]:
                            dir = 1
                        else:
                            dir = -1
                        bracket[np.argmin(keffs)] = bracket[np.argmax(keffs)]
                        bracket[np.argmax(keffs)] += (
                            slope * (self.target - max(keffs).n) * dir
                        )
                    else:
                        if guesses[np.argmax(keffs)] > guesses[np.argmin(keffs)]:
                            dir = -1
                        else:
                            dir = 1
                        bracket[np.argmax(keffs)] = bracket[np.argmin(keffs)]
                        bracket[np.argmin(keffs)] += (
                            slope * (min(keffs).n - self.target) * dir
                        )

                    # check if adapted bracket lies completely outside of limits
                    msg = (
                        "WARNING: Adaptive iterative bracket {} went off "
                        "bracket limits. Set root to {:.2f} and continue."
                    )
                    if all(np.array(bracket) + val <= self.bracket_limit[0]):
                        warn(msg.format(bracket, self.bracket_limit[0]))
                        root = self.bracket_limit[0]

                    if all(np.array(bracket) + val >= self.bracket_limit[1]):
                        warn(msg.format(bracket, self.bracket_limit[1]))
                        root = self.bracket_limit[1]

                    # check if adapted bracket ends are outside bracketing limits
                    if bracket[1] + val > self.bracket_limit[1]:
                        bracket[1] = self.bracket_limit[1] - val

                    if bracket[0] + val < self.bracket_limit[0]:
                        bracket[0] = self.bracket_limit[0] - val

                else:
                    # Set res with closest limit and continue
                    arg_min = abs(np.array(self.bracket_limit) - guesses).argmin()
                    warn(
                        "Search_for_keff returned values off "
                        "bracket limits. Set root to {:.2f} and continue.".format(
                            self.bracket_limit[arg_min]
                        )
                    )
                    root = self.bracket_limit[arg_min]

            else:
                raise ValueError("ERROR: Search_for_keff output is not valid")

        return root

    def _update_materials(self, x):
        """Update number density and material compositions in OpenMC on all
        processes.

        Parameters
        ----------
        x : list of numpy.ndarray
            Total atom concentrations
        """
        self.operator.number.set_density(x)

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
        """Update x vector with new volumes, before assign them to AtomNumber.

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
                        x[mat_idx][nuc_idx] *= number_i.get_mat_volume(mat) / res_vol

                # update all concentration data with the new updated volumes
                for nuc, dens in zip(
                    openmc.lib.materials[int(mat)].nuclides,
                    openmc.lib.materials[int(mat)].densities,
                ):
                    nuc_idx = number_i.index_nuc[nuc]
                    n_of_atoms = dens / 1.0e-24 * res_vol

                    if nuc in number_i.burnable_nuclides:
                        # convert [#atoms/b-cm] into [#atoms]
                        x[mat_idx][nuc_idx] = n_of_atoms
                    # when the nuclide is not in depletion chain update the AtomNumber
                    else:
                        # Atom density needs to be in [#atoms/cm3]
                        number_i[mat, nuc] = n_of_atoms

                # Now we can update the new volume in AtomNumber
                number_i.volume[mat_idx] = res_vol
        return x

class CellReactivityController(ReactivityController):
    """Abstract class holding reactivity control cell-based functions.

    .. versionadded:: 0.14.1

    Parameters
    ----------
    cell : openmc.Cell or int or str
        OpenMC Cell identifier to where apply a reactivty control
    operator : openmc.deplete.Operator
        OpenMC operator object
    attribute : str
        openmc.lib.cell attribute. Only support 'translation' and 'rotation'
    bracket : list of float
        Initial bracketing interval to search for the solution, relative to the
        solution at previous step.
    bracket_limit : list of float
        Absolute bracketing interval lower and upper; if during the adaptive
        algorithm the search_for_keff solution lies off these limits the closest
        limit will be set as new result.
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
    attribute : str
        openmc.lib.cell attribute. Only support 'translation' and 'rotation'
    depletable_cells : list of openmc.Cell
        Cells that fill the openmc.Universe that fills the main cell
        to where apply batch wise scheme, if cell materials are set as
        depletable.
    lib_cell : openmc.lib.Cell
        Corresponding openmc.lib.Cell once openmc.lib is initialized
    axis : int {0,1,2}
        Directional axis for geometrical parametrization, where 0, 1 and 2 stand
        for 'x', 'y' and 'z', respectively.

    """

    def __init__(
        self,
        cell,
        operator,
        attribute,
        bracket,
        bracket_limit,
        axis,
        bracketed_method="brentq",
        tol=0.01,
        target=1.0,
        samples=1000000,
        print_iterations=True,
        search_for_keff_output=True,
    ):

        super().__init__(
            operator,
            bracket,
            bracket_limit,
            bracketed_method,
            tol,
            target,
            print_iterations,
            search_for_keff_output,
        )

        self.cell = self._get_cell(cell)

        # Only lib cell attributes are valid
        check_value("attribute", attribute, ('translation', 'rotation'))
        self.attribute = attribute

        # Initialize to None, as openmc.lib is not initialized yet here
        self.lib_cell = None

        # A cell translation or roation must correspond to a geometrical
        # variation of its filled materials, thus at least 2 are necessary
        if isinstance(self.cell.fill, openmc.Universe):
            # check if cell fill universe contains at least 2 cells
            check_greater_than("universe cells",
                len(self.cell.fill.cells), 2, True)

        if isinstance(self.cell.fill, openmc.Lattice):
            # check if cell fill lattice contains at least 2 universes
            check_greater_than("lattice universes",
                len(self.cell.fill.get_unique_universes()), 2, True)

        self.depletable_cells = [
                cell
                for cell in self.cell.fill.cells.values()
                if cell.fill.depletable
            ]

        # index of cell directional axis
        check_value("axis", axis, (0, 1, 2))
        self.axis = axis

        check_type("samples", samples, int)
        self.samples = samples

        if self.depletable_cells:
            self._initialize_volume_calc()

    def _get_cell_attrib(self):
        """Get cell attribute coefficient.

        Returns
        -------
        coeff : float
            cell coefficient
        """
        return getattr(self.lib_cell, self.attribute)[self.axis]

    def _set_cell_attrib(self, val):
        """Set cell attribute to the cell instance.

        Parameters
        ----------
        val : float
            cell coefficient to set
        """
        attr = getattr(self.lib_cell, self.attribute)
        attr[self.axis] = val
        setattr(self.lib_cell, self.attribute, attr)

    def _set_lib_cell(self):
        """Set openmc.lib.cell cell to self.lib_cell attribute
        """
        self.lib_cell = [
            cell for cell in openmc.lib.cells.values() if cell.id == self.cell.id
        ][0]

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
        cell_bundles = [
            (cell.id, cell.name, cell)
            for cell in self.geometry.get_all_cells().values()
        ]

        if isinstance(val, Cell):
            check_value(
                "Cell exists", val, [cell_bundle[2] for cell_bundle in cell_bundles]
            )

        elif isinstance(val, str):
            if val.isnumeric():
                check_value(
                    "Cell id exists",
                    int(val),
                    [cell_bundle[0] for cell_bundle in cell_bundles],
                )

                val = [
                    cell_bundle[2]
                    for cell_bundle in cell_bundles
                    if cell_bundle[0] == int(val)
                ][0]

            else:
                check_value(
                    "Cell name exists",
                    val,
                    [cell_bundle[1] for cell_bundle in cell_bundles],
                )

                val = [
                    cell_bundle[2]
                    for cell_bundle in cell_bundles
                    if cell_bundle[1] == val
                ][0]

        elif isinstance(val, int):
            check_value(
                "Cell id exists", val, [cell_bundle[0] for cell_bundle in cell_bundles]
            )

            val = [
                cell_bundle[2] for cell_bundle in cell_bundles if cell_bundle[0] == val
            ][0]

        else:
            ValueError(f"Cell: {val} is not supported")

        return val

    def _initialize_volume_calc(self):
        """Set volume calculation model settings of depletable materials filling
        the parametric Cell.
        """
        ll, ur = self.geometry.bounding_box
        mat_vol = openmc.VolumeCalculation(self.depletable_cells, self.samples, ll, ur)
        self.model.settings.volume_calculations = mat_vol

    def _adjust_volumes(self):
        """Perform stochastic volume calculation and return new volumes.

        Returns
        -------
        volumes : dict
            Dictionary of calculated volumes, where key is mat id and value
            material volume, in cm3
        """
        openmc.lib.calculate_volumes()
        volumes = {}
        if comm.rank == 0:
            res = openmc.VolumeCalculation.from_hdf5("volume_1.h5")
            for cell in self.depletable_cells:
                mat_id = cell.fill.id
                volumes[str(mat_id)] = res.volumes[cell.id].n
        volumes = comm.bcast(volumes)
        return volumes

    def _model_builder(self, param):
        """Builds the parametric model to be passed to the
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

        return self.model

    def search_for_keff(self, x, step_index):
        """Perform the criticality search on the parametric cell coefficient and
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
        root : float
             Search_for_keff returned root value
        """
        # set _cell argument, once openmc.lib is initialized
        if self.lib_cell is None:
            self._set_lib_cell()

        # Get cell attribute from previous iteration
        val = self._get_cell_attrib()
        check_type("Cell coeff", val, Real)

        # Update all material densities from concentration vectors
        # before performing the search_for_keff. This is needed before running
        # the  transport equations in the search_for_keff algorithm
        self._update_materials(x)

        # Calculate new cell attribute
        root = self._search_for_keff(val)

        # set results value as attribute in the geometry
        self._set_cell_attrib(root)

        # if at least one of the cell fill materials is depletable, assign new
        # volume to the material and update x and number accordingly
        # new volume
        if self.depletable_cells:
            volumes = self._adjust_volumes()
            x = self._update_x_and_set_volumes(x, volumes)

        return x, root

class MaterialReactivityController(ReactivityController):
    """Abstract class holding reactivity control material-based functions.

    .. versionadded:: 0.14.1

    Parameters
    ----------
    material : openmc.Material or int or str
        OpenMC Material identifier to where apply reactivity control
    operator : openmc.deplete.Operator
        OpenMC operator object
    material_to_mix : openmc.Material
        OpenMC Material to mix with main material for reactivity control.
    bracket : list of float
        Initial bracketing interval to search for the solution, relative to the
        solution at previous step.
    bracket_limit : list of float
        Absolute bracketing interval lower and upper; if during the adaptive
        algorithm the search_for_keff solution lies off these limits the closest
        limit will be set as new result.
    units : {'grams', 'atoms', 'cc', 'cm3'} str
        Units of material parameter to be mixed.
        Default to 'cc'
    dilute : bool
        Whether or not to update material volume after a material mix.
        Default to False
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
    material_to_mix : openmc.Material
        OpenMC Material to mix with main material for reactivity control.
    units : {'grams', 'atoms', 'cc', 'cm3'} str
        Units of material parameter to be mixed.
        Default to 'cc'
    dilute : bool
        Whether or not to update material volume after a material mix.
        Default to False

    """

    def __init__(
        self,
        material,
        operator,
        material_to_mix,
        bracket,
        bracket_limit,
        units='cc',
        dilute=False,
        bracketed_method="brentq",
        tol=0.01,
        target=1.0,
        print_iterations=True,
        search_for_keff_output=True,
    ):

        super().__init__(
            operator,
            bracket,
            bracket_limit,
            bracketed_method,
            tol,
            target,
            print_iterations,
            search_for_keff_output,
        )

        self.material = self._get_material(material)

        check_type("material to mix", material_to_mix, openmc.Material)
        for nuc in material_to_mix.get_nuclides():
            check_value("check nuclide exists", nuc, self.operator.nuclides_with_data)
        self.material_to_mix = material_to_mix

        check_value('units', units, ('grams', 'atoms', 'cc', 'cm3'))
        self.units = units

        check_type("dilute", dilute, bool)
        self.dilute = dilute

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
            check_value("Material", str(val.id), self.burn_mats)

        elif isinstance(val, str):
            if val.isnumeric():
                check_value("Material id", val, self.burn_mats)
                val = [mat for mat in self.model.materials if mat.id == int(val)][0]

            else:
                check_value(
                    "Material name",
                    val,
                    [mat.name for mat in self.model.materials if mat.depletable],
                )
                val = [mat for mat in self.model.materials if mat.name == val][0]

        elif isinstance(val, int):
            check_value("Material id", str(val), self.burn_mats)
            val = [mat for mat in self.model.materials if mat.id == val][0]

        else:
            ValueError(f"Material: {val} is not supported")

        return val

    def search_for_keff(self, x, step_index):
        """Perform the criticality search on parametric material variable.

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
        root : float
             Search_for_keff returned root value
        """
        # Update AtomNumber with new conc vectors x. Materials are also updated
        # even though they are also re-calculated when running the search_for_kef
        self._update_materials(x)

        # Set initial material addition to 0 and let program calculate the
        # right amount
        root = self._search_for_keff(0)

        # Update concentration vector and volumes with new value for depletion
        volumes = self._adjust_volumes(root)
        x = self._update_x_and_set_volumes(x, volumes)

        return x, root

    def _set_mix_material_volume(self, param):
        """
        Set volume in cc to `self.material_to_mix` as a function of `param` after
        converion, based on `self.units` attribute.

        Parameters
        ----------
        param : float
            grams or atoms or cc of material_to_mix to convert to cc

        """
        if self.units == 'grams':
            multiplier = 1 / self.material_to_mix.get_mass_density()
        elif self.units == 'atoms':
            multiplier = (
                self.material_to_mix.average_molar_mass
                / AVOGADRO
                / self.material_to_mix.get_mass_density()
            )
        elif self.units == 'cc' or self.units == 'cm3':
            multiplier = 1

        self.material_to_mix.volume = param * multiplier


    def _model_builder(self, param):
        """Callable function which builds a model according to a passed
        parameter. This function must return an openmc.model.Model object.

        Builds the parametric model to be passed to the
        :meth:`openmc.search.search_for_keff` method.

        Parameters
        ----------
        param : float
            Model function variable, cc of material_to_mix

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

                if int(mat_id) == self.material.id:
                    self._set_mix_material_volume(param)

                    if self.dilute:
                        # increase the total volume of the material
                        # assuming ideal materials mixing
                        vol = (
                            number_i.get_mat_volume(mat_id)
                            + self.material_to_mix.volume
                        )
                    else:
                        # we assume the volume after mix won't change
                        vol = number_i.get_mat_volume(mat_id)

                    # Total atom concentration in [#atoms/cm-b]
                    #mat_index = number_i.index_mat[mat_id]
                    #tot_conc = 1.0e-24 * sum(number_i.number[mat_index]) / vol

                    for nuc in number_i.index_nuc:
                        # check only nuclides with cross sections data
                        if nuc in self.operator.nuclides_with_data:
                            atoms = number_i[mat_id, nuc]
                            if nuc in self.material_to_mix.get_nuclides():
                                # update atoms number
                                atoms += self.material_to_mix.get_nuclide_atoms()[nuc]

                            atoms_per_bcm = 1.0e-24 * atoms / vol

                            if atoms_per_bcm > 0.0:
                                nuclides.append(nuc)
                                densities.append(atoms_per_bcm)

                else:
                    # for all other materials, still check atom density limits
                    for nuc in number_i.nuclides:
                        if nuc in self.operator.nuclides_with_data:
                            # get normalized atoms density in [atoms/b-cm]
                            atoms_per_bcm = 1.0e-24 * number_i.get_atom_density(mat_id, nuc)

                            if atoms_per_bcm > 0.0:
                                nuclides.append(nuc)
                                densities.append(atoms_per_bcm)

                # assign nuclides and densities to the in-memory model
                openmc.lib.materials[int(mat_id)].set_densities(nuclides, densities)

        # always need to return a model
        return self.model

    def _adjust_volumes(self, res):
        """Uses :meth:`openmc.deplete.reactivity_control._search_for_keff`
        solution as cc to update depletable volume if `self.dilute` attribute
        is set to True.

        Parameters
        ----------
        res : float
            Solution in cc of material_to_mix, coming from
            :meth:`openmc.deplete.reactivity_control._search_for_keff`

        Returns
        -------
        volumes : dict
            Dictionary of calculated volume, where key mat id and value
            material volume, in cm3
        """
        number_i = self.operator.number
        volumes = {}

        for mat_id in self.local_mats:
            if int(mat_id) == self.material.id:
                if self.dilute:
                    volumes[mat_id] = number_i.get_mat_volume(mat_id) + res
                else:
                    volumes[mat_id] = number_i.get_mat_volume(mat_id)
        return volumes
