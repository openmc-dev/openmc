"""Transport-coupled transport operator for depletion.

This module implements a transport operator coupled to OpenMC's transport solver
so that it can be used by depletion integrators. The implementation makes use of
the Python bindings to OpenMC's C API so that reading tally results and updating
material number densities is all done in-memory instead of through the
filesystem.

"""

import copy
from warnings import warn

import numpy as np
from uncertainties import ufloat

import openmc
from openmc.checkvalue import check_value
from openmc.data import DataLibrary
from openmc.exceptions import DataError
import openmc.lib
from openmc.mpi import comm
from .abc import OperatorResult
from .openmc_operator import OpenMCOperator
from .pool import _distribute
from .results import Results
from .helpers import (
    DirectReactionRateHelper, ChainFissionHelper, ConstantFissionYieldHelper,
    FissionYieldCutoffHelper, AveragedFissionYieldHelper, EnergyScoreHelper,
    SourceRateHelper, FluxCollapseHelper)


__all__ = ["CoupledOperator", "Operator", "OperatorResult"]


def _find_cross_sections(model: str | None = None):
    """Determine cross sections to use for depletion

    Parameters
    ----------
    model : openmc.model.Model, optional
        Reactor model

    """
    if model:
        if model.materials and model.materials.cross_sections is not None:
            # Prefer info from Model class if available
            return model.materials.cross_sections

    # otherwise fallback to environment variable
    cross_sections = openmc.config.get("cross_sections")
    if cross_sections is None:
        raise DataError(
            "Cross sections were not specified in Model.materials and "
            "openmc.config['cross_sections'] is not set."
        )
    return cross_sections


def _get_nuclides_with_data(cross_sections):
    """Loads cross_sections.xml file to find nuclides with neutron data

    Parameters
    ----------
    cross_sections : str
        Path to cross_sections.xml file

    Returns
    -------
    nuclides : set of str
        Set of nuclide names that have cross section data

    """
    nuclides = set()
    data_lib = DataLibrary.from_xml(cross_sections)
    for library in data_lib.libraries:
        if library['type'] != 'neutron':
            continue
        for name in library['materials']:
            if name not in nuclides:
                nuclides.add(name)

    return nuclides


class CoupledOperator(OpenMCOperator):
    """Transport-coupled transport operator.

    Instances of this class can be used to perform transport-coupled depletion
    using OpenMC's transport solver. Normally, a user needn't call methods of
    this class directly. Instead, an instance of this class is passed to an
    integrator class, such as :class:`openmc.deplete.CECMIntegrator`.

    .. versionchanged:: 0.13.0
        The geometry and settings parameters have been replaced with a
        model parameter that takes a :class:`~openmc.model.Model` object

    .. versionchanged:: 0.13.1
        Name changed from ``Operator`` to ``CoupledOperator``

    Parameters
    ----------
    model : openmc.model.Model
        OpenMC model object
    chain_file : str | Chain, optional
        Path to the depletion chain XML file or instance of openmc.deplete.Chain.
        Defaults to ``openmc.config['chain_file']``.
    prev_results : Results, optional
        Results from a previous depletion calculation. If this argument is
        specified, the depletion calculation will start from the latest state
        in the previous results.
    diff_burnable_mats : bool, optional
        Whether to differentiate burnable materials with multiple instances.
        Volumes are divided equally from the original material volume.
    normalization_mode : {"energy-deposition", "fission-q", "source-rate"}
        Indicate how tally results should be normalized. ``"energy-deposition"``
        computes the total energy deposited in the system and uses the ratio of
        the power to the energy produced as a normalization factor.
        ``"fission-q"`` uses the fission Q values from the depletion chain to
        compute the  total energy deposited. ``"source-rate"`` normalizes
        tallies based on the source rate (for fixed source calculations).
    fission_q : dict, optional
        Dictionary of nuclides and their fission Q values [eV]. If not given,
        values will be pulled from the ``chain_file``. Only applicable
        if ``"normalization_mode" == "fission-q"``
    fission_yield_mode : {"constant", "cutoff", "average"}
        Key indicating what fission product yield scheme to use. The
        key determines what fission energy helper is used:

        * "constant": :class:`~openmc.deplete.helpers.ConstantFissionYieldHelper`
        * "cutoff": :class:`~openmc.deplete.helpers.FissionYieldCutoffHelper`
        * "average": :class:`~openmc.deplete.helpers.AveragedFissionYieldHelper`

        The documentation on these classes describe their methodology
        and differences. Default: ``"constant"``
    fission_yield_opts : dict of str to option, optional
        Optional arguments to pass to the helper determined by
        ``fission_yield_mode``. Will be passed directly on to the
        helper. Passing a value of None will use the defaults for
        the associated helper.
    reaction_rate_mode : {"direct", "flux"}, optional
        Indicate how one-group reaction rates should be calculated. The "direct"
        method tallies transmutation reaction rates directly. The "flux" method
        tallies a multigroup flux spectrum and then collapses one-group reaction
        rates after a transport solve (with an option to tally some reaction
        rates directly).

        .. versionadded:: 0.12.1
    reaction_rate_opts : dict, optional
        Keyword arguments that are passed to the reaction rate helper class.
        When ``reaction_rate_mode`` is set to "flux", energy group boundaries
        can be set using the "energies" key. See the
        :class:`~openmc.deplete.helpers.FluxCollapseHelper` class for all
        options.

        .. versionadded:: 0.12.1
    reduce_chain_level : int, optional
        Depth of the search when reducing the depletion chain. The default
        value of ``None`` implies no limit on the depth.

        .. versionadded:: 0.12
    diff_volume_method : str
        Specifies how the volumes of the new materials should be found. Default
        is to 'divide equally' which divides the original material volume
        equally between the new materials, 'match cell' sets the volume of the
        material to volume of the cell they fill.

        .. versionadded:: 0.14.0

    Attributes
    ----------
    model : openmc.model.Model
        OpenMC model object
    output_dir : pathlib.Path
        Path to output directory to save results.
    round_number : bool
        Whether or not to round output to OpenMC to 8 digits.
        Useful in testing, as OpenMC is incredibly sensitive to exact values.
    number : openmc.deplete.AtomNumber
        Total number of atoms in simulation.
    nuclides_with_data : set of str
        A set listing all unique nuclides available from cross_sections.xml.
    chain : openmc.deplete.Chain
        The depletion chain information necessary to form matrices and tallies.
    reaction_rates : openmc.deplete.ReactionRates
        Reaction rates from the last operator step.
    burnable_mats : list of str
        All burnable material IDs
    heavy_metal : float
        Initial heavy metal inventory [g]
    local_mats : list of str
        All burnable material IDs being managed by a single process
    prev_res : Results or None
        Results from a previous depletion calculation. ``None`` if no
        results are to be used.
    cleanup_when_done : bool
        Whether to finalize and clear the shared library memory when the
        depletion operation is complete. Defaults to clearing the library.
    """
    _fission_helpers = {
        "average": AveragedFissionYieldHelper,
        "constant": ConstantFissionYieldHelper,
        "cutoff": FissionYieldCutoffHelper,
    }

    def __init__(self, model, chain_file=None, prev_results=None,
                 diff_burnable_mats=False, diff_volume_method="divide equally",
                 normalization_mode="fission-q", fission_q=None,
                 fission_yield_mode="constant", fission_yield_opts=None,
                 reaction_rate_mode="direct", reaction_rate_opts=None,
                 reduce_chain_level=None):

        # check for old call to constructor
        if isinstance(model, openmc.Geometry):
            msg = "As of version 0.13.0 openmc.deplete.CoupledOperator " \
                "requires an openmc.Model object rather than the " \
                "openmc.Geometry and openmc.Settings parameters. Please use " \
                "the geometry and settings objects passed here to create a " \
                " model with which to generate the transport Operator."
            raise TypeError(msg)

        # Determine cross sections
        cross_sections = _find_cross_sections(model)

        check_value('fission yield mode', fission_yield_mode,
                    self._fission_helpers.keys())
        check_value('normalization mode', normalization_mode,
                    ('energy-deposition', 'fission-q', 'source-rate'))
        if normalization_mode != "fission-q":
            if fission_q is not None:
                warn("Fission Q dictionary will not be used")
                fission_q = None
        self.model = model

        # determine set of materials in the model
        if not model.materials:
            model.materials = openmc.Materials(
                model.geometry.get_all_materials().values()
            )

        self.cleanup_when_done = True

        if reaction_rate_opts is None:
            reaction_rate_opts = {}
        if fission_yield_opts is None:
            fission_yield_opts = {}
        helper_kwargs = {
            'reaction_rate_mode': reaction_rate_mode,
            'normalization_mode': normalization_mode,
            'fission_yield_mode': fission_yield_mode,
            'reaction_rate_opts': reaction_rate_opts,
            'fission_yield_opts': fission_yield_opts
        }

        # Records how many times the operator has been called
        self._n_calls = 0

        super().__init__(
            materials=model.materials,
            cross_sections=cross_sections,
            chain_file=chain_file,
            prev_results=prev_results,
            diff_burnable_mats=diff_burnable_mats,
            diff_volume_method=diff_volume_method,
            fission_q=fission_q,
            helper_kwargs=helper_kwargs,
            reduce_chain_level=reduce_chain_level)

    def _differentiate_burnable_mats(self):
        """Assign distribmats for each burnable material"""

        self.model.differentiate_depletable_mats(
            diff_volume_method=self.diff_volume_method
        )

    def _load_previous_results(self):
        """Load results from a previous depletion simulation"""
        # Reload volumes into geometry
        self.prev_res[-1].transfer_volumes(self.model)

        # Store previous results in operator
        # Distribute reaction rates according to those tracked
        # on this process
        if comm.size != 1:
            prev_results = self.prev_res
            self.prev_res = Results(filename=None)
            mat_indexes = _distribute(range(len(self.burnable_mats)))
            for res_obj in prev_results:
                new_res = res_obj.distribute(self.local_mats, mat_indexes)
                self.prev_res.append(new_res)

    def _get_nuclides_with_data(self, cross_sections):
        """Loads cross_sections.xml file to find nuclides with neutron data

        Parameters
        ----------
        cross_sections : str
            Path to cross_sections.xml file

        Returns
        -------
        nuclides : set of str
            Set of nuclide names that have cross secton data

        """
        return _get_nuclides_with_data(cross_sections)

    def _get_helper_classes(self, helper_kwargs):
        """Create the ``_rate_helper``, ``_normalization_helper``, and
        ``_yield_helper`` objects.

        Parameters
        ----------
        helper_kwargs : dict
            Keyword arguments for helper classes

        """
        reaction_rate_mode = helper_kwargs['reaction_rate_mode']
        normalization_mode = helper_kwargs['normalization_mode']
        fission_yield_mode = helper_kwargs['fission_yield_mode']
        reaction_rate_opts = helper_kwargs['reaction_rate_opts']
        fission_yield_opts = helper_kwargs['fission_yield_opts']

        # Get classes to assist working with tallies
        if reaction_rate_mode == "direct":
            self._rate_helper = DirectReactionRateHelper(
                self.reaction_rates.n_nuc, self.reaction_rates.n_react)
        elif reaction_rate_mode == "flux":
            # Ensure energy group boundaries were specified
            if 'energies' not in reaction_rate_opts:
                raise ValueError(
                    "Energy group boundaries must be specified in the "
                    "reaction_rate_opts argument when reaction_rate_mode is"
                    "set to 'flux'.")

            self._rate_helper = FluxCollapseHelper(
                self.reaction_rates.n_nuc,
                self.reaction_rates.n_react,
                **reaction_rate_opts
            )
        else:
            raise ValueError("Invalid reaction rate mode.")

        if normalization_mode == "fission-q":
            self._normalization_helper = ChainFissionHelper()
        elif normalization_mode == "energy-deposition":
            score = "heating" if self.model.settings.photon_transport else "heating-local"
            self._normalization_helper = EnergyScoreHelper(score)
        else:
            self._normalization_helper = SourceRateHelper()

        # Select and create fission yield helper
        fission_helper = self._fission_helpers[fission_yield_mode]
        self._yield_helper = fission_helper.from_operator(
            self, **fission_yield_opts)

    def initial_condition(self):
        """Performs final setup and returns initial condition.

        Returns
        -------
        list of numpy.ndarray
            Total density for initial conditions.

        """

        # Create XML files
        if comm.rank == 0:
            self.model.geometry.export_to_xml()
            self.model.settings.export_to_xml()
            if self.model.plots:
                self.model.plots.export_to_xml()
            if self.model.tallies:
                self.model.tallies.export_to_xml()
            self._generate_materials_xml()

        # Initialize OpenMC library
        comm.barrier()
        if not openmc.lib.is_initialized:
            openmc.lib.init(intracomm=comm)

        # Generate tallies in memory
        materials = [openmc.lib.materials[int(i)] for i in self.burnable_mats]

        return super().initial_condition(materials)

    def _generate_materials_xml(self):
        """Creates materials.xml from self.number.

        Due to uncertainty with how MPI interacts with OpenMC API, this
        constructs the XML manually.  The long term goal is to do this
        through direct memory writing.

        """
        # Sort nuclides according to order in AtomNumber object
        nuclides = list(self.number.nuclides)
        for mat in self.materials:
            mat._nuclides.sort(key=lambda x: nuclides.index(x[0]))

        self.materials.export_to_xml(nuclides_to_ignore=self._decay_nucs)

    def __call__(self, vec, source_rate):
        """Runs a simulation.

        Simulation will abort under the following circumstances:

            1) No energy is computed using OpenMC tallies.

        Parameters
        ----------
        vec : list of numpy.ndarray
            Total atoms to be used in function.
        source_rate : float
            Power in [W] or source rate in [neutron/sec]

        Returns
        -------
        openmc.deplete.OperatorResult
            Eigenvalue and reaction rates resulting from transport operator

        """
        # Reset results in OpenMC
        openmc.lib.reset()

        # The timers are reset only if the operator has been called before.
        # This is because we call this method after loading cross sections, and
        # no transport has taken place yet. As a result, we only reset the
        # timers after the first step so as to correctly report the time spent
        # reading cross sections in the first depletion step, and from there
        # correctly report all particle tracking rates in multistep depletion
        # solvers.
        if self._n_calls > 0:
            openmc.lib.reset_timers()

        self._update_materials_and_nuclides(vec)

        # If the source rate is zero, return zero reaction rates without running
        # a transport solve
        if source_rate == 0.0:
            rates = self.reaction_rates.copy()
            rates.fill(0.0)
            return OperatorResult(ufloat(0.0, 0.0), rates)

        # Run OpenMC
        openmc.lib.run()

        # Extract results
        rates = self._calculate_reaction_rates(source_rate)

        # Get k and uncertainty
        keff = ufloat(*openmc.lib.keff())

        op_result = OperatorResult(keff, rates)

        self._n_calls += 1

        return copy.deepcopy(op_result)

    def _update_materials(self):
        """Updates material compositions in OpenMC on all processes."""

        for rank in range(comm.size):
            number_i = comm.bcast(self.number, root=rank)

            for mat in number_i.materials:
                nuclides = []
                densities = []
                for nuc in number_i.nuclides:
                    if nuc in self.nuclides_with_data:
                        val = 1.0e-24 * number_i.get_atom_density(mat, nuc)

                        # If nuclide is zero, do not add to the problem.
                        if val > 0.0:
                            if self.round_number:
                                val_magnitude = np.floor(np.log10(val))
                                val_scaled = val / 10**val_magnitude
                                val_round = round(val_scaled, 8)

                                val = val_round * 10**val_magnitude

                            nuclides.append(nuc)
                            densities.append(val)
                        else:
                            # Only output warnings if values are significantly
                            # negative. CRAM does not guarantee positive
                            # values.
                            if val < -1.0e-21:
                                print(f'WARNING: nuclide {nuc} in material'
                                      f'{mat} is negative (density = {val}'

                                      ' atom/b-cm)')

                                number_i[mat, nuc] = 0.0

                # Update densities on C API side
                mat_internal = openmc.lib.materials[int(mat)]
                mat_internal.set_densities(nuclides, densities)

                # TODO Update densities on the Python side, otherwise the
                # summary.h5 file contains densities at the first time step

    @staticmethod
    def write_bos_data(step):
        """Write a state-point file with beginning of step data

        Parameters
        ----------
        step : int
            Current depletion step including restarts

        """
        openmc.lib.statepoint_write(
            f"openmc_simulation_n{step}.h5",
            write_source=False)

    def finalize(self):
        """Finalize a depletion simulation and release resources."""
        if self.cleanup_when_done:
            openmc.lib.finalize()

    # The next few class variables and methods should be removed after one
    # release cycle or so. For now, we will provide compatibility to
    # accessing CoupledOperator.settings and CoupledOperator.geometry. In
    # the future these should stay on the Model class.

    var_warning_msg = "The CoupledOperator.{0} variable should be \
accessed through CoupledOperator.model.{0}."
    geometry_warning_msg = var_warning_msg.format("geometry")
    settings_warning_msg = var_warning_msg.format("settings")

    @property
    def settings(self):
        warn(self.settings_warning_msg, FutureWarning)
        return self.model.settings

    @settings.setter
    def settings(self, new_settings):
        warn(self.settings_warning_msg, FutureWarning)
        self.model.settings = new_settings

    @property
    def geometry(self):
        warn(self.geometry_warning_msg, FutureWarning)
        return self.model.geometry

    @geometry.setter
    def geometry(self, new_geometry):
        warn(self.geometry_warning_msg, FutureWarning)
        self.model.geometry = new_geometry


# Retain deprecated name for the time being
def Operator(*args, **kwargs):
    # warn of name change
    warn(
        "The Operator(...) class has been renamed and will "
        "be removed in a future version of OpenMC. Use "
        "CoupledOperator(...) instead.",
        FutureWarning
    )
    return CoupledOperator(*args, **kwargs)
