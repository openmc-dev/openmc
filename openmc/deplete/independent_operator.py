"""Transport-independent transport operator for depletion.

This module implements a transport operator that runs independently of any
transport solver by using user-provided one-group cross sections.

"""

from __future__ import annotations
from collections.abc import Iterable
import copy

import numpy as np
from uncertainties import ufloat

import openmc
from openmc.checkvalue import check_type
from openmc.mpi import comm
from .abc import ReactionRateHelper, OperatorResult
from .openmc_operator import OpenMCOperator
from .pool import _distribute
from .microxs import MicroXS
from .results import Results
from .helpers import ChainFissionHelper, ConstantFissionYieldHelper, SourceRateHelper


class IndependentOperator(OpenMCOperator):
    """Transport-independent transport operator that uses one-group cross
    sections to calculate reaction rates.

    Instances of this class can be used to perform depletion using one-group
    cross sections and constant flux or constant power. Normally, a user needn't
    call methods of this class directly. Instead, an instance of this class is
    passed to an integrator class, such as
    :class:`openmc.deplete.CECMIntegrator`.

    Note that passing an empty :class:`~openmc.deplete.MicroXS` instance to the
    ``micro_xs`` argument allows a decay-only calculation to be run.

    .. versionadded:: 0.13.1

    .. versionchanged:: 0.14.0
        Arguments updated to include list of fluxes and microscopic cross
        sections.

    Parameters
    ----------
    materials : iterable of openmc.Material
        Materials to deplete.
    fluxes : list of numpy.ndarray
        Flux in each group in [n-cm/src] for each domain
    micros : list of MicroXS
        Cross sections in [b] for each domain. If the
        :class:`~openmc.deplete.MicroXS` object is empty, a decay-only
        calculation will be run.
    chain_file : str
        Path to the depletion chain XML file. Defaults to
        ``openmc.config['chain_file']``.
    keff : 2-tuple of float, optional
       keff eigenvalue and uncertainty from transport calculation.
       Default is None.
    prev_results : Results, optional
        Results from a previous depletion calculation.
    normalization_mode : {"fission-q", "source-rate"}
        Indicate how reaction rates should be calculated.
        ``"fission-q"`` uses the fission Q values from the depletion chain to
        compute the flux based on the power. ``"source-rate"`` uses a the
        source rate (assumed to be neutron flux) to calculate the
        reaction rates.
    fission_q : dict, optional
        Dictionary of nuclides and their fission Q values [eV]. If not given,
        values will be pulled from the ``chain_file``. Only applicable
        if ``"normalization_mode" == "fission-q"``.
    reduce_chain : bool, optional
        If True, use :meth:`openmc.deplete.Chain.reduce` to reduce the
        depletion chain up to ``reduce_chain_level``.
    reduce_chain_level : int, optional
        Depth of the search when reducing the depletion chain. Only used
        if ``reduce_chain`` evaluates to true. The default value of
        ``None`` implies no limit on the depth.
    fission_yield_opts : dict of str to option, optional
        Optional arguments to pass to the
        :class:`openmc.deplete.helpers.FissionYieldHelper` object. Will be
        passed directly on to the helper. Passing a value of None will use
        the defaults for the associated helper.

    Attributes
    ----------
    materials : openmc.Materials
        All materials present in the model
    cross_sections : MicroXS
        Object containing one-group cross-sections in [cm^2].
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

    """

    def __init__(self,
                 materials,
                 fluxes,
                 micros,
                 chain_file=None,
                 keff=None,
                 normalization_mode='fission-q',
                 fission_q=None,
                 prev_results=None,
                 reduce_chain=False,
                 reduce_chain_level=None,
                 fission_yield_opts=None):
        # Validate micro-xs parameters
        check_type('materials', materials, Iterable, openmc.Material)
        check_type('micros', micros, Iterable, MicroXS)
        materials = openmc.Materials(materials)

        if not (len(fluxes) == len(micros) == len(materials)):
            msg = (f'The length of fluxes ({len(fluxes)}) should be equal to '
                   f'the length of micros ({len(micros)}) and the length of '
                   f'materials ({len(materials)}).')
            raise ValueError(msg)

        if keff is not None:
            check_type('keff', keff, tuple, float)
            keff = ufloat(*keff)

        self._keff = keff

        if fission_yield_opts is None:
            fission_yield_opts = {}
        helper_kwargs = {'normalization_mode': normalization_mode,
                         'fission_yield_opts': fission_yield_opts}

        # Sort fluxes and micros in same order that materials get sorted
        index_sort = np.argsort([mat.id for mat in materials])
        fluxes = [fluxes[i] for i in index_sort]
        micros = [micros[i] for i in index_sort]

        self.fluxes = fluxes
        super().__init__(
            materials=materials,
            cross_sections=micros,
            chain_file=chain_file,
            prev_results=prev_results,
            fission_q=fission_q,
            helper_kwargs=helper_kwargs,
            reduce_chain=reduce_chain,
            reduce_chain_level=reduce_chain_level)

    @classmethod
    def from_nuclides(cls, volume, nuclides,
                      flux,
                      micro_xs,
                      chain_file=None,
                      nuc_units='atom/b-cm',
                      keff=None,
                      normalization_mode='fission-q',
                      fission_q=None,
                      prev_results=None,
                      reduce_chain=False,
                      reduce_chain_level=None,
                      fission_yield_opts=None):
        """
        Alternate constructor from a dictionary of nuclide concentrations

        volume : float
            Volume of the material being depleted in [cm^3]
        nuclides : dict of str to float
            Dictionary with nuclide names as keys and nuclide concentrations as
            values.
        flux : numpy.ndarray
            Flux in each group in [n-cm/src]
        micro_xs : MicroXS
            Cross sections in [b]. If the :class:`~openmc.deplete.MicroXS`
            object is empty, a decay-only calculation will be run.
        chain_file : str, optional
            Path to the depletion chain XML file. Defaults to
            ``openmc.config['chain_file']``.
        nuc_units : {'atom/cm3', 'atom/b-cm'}, optional
            Units for nuclide concentration.
        keff : 2-tuple of float, optional
           keff eigenvalue and uncertainty from transport calculation.
           Default is None.
        normalization_mode : {"fission-q", "source-rate"}
            Indicate how reaction rates should be calculated.
            ``"fission-q"`` uses the fission Q values from the depletion
            chain to compute the flux based on the power. ``"source-rate"`` uses
            the source rate (assumed to be neutron flux) to calculate the
            reaction rates.
        fission_q : dict, optional
            Dictionary of nuclides and their fission Q values [eV]. If not
            given, values will be pulled from the ``chain_file``. Only
            applicable if ``"normalization_mode" == "fission-q"``.
        prev_results : Results, optional
            Results from a previous depletion calculation.
        reduce_chain : bool, optional
            If True, use :meth:`openmc.deplete.Chain.reduce` to reduce the
            depletion chain up to ``reduce_chain_level``. Default is False.
        reduce_chain_level : int, optional
            Depth of the search when reducing the depletion chain. Only used
            if ``reduce_chain`` evaluates to true. The default value of
            ``None`` implies no limit on the depth.
        fission_yield_opts : dict of str to option, optional
            Optional arguments to pass to the
            :class:`openmc.deplete.helpers.FissionYieldHelper` class. Will be
            passed directly on to the helper. Passing a value of None will use
            the defaults for the associated helper.

        """
        check_type('nuclides', nuclides, dict, str)
        materials = cls._consolidate_nuclides_to_material(nuclides, nuc_units, volume)
        fluxes = [flux]
        micros = [micro_xs]
        return cls(materials,
                   fluxes,
                   micros,
                   chain_file,
                   keff=keff,
                   normalization_mode=normalization_mode,
                   fission_q=fission_q,
                   prev_results=prev_results,
                   reduce_chain=reduce_chain,
                   reduce_chain_level=reduce_chain_level,
                   fission_yield_opts=fission_yield_opts)

    @staticmethod
    def _consolidate_nuclides_to_material(nuclides, nuc_units, volume):
        """Puts nuclide list into an openmc.Materials object.

        """
        mat = openmc.Material()
        if nuc_units == 'atom/b-cm':
            for nuc, conc in nuclides.items():
                mat.add_nuclide(nuc, conc)
        elif nuc_units == 'atom/cm3':
            for nuc, conc in nuclides.items():
                mat.add_nuclide(nuc, conc * 1e-24)  # convert to at/b-cm
        else:
            raise ValueError(f"Unit '{nuc_units}' is invalid.")

        mat.volume = volume
        mat.depletable = True

        return openmc.Materials([mat])

    def _load_previous_results(self):
        """Load results from a previous depletion simulation"""
        # Reload volumes into geometry
        model = openmc.Model(materials=self.materials)
        self.prev_res[-1].transfer_volumes(model)
        self.materials = model.materials

        # Store previous results in operator
        # Distribute reaction rates according to those tracked
        # on this process
        if comm.size != 1:
            prev_results = self.prev_res
            self.prev_res = Results()
            mat_indexes = _distribute(range(len(self.burnable_mats)))
            for res_obj in prev_results:
                new_res = res_obj.distribute(self.local_mats, mat_indexes)
                self.prev_res.append(new_res)

    def _get_nuclides_with_data(self, cross_sections: list[MicroXS]) -> set[str]:
        """Finds nuclides with cross section data"""
        return set(cross_sections[0].nuclides)

    class _IndependentRateHelper(ReactionRateHelper):
        """Class for generating one-group reaction rates with flux and
        one-group cross sections.

        This class does not generate tallies, and instead stores cross sections
        for each nuclide and transmutation reaction relevant for a depletion
        calculation. The reaction rate is calculated by multiplying the flux by
        the cross sections.

        Parameters
        ----------
        op : openmc.deplete.IndependentOperator
            Reference to the object encapsulate _IndependentRateHelper.
            We pass this so we don't have to duplicate the
            :attr:`IndependentOperator.number` object.

        Attributes
        ----------
        nuc_ind_map : dict of int to str
            Dictionary mapping the nuclide index to nuclide name
        rx_ind_map : dict of int to str
            Dictionary mapping reaction index to reaction name

        """

        def __init__(self, op: IndependentOperator):
            rates = op.reaction_rates
            super().__init__(rates.n_nuc, rates.n_react)

            self.nuc_ind_map = {ind: nuc for nuc, ind in rates.index_nuc.items()}
            self.rx_ind_map = {ind: rxn for rxn, ind in rates.index_rx.items()}
            self._op = op

        def generate_tallies(self, materials, scores):
            """Unused in this case"""
            pass

        def reset_tally_means(self):
            """Unused in this case"""
            pass

        def get_material_rates(self, mat_index, nuc_index, react_index):
            """Return 2D array of [nuclide, reaction] reaction rates

            Parameters
            ----------
            mat_index : int
                Index for the material
            nuc_index : list of str
                Ordering of desired nuclides
            react_index : list of str
                Ordering of reactions
            """
            self._results_cache.fill(0.0)

            # Get flux and microscopic cross sections from operator
            flux = self._op.fluxes[mat_index]
            xs = self._op.cross_sections[mat_index]

            for i_nuc in nuc_index:
                nuc = self.nuc_ind_map[i_nuc]
                for i_rx in react_index:
                    rx = self.rx_ind_map[i_rx]

                    # Determine reaction rate by multiplying xs in [b] by flux
                    # in [n-cm/src] to give [(reactions/src)*b-cm/atom]
                    self._results_cache[i_nuc, i_rx] = (xs[nuc, rx] * flux).sum()

            return self._results_cache

    def _get_helper_classes(self, helper_kwargs):
        """Get helper classes for calculating reaction rates and fission yields

        Parameters
        ----------
        helper_kwargs : dict
            Keyword arguments for helper classes

        """

        normalization_mode = helper_kwargs['normalization_mode']
        fission_yield_opts = helper_kwargs['fission_yield_opts']

        self._rate_helper = self._IndependentRateHelper(self)
        if normalization_mode == "fission-q":
            self._normalization_helper = ChainFissionHelper()
        else:
            self._normalization_helper = SourceRateHelper()

        # Select and create fission yield helper
        fission_helper = ConstantFissionYieldHelper
        self._yield_helper = fission_helper.from_operator(
            self, **fission_yield_opts)

    def initial_condition(self):
        """Performs final setup and returns initial condition.

        Returns
        -------
        list of numpy.ndarray
            Total density for initial conditions.
        """

        # Return number density vector
        return super().initial_condition(self.materials)

    def __call__(self, vec, source_rate):
        """Obtain the reaction rates

        Parameters
        ----------
        vec : list of numpy.ndarray
            Total atoms to be used in function.
        source_rate : float
            Power in [W] or flux in [neutron/cm^2-s]

        Returns
        -------
        openmc.deplete.OperatorResult
            Eigenvalue and reaction rates resulting from transport operator

        """

        self._update_materials_and_nuclides(vec)

        # If the source rate is zero, return zero reaction rates
        if source_rate == 0.0:
            rates = self.reaction_rates.copy()
            rates.fill(0.0)
            return OperatorResult(ufloat(0.0, 0.0), rates)

        rates = self._calculate_reaction_rates(source_rate)
        keff = self._keff

        op_result = OperatorResult(keff, rates)
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
