"""Transport-independent transport operator for depletion.

This module implements a transport operator that runs independently of any
transport solver by using user-provided multigroup fluxes and cross sections.

"""

from __future__ import annotations
from collections.abc import Iterable
import copy
from warnings import warn

import numpy as np
from uncertainties import ufloat

import openmc
from openmc.checkvalue import check_type
from openmc.mpi import comm
from .abc import ReactionRateHelper, OperatorResult
from .chain import REACTIONS
from .openmc_operator import OpenMCOperator
from .pool import _distribute
from .microxs import MicroXS
from .results import Results
from .helpers import ChainFissionHelper, ConstantFissionYieldHelper, SourceRateHelper


class IndependentOperator(OpenMCOperator):
    """Transport-independent transport operator based on multigroup data.

    Instances of this class can be used to perform depletion using multigroup
    cross sections and multigroup fluxes. Normally, a user needn't call methods
    of this class directly. Instead, an instance of this class is passed to an
    integrator class, such as :class:`openmc.deplete.CECMIntegrator`.

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
    chain_file : PathLike or Chain, optional
        Path to the depletion chain XML file or instance of openmc.deplete.Chain.
        Defaults to ``openmc.config['chain_file']``.
    keff : 2-tuple of float, optional
       keff eigenvalue and uncertainty from transport calculation. When not
       provided and every :class:`~openmc.deplete.MicroXS` instance contains
       both 'fission' and 'nu-fission' cross sections as well as all
       transmutation reactions defined by the depletion chain, the infinite
       multiplication factor of the depletable materials is estimated
       automatically from the material compositions and multigroup cross
       sections at each depletion step.

       .. versionchanged:: 0.16.1
           k-infinity is now estimated automatically when ``keff`` is not
           given and the required cross sections are present.
    prev_results : Results, optional
        Results from a previous depletion calculation.
    normalization_mode : {"fission-q", "source-rate"}
        Indicate how reaction rates should be calculated. ``"fission-q"`` uses
        the fission Q values from the depletion chain to compute the flux based
        on the power. ``"source-rate"`` uses a the source rate (assumed to be
        neutron flux) to calculate the reaction rates.
    fission_q : dict, optional
        Dictionary of nuclides and their fission Q values [eV]. If not given,
        values will be pulled from the ``chain_file``. Only applicable if
        ``"normalization_mode" == "fission-q"``.
    reduce_chain_level : int, optional
        Depth of the search when reducing the depletion chain. The default
        value of ``None`` implies no limit on the depth.
    fission_yield_opts : dict of str to option, optional
        Optional arguments to pass to the
        :class:`openmc.deplete.helpers.FissionYieldHelper` object. Will be
        passed directly on to the helper. Passing a value of None will use the
        defaults for the associated helper.

    Attributes
    ----------
    materials : openmc.Materials
        All materials present in the model
    cross_sections : list of MicroXS
        Object containing multigroup cross-sections in [b] for each material.
    output_dir : pathlib.Path
        Path to output directory to save results.
    round_number : bool
        Whether or not to round output to OpenMC to 8 digits. Useful in testing,
        as OpenMC is incredibly sensitive to exact values.
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
        Results from a previous depletion calculation. ``None`` if no results
        are to be used.

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

        # Auto-detect k-infinity capability: estimate kinf when keff is not
        # provided and every MicroXS contains fission + nu-fission data.
        self._calculate_kinf = (
            keff is None
            and len(micros) > 0
            and all('fission' in m.reactions and 'nu-fission' in m.reactions
                    for m in micros)
        )

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
            reduce_chain_level=reduce_chain_level)

        # The k-infinity estimate divides the neutron production rate by the
        # neutron loss rate, so the loss term must include every absorption
        # channel that the depletion chain will use. If a MicroXS is missing
        # some of the chain's transmutation reactions (e.g., only 'fission'
        # and 'nu-fission' were tallied), the ratio would silently degenerate
        # toward nu-bar rather than k-infinity. Note that self.chain only
        # exists after the super().__init__() call above (which also applies
        # any chain reduction), so this check must come here.
        if self._calculate_kinf:
            chain_rxns = set(self.chain.reactions)
            for m in micros:
                if not chain_rxns <= set(m.reactions):
                    missing = chain_rxns - set(m.reactions)
                    warn(f'Disabling k-infinity estimate: MicroXS is missing '
                         f'chain reactions {missing}. The estimate requires '
                         f'all absorption channels to be present.')
                    self._calculate_kinf = False
                    break

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
        chain_file : PathLike or Chain, optional
            Path to the depletion chain XML file or instance of
            openmc.deplete.Chain. Defaults to ``openmc.config['chain_file']``.
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
        reduce_chain_level : int, optional
            Depth of the search when reducing the depletion chain. The default
            value of ``None`` implies no limit on the depth.
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
        """Finds nuclides with cross section data

        Parameters
        ----------
        cross_sections : iterable of :class`~openmc.deplete.MicroXS`
            List of multigroup cross-section data.

        Returns
        -------
        nuclides : set of str
            Set of nuclide names that have cross section data

        """
        return set(cross_sections[0].nuclides)

    class _IndependentRateHelper(ReactionRateHelper):
        """Class for generating reaction rates with multigroup fluxes and
        multigroup cross sections.

        This class does not generate tallies and instead stores cross sections
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
                if nuc not in xs._index_nuc:
                    continue
                for i_rx in react_index:
                    rx = self.rx_ind_map[i_rx]
                    if rx not in xs._index_rx:
                        continue

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

    def __call__(self, vec, source_rate) -> OperatorResult:
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
            if self._calculate_kinf:
                keff = self._estimate_k_inf()
            else:
                keff = ufloat(0.0, 0.0)
            return OperatorResult(keff, rates)

        rates = self._calculate_reaction_rates(source_rate)
        if self._calculate_kinf:
            keff = self._estimate_k_inf()
        else:
            keff = self._keff

        op_result = OperatorResult(keff, rates)
        return copy.deepcopy(op_result)

    def _estimate_k_inf(self):
        r"""Estimate the infinite multiplication factor of the depletable
        materials.

        The estimate is computed as the ratio of the neutron production rate
        to the neutron loss rate summed over the *depletable materials only*:

        .. math::
            k_\infty = \frac{\displaystyle\sum_m \frac{1}{V_m} \sum_i N_{m,i}
                             \sum_g (\nu\sigma_f)_{m,i,g}\, \phi_{m,g}}
                            {\displaystyle\sum_m \frac{1}{V_m} \sum_i N_{m,i}
                             \sum_j (1 - x_j) \sum_g \sigma_{m,i,j,g}\,
                             \phi_{m,g}}

        where the index :math:`m` runs over the depletable materials,
        :math:`i` over the nuclides with cross-section data, :math:`j` over
        the transmutation reactions, and :math:`g` over the energy groups.
        :math:`N_{m,i}` is the number of atoms of nuclide :math:`i` in
        material :math:`m`, :math:`V_m` is the material volume,
        :math:`\phi_{m,g}` is the volume-integrated multigroup flux from the
        transport run, :math:`(\nu\sigma_f)_{m,i,g}` is the fission neutron
        production cross section, :math:`\sigma_{m,i,j,g}` is the cross
        section of transmutation reaction :math:`j`, and :math:`x_j` is the
        number of neutrons emitted by reaction :math:`j`.

        **This is not k-eff.** The balance above contains no leakage term, so
        it relates to the effective multiplication factor as
        :math:`k_\infty = k_\text{eff} / (1 - L)` where :math:`L` is the
        leakage fraction. Moreover, only depletable materials contribute to
        the loss term: for models that also contain non-depletable materials
        (moderator, cladding, reflector, ...), absorption in those materials
        is not accounted for, and the estimate will be *higher* than the true
        k-infinity of the full system. In other words, the estimate assumes
        that all relevant absorption happens in the depletable materials.

        The treatment of (n,xn) reactions follows from writing the
        multiplication factor as

        .. math::
            k_\text{eff} = \frac{P}{A + L - X}

        where :math:`P` is the fission neutron production rate, :math:`A` the
        absorption rate, :math:`L` the leakage rate, and :math:`X` the net
        neutron production rate from (n,xn) reactions. The denominator uses
        "reduced absorption" :math:`A - X`, which is exactly the convention
        used by OpenMC's k-eff estimators: neutrons produced in (n,xn)
        reactions are not counted as production; instead each (n,xn) reaction
        with :math:`x` neutrons out contributes :math:`(1 - x)` times its
        rate to the loss term, giving :math:`A - X` in a single pass over the
        reactions.

        Assumptions: the multigroup fluxes are those obtained from the
        transport run (and are not recomputed as the compositions change),
        non-depletable materials do not deplete, and absorption outside the
        depletable materials is ignored entirely. The bias introduced by that
        last assumption stays roughly constant over the depletion only to the
        extent that the background absorption itself does.

        Returns
        -------
        uncertainties.UFloat
            Estimated k-infinity with zero uncertainty

        """
        production = 0.0
        loss = 0.0
        for mat in self.local_mats:
            i_mat = self._mat_index_map[mat]
            flux = self.fluxes[i_mat]
            micro_xs = self.cross_sections[i_mat]

            # Convert total atoms and volume-integrated flux to rates per
            # unit volume, consistent with _calculate_reaction_rates
            volume_b_cm = 1e24 * self.number.get_mat_volume(mat)

            for nuc in micro_xs.nuclides:
                if nuc not in self.number.index_nuc:
                    continue
                atoms = self.number[mat, nuc]
                if atoms <= 0.0:
                    continue
                for rxn in micro_xs.reactions:
                    rate = (atoms * (micro_xs[nuc, rxn] * flux).sum()
                            / volume_b_cm)
                    if rxn == 'nu-fission':
                        production += rate
                    elif rxn == 'damage-energy':
                        pass
                    elif rxn in REACTIONS:
                        n_out = REACTIONS[rxn].neutrons_out
                        loss += (1 - n_out) * rate
                    else:
                        loss += rate

        # Sum contributions over all MPI processes
        production = comm.allreduce(production)
        loss = comm.allreduce(loss)

        if loss <= 0.0:
            warn('Unable to estimate k-infinity because the total neutron '
                 'loss rate is zero. Check that the supplied MicroXS data '
                 'contains absorption reactions for the nuclides present.')
            return ufloat(0.0, 0.0)
        return ufloat(production / loss, 0.0)

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
