"""Pure depletion operator

This module implements a pure depletion operator that uses user- provided fluxes
and one-group cross sections.

"""

import copy
from collections import OrderedDict
from warnings import warn
from itertools import product

import numpy as np
import pandas as pd
from uncertainties import ufloat

import openmc
from openmc.checkvalue import check_type, check_value, check_iterable_type
from openmc.mpi import comm
from openmc.mgxs import EnergyGroups, ArbitraryXS, FissionXS
from .abc import ReactionRateHelper, OperatorResult
from .chain import REACTIONS
from .openmc_operator import OpenMCOperator
from .helpers import ChainFissionHelper, ConstantFissionYieldHelper, SourceRateHelper

_valid_rxns = list(REACTIONS)
_valid_rxns.append('fission')


class FluxDepletionOperator(OpenMCOperator):
    """Depletion operator that uses one-group
    cross sections to calculate reaction rates.

    Instances of this class can be used to perform depletion using one-group
    cross sections and constant flux or constant power. Normally, a user needn't call methods of
    this class directly. Instead, an instance of this class is passed to an
    integrator class, such as :class:`openmc.deplete.CECMIntegrator`.

    Parameters
    ----------
    materials : openmc.Materials
        Materials to deplete.
    micro_xs : pandas.DataFrame
        DataFrame with nuclides names as index and microscopic cross section
        data in the columns. Cross section units are [cm^-2].
    chain_file : str
        Path to the depletion chain XML file.
    keff : 2-tuple of float, optional
       keff eigenvalue and uncertainty from transport calculation.
       Default is None.
    prev_results : Results, optional
        Results from a previous depletion calculation.
    normalization_mode : {"constant-power", "constant-flux"}
        Indicate how reaction rates should be calculated.
        ``"constant-power"`` uses the fission Q values from the depletion chain to
        compute the flux based on the power. ``"constant-flux"`` uses the value stored in `_normalization_helper` as the flux.
    fission_q : dict, optional
        Dictionary of nuclides and their fission Q values [eV]. If not given,
        values will be pulled from the ``chain_file``. Only applicable
        if ``"normalization_mode" == "constant-power"``.
    reduce_chain : bool, optional
        If True, use :meth:`openmc.deplete.Chain.reduce` to reduce the
        depletion chain up to ``reduce_chain_level``. Default is False.
    reduce_chain_level : int, optional
        Depth of the search when reducing the depletion chain. Only used
        if ``reduce_chain`` evaluates to true. The default value of
        ``None`` implies no limit on the depth.
    fission_yield_opts : dict of str to option, optional
        Optional arguments to pass to the `FissionYieldHelper`. Will be
        passed directly on to the helper. Passing a value of None will use
        the defaults for the associated helper.


    Attributes
    ----------
    materials : openmc.Materials
        All materials present in the model
    cross_sections : pandas.DataFrame
        Object containing one-group cross-sections.
    dilute_initial : float
        Initial atom density [atoms/cm^3] to add for nuclides that
        are zero in initial condition to ensure they exist in the decay
        chain. Only done for nuclides with reaction rates.
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
                 micro_xs,
                 chain_file,
                 keff=None,
                 normalization_mode = 'constant-flux',
                 fission_q=None,
                 prev_results=None,
                 reduce_chain=False,
                 reduce_chain_level=None,
                 fission_yield_opts=None):
        # Validate micro-xs parameters
        check_type('materials', materials, openmc.Materials)
        check_type('micro_xs', micro_xs, pd.DataFrame)
        if keff is not None:
            check_type('keff', keff, tuple, float)
            keff = ufloat(*keff)

        self._keff = keff

        helper_kwargs = dict()
        helper_kwargs = {'normalization_mode': normalization_mode,
                         'fission_yield_opts': fission_yield_opts}

        super().__init__(
            materials,
            micro_xs,
            chain_file,
            prev_results,
            fission_q=fission_q,
            helper_kwargs=helper_kwargs,
            reduce_chain=reduce_chain,
            reduce_chain_level=reduce_chain_level)

    @classmethod
    def from_nuclides(cls, volume, nuclides, micro_xs,
                      chain_file,
                      keff=None,
                      normalization_mode='constant-flux',
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
            ,Dictionary with nuclide names as keys and nuclide concentrations as
            values. Nuclide concentration units are [atom/cm^3].
        micro_xs : pandas.DataFrame
            DataFrame with nuclides names as index and microscopic cross section
            data in the columns. Cross section units are [cm^-2].
        chain_file : str
            Path to the depletion chain XML file.
        keff : 2-tuple of float, optional
           keff eigenvalue and uncertainty from transport calculation.
           Default is None.
        normalization_mode : {"constant-power", "constant-flux"}
            Indicate how reaction rates should be calculated.
            ``"constant-power"`` uses the fission Q values from the depletion
            chain to compute the flux based on the power. ``"constant-flux"``
            uses the value stored in `_normalization_helper` as the flux.
        fission_q : dict, optional
            Dictionary of nuclides and their fission Q values [eV]. If not
            given, values will be pulled from the ``chain_file``. Only
            applicable if ``"normalization_mode" == "constant-power"``.
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
            Optional arguments to pass to the `FissionYieldHelper`. Will be
            passed directly on to the helper. Passing a value of None will use
            the defaults for the associated helper.

        """
        check_type('nuclides', nuclides, dict, str)
        materials = cls._consolidate_nuclides_to_material(nuclides, volume)
        return cls(materials,
                   micro_xs,
                   chain_file,
                   keff=keff,
                   normalization_mode=normalization_mode,
                   fission_q=fission_q,
                   prev_results=prev_results,
                   reduce_chain=reduce_chain,
                   reduce_chain_level=reduce_chain_level,
                   fission_yield_opts=fission_yield_opts)


    @staticmethod
    def _consolidate_nuclides_to_material(nuclides, volume):
        """Puts nuclide list into an openmc.Materials object.

        """
        openmc.reset_auto_ids()
        mat = openmc.Material()
        for nuc, conc in nuclides.items():
            mat.add_nuclide(nuc, conc * 1e-24)  # convert to at/b-cm

        mat.volume = volume
        mat.depletable = True

        return openmc.Materials([mat])

    def _load_previous_results(self):
        """Load results from a previous depletion simulation"""
        # Reload volumes into geometry
        self.prev_res[-1].transfer_volumes(self.materials)

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


    def _get_nuclides_with_data(self, cross_sections):
        """Finds nuclides with cross section data"""
        return set(cross_sections.index)

    class _FluxDepletionNormalizationHelper(ChainFissionHelper):
        """Class for calculating one-group flux
        based on a power.

        flux = Power / X, where X = volume * sum_i(Q_i * fission_micro_xs_i * density_i)

        Parameters
        ----------
        op : openmc.deplete.FluxDepletionOperator
            Reference to the object encapsulate _FluxDepletionNormalizationHelper.
            We pass this so we don't have to duplicate the ``number`` object.

        """

        def __init__(self, op):
            rates = op.reaction_rates
            self.nuc_ind_map = {ind: nuc for nuc, ind in rates.index_nuc.items()}
            self._op = op
            super().__init__()

        def update(self, fission_rates, mat_index=None):
            """Update 'energy' produced with fission rates in a material. What
            this actually calculates is the quantity X.

            Parameters
            ----------
            fission_rates : numpy.ndarray
                fission reaction rate for each isotope in the specified
                material. Should be ordered corresponding to initial
                ``rate_index`` used in :meth:`prepare`
            mat_index : int
                Material index

        """
            volume = self._op.number.get_mat_volume(mat_index)
            densities = np.empty(np.shape(fission_rates))
            for i_nuc in self.nuc_ind_map:
                nuc = self.nuc_ind_map[i_nuc]
                densities[i_nuc] = self._op.number.get_atom_density(mat_index, nuc)
            fission_rates = fission_rates * volume * densities

            super().update(fission_rates)

    class _FluxDepletionRateHelper(ReactionRateHelper):
        """Class for generating one-group reaction rates with flux and
        one-group cross sections.

        This class does not generate tallies, and instead stores cross sections
        for each nuclide and transmutation reaction relevant for a depletion
        calculation. The reaction rate is calculated by multiplying the flux by the
        cross sections.

        Parameters
        ----------
        n_nucs : int
            Number of burnable nuclides tracked by :class:`openmc.deplete.Operator`
        n_react : int
            Number of reactions tracked by :class:`openmc.deplete.Operator`
        op : openmc.deplete.FluxDepletionOperator
            Reference to the object encapsulate _FluxDepletionRateHelper.
            We pass this so we don't have to duplicate the ``number`` object.


        Attributes
        ----------
        nuc_ind_map : dict of int to str
            Dictionary mapping the nuclide index to nuclide name
        rxn_ind_map : dict of int to str
            Dictionary mapping reaction index to reaction name

        """

        def __init__(self, n_nuc, n_react, op):
            super().__init__(n_nuc, n_react)
            rates = op.reaction_rates

            self.nuc_ind_map = {ind: nuc for nuc, ind in rates.index_nuc.items()}
            self.rxn_ind_map = {ind: rxn for rxn, ind in rates.index_rx.items()}
            self._op = op

        def generate_tallies(self, materials, scores):
            """Unused in this case"""
            pass

        def get_material_rates(self, mat_id, nuc_index, react_index):
            """Return 2D array of [nuclide, reaction] reaction rates

            Parameters
            ----------
            mat_id : int
                Unique ID for the requested material
            nuc_index : list of str
                Ordering of desired nuclides
            react_index : list of str
                Ordering of reactions
            """
            self._results_cache.fill(0.0)

            volume = self._op.number.get_mat_volume(mat_id)
            for i_nuc, i_react in product(nuc_index, react_index):
                nuc = self.nuc_ind_map[i_nuc]
                rxn = self.rxn_ind_map[i_react]
                density = self._op.number.get_atom_density(mat_id, nuc)
                self._results_cache[i_nuc,
                                    i_react] = self._op.cross_sections[rxn][nuc] * density * volume

            return self._results_cache

    def _get_helper_classes(self, helper_kwargs):
        """Get helper classes for calculating reation rates and fission yields

        Parameters
        ----------
        helper_kwargs : dict
            Keyword arguments for helper classes

        """

        normalization_mode = helper_kwargs['normalization_mode']
        fission_yield_opts = helper_kwargs.get('fission_yield_opts', {})

        self._rate_helper = self._FluxDepletionRateHelper(
            self.reaction_rates.n_nuc, self.reaction_rates.n_react, self)
        if normalization_mode == "constant-power":
            self._normalization_helper = self._FluxDepletionNormalizationHelper(self)
        else:
            self._normalization_helper = SourceRateHelper()

        # Select and create fission yield helper
        fission_helper = ConstantFissionYieldHelper
        if fission_yield_opts is None:
            fission_yield_opts = {}
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
            Power in [W] or source rate in [neutron/sec]

        Returns
        -------
        openmc.deplete.OperatorResult
            Eigenvalue and reaction rates resulting from transport operator

        """

        self._update_materials_and_nuclides(vec)

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

                                      ' at/barn-cm)')
                            number_i[mat, nuc] = 0.0

    @staticmethod
    def create_micro_xs_from_data_array(nuclides, reactions, data, units='barn'):
        """
        Creates a ``micro_xs`` parameter from a dictionary.

        Parameters
        ----------
        nuclides : list of str
            List of nuclide symbols for that have data for at least one
            reaction.
        reactions : list of str
            List of reactions. All reactions must match those in ``chain.REACTONS``
        data : ndarray of floats
            Array containing one-group microscopic cross section information for each
            nuclide and reaction.
        units : {'barn', 'cm^2'}, optional
            Units for microscopic cross section data. Defaults to ``barn``.

        Returns
        -------
        micro_xs : pandas.DataFrame
            A DataFrame object correctly formatted for use in ``FluxOperator``
        """

        # Validate inputs
        if data.shape != (len(nuclides), len(reactions)):
            raise ValueError(
                f'Nuclides list of length {len(nuclides)} and '
                f'reactions array of length {len(reactions)} do not '
                f'match dimensions of data array of shape {data.shape}')

        FluxDepletionOperator._validate_micro_xs_inputs(
            nuclides, reactions, data)

        # Convert to cm^2
        if units == 'barn':
            data *= 1e-24

        return pd.DataFrame(index=nuclides, columns=reactions, data=data)

    @staticmethod
    def create_micro_xs_from_csv(csv_file, units='barn'):
        """
        Create the ``micro_xs`` parameter from a ``.csv`` file.

        Parameters
        ----------
        csv_file : str
            Relative path to csv-file containing microscopic cross section
            data.
        units : {'barn', 'cm^2'}, optional
            Units for microscopic cross section data. Defaults to ``barn``.

        Returns
        -------
        micro_xs : pandas.DataFrame
            A DataFrame object correctly formatted for use in ``FluxOperator``

        """
        micro_xs = pd.read_csv(csv_file, index_col=0)

        FluxDepletionOperator._validate_micro_xs_inputs(list(micro_xs.index),
                                                        list(micro_xs.columns),
                                                        micro_xs.to_numpy())

        if units == 'barn':
            micro_xs *= 1e-24

        return micro_xs

    # Convenience function for the micro_xs static methods
    @staticmethod
    def _validate_micro_xs_inputs(nuclides, reactions, data):
        check_iterable_type('nuclides', nuclides, str)
        check_iterable_type('reactions', reactions, str)
        check_type('data', data, np.ndarray, expected_iter_type=float)
        for reaction in reactions:
            check_value('reactions', reaction, _valid_rxns)


    @staticmethod
    def generate_1g_cross_sections(model, reaction_domain, reactions=['(n,gamma)', '(n,2n)', '(n,p)', '(n,a)', '(n,3n)', '(n,4n)', 'fission'],  energy_bounds=(0, 20e6), write_to_csv=True, filename='micro_xs.csv'):
        """Helper function to generate a one-group cross-section dataframe
        using OpenMC. Note that the ``openmc`` C executable must be compiled.

        Parameters
        ----------
        model : openmc.model.Model
            OpenMC model object. Must contain geometry, materials, and settings.
        reaction_domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.RegularMesh
            Domain in which to tally reaction rates.
        reactions : list of str, optional
            Reaction names to tally
        energy_bound : 2-tuple of float, optional
            Bounds for the energy group.
        write_to_csv : bool, optional
            Option to write the DataFrame to a `.csv` file. If `False`,
            returns the dataframe object.
        filename : str
            Name for csv file. Only applicable if ``write_to_csv == True``

        Returns
        -------
        None or pandas.DataFrame

        """
        groups = EnergyGroups()
        groups.group_edges = np.array(list(energy_bounds))

        # Set up the reaction tallies
        tallies = openmc.Tallies()
        xs = dict()
        for rxn in reactions:
            if rxn == 'fission':
                xs[rxn] = FissionXS(domain=reaction_domain, groups=groups, by_nuclide=True)
            else:
                xs[rxn] = ArbitraryXS(rxn, domain=reaction_domain, groups=groups, by_nuclide=True)
            tallies += xs[rxn].tallies.values()

        model.tallies = tallies
        statepoint_path = model.run()

        sp = openmc.StatePoint(statepoint_path)

        for rxn in xs:
            xs[rxn].load_from_statepoint(sp)

        sp.close()

        # Build the DataFrame
        micro_xs = pd.DataFrame()
        for rxn in xs:
            df = xs[rxn].get_pandas_dataframe(xs_type='micro')
            df.index = df['nuclide']
            df.drop(['nuclide', xs[rxn].domain_type, 'group in', 'std. dev.'], axis=1, inplace=True)
            df.rename({'mean':rxn}, axis=1, inplace=True)
            micro_xs = pd.concat([micro_xs, df], axis=1)

        if write_to_csv:
            micro_xs.to_csv(filename)
            return None
        else:
            return micro_xs
