"""Pure depletion operator

This module implements a pure depletion operator that uses user provided fluxes
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
from .abc import TransportOperator, ReactionRateHelper, OperatorResult
from .atom_number import AtomNumber
from .chain import REACTIONS
from .openmc_operator import OpenMCOperator
from .reaction_rates import ReactionRates
from .helpers import ConstantFissionYieldHelper, SourceRateHelper

valid_rxns = list(REACTIONS)
valid_rxns.append('fission')

class FluxDepletionOperator(OpenMCOperator):
    """Depletion operator that uses a user-provided flux spectrum and one-group
    cross sections to calculate reaction rates.

    Instances of this class can be used to perform depletion using one group
    cross sections and constant flux. Normally, a user needn't call methods of
    this class directly. Instead, an instance of this class is passed to an
    integrator class, such as :class:`openmc.deplete.CECMIntegrator`

    Parameters
    ----------
    materials : openmc.Materials
        Materials to deplete.
    micro_xs : pandas.DataFrame
        DataFrame with nuclides names as index and microscopic cross section
        data in the columns. Cross section units are [cm^-2].
    flux_spectra : float
        Flux spectrum [n cm^-2 s^-1]
    chain_file : str
        Path to the depletion chain XML file.
    keff : 2-tuple of float, optional
       keff eigenvalue and uncertainty from transport calculation.
       Default is None.
    fission_q : dict, optional
        Dictionary of nuclides and their fission Q values [eV]. If not given,
        values will be pulled from the ``chain_file``.
    prev_results : Results, optional
        Results from a previous depletion calculation.
    reduce_chain : bool, optional
        If True, use :meth:`openmc.deplete.Chain.reduce` to reduce the
        depletion chain up to ``reduce_chain_level``. Default is False.
    reduce_chain_level : int, optional
        Depth of the search when reducing the depletion chain. Only used
        if ``reduce_chain`` evaluates to true. The default value of
        ``None`` implies no limit on the depth.


    Attributes
    ----------
    round_number : bool
        Whether or not to round output to OpenMC to 8 digits.
        Useful in testing, as OpenMC is incredibly sensitive to exact values.
    prev_res : Results or None
        Results from a previous depletion calculation. ``None`` if no
        results are to be used.
    number : openmc.deplete.AtomNumber
        Total number of atoms in simulation.
    nuclides_with_data : set of str
        A set listing all unique nuclides available from cross_sections.xml.
    chain : openmc.deplete.Chain
        The depletion chain information necessary to form matrices and tallies.
    reaction_rates : openmc.deplete.ReactionRates
        Reaction rates from the last operator step.
    prev_res : Results or None
        Results from a previous depletion calculation. ``None`` if no
        results are to be used.
    """

    @classmethod
    def from_nuclides(cls, volume, nuclides, micro_xs,
                 flux_spectra,
                 chain_file,
                 keff=None,
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
            values. Nuclide concentration units are [atom/cm^3].
        micro_xs : pandas.DataFrame
            DataFrame with nuclides names as index and microscopic cross section
            data in the columns. Cross section units are [cm^-2].
        flux_spectra : float
            Flux spectrum [n cm^-2 s^-1]
        chain_file : str
            Path to the depletion chain XML file.
        keff : 2-tuple of float, optional
           keff eigenvalue and uncertainty from transport calculation.
           Default is None.
        fission_q : dict, optional
            Dictionary of nuclides and their fission Q values [eV]. If not given,
            values will be pulled from the ``chain_file``.
        prev_results : Results, optional
            Results from a previous depletion calculation.
        reduce_chain : bool, optional
            If True, use :meth:`openmc.deplete.Chain.reduce` to reduce the
            depletion chain up to ``reduce_chain_level``. Default is False.
        reduce_chain_level : int, optional
            Depth of the search when reducing the depletion chain. Only used
            if ``reduce_chain`` evaluates to true. The default value of
            ``None`` implies no limit on the depth.
        """
        check_type('nuclides', nuclides, dict, str)
        materials = cls._consolidate_nuclides_to_material(nuclides, volume)
        return cls(materials,
                     micro_xs,
                     flux_spectra,
                     chain_file,
                     keff,
                     fission_q,
                     prev_results,
                     reduce_chain,
                     reduce_chain_level,
                     fission_yield_opts)

    def __init__(self,
                 materials,
                 micro_xs,
                 flux_spectra,
                 chain_file,
                 keff=None,
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
            keff = ufloat(keff)

        self._keff = keff
        self.flux_spectra = flux_spectra

        diff_burnable_mats=False
        helper_kwargs = dict()
        helper_kwargs['fission_yield_opts'] = fission_yield_opts

        super().__init__(
            materials,
            micro_xs,
            chain_file,
            prev_results,
            diff_burnable_mats,
            fission_q,
            0.0,
            helper_kwargs,
            reduce_chain,
            reduce_chain_level)

    @staticmethod
    def _consolidate_nuclides_to_material(nuclides, volume):
        """Puts nuclide list into an openmc.Materials object.

        """
        openmc.reset_auto_ids()
        mat = openmc.Material()
        for nuc, conc in nuclides.items():
            mat.add_nuclide(nuc, conc / 1e24) #convert to at/b-cm

        mat.volume = volume
        mat.depleteable = True

        return openmc.Materials([mat])

    def _differentiate_burnable_mats():
        """Assign distribmats for each burnable material"""
        pass

    def _load_previous_results():
        """Load in results from a previous depletion calculation."""
        pass

    def _get_nuclides_with_data(self, cross_sections):
        """Finds nuclides with cross section data
        """
        return set(cross_sections.index)

    class FluxTimesXSHelper(ReactionRateHelper):
        """Class for generating one-group reaction rates with flux and
        one-group cross sections.

        This class does not generate tallies, and instead stores cross sections
        for each nuclide and transmutation reaction relevant for a depletion
        calculation. The reaction rate is calculated by multiplying the flux by the
        cross sections.

        Parameters
        ----------
        outer : openmc.deplete.FluxDepletionOperator
            Reference to the object encapsulate FluxTimesXSHelper.
            We pass this so we don't have to duplicate the ``number`` object.
        n_nucs : int
            Number of burnable nuclides tracked by :class:`openmc.deplete.Operator`
        n_react : int
            Number of reactions tracked by :class:`openmc.deplete.Operator`

        Attributes
        ----------
        nuc_ind_map : dict of int to str
            Dictionary mapping the nuclide index to nuclide name
        rxn_ind_map : dict of int to str
            Dictionary mapping reaction index to reaction name

        """
        def __init__(self, n_nuc, n_react, outer):
            super().__init__(n_nuc, n_react)
            self.outer = outer
            self.nuc_ind_map = None
            self.rxn_ind_map = None

        def generate_tallies(self, materials, scores):
            """Unused in this case"""

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
            for i, (i_nuc, i_react) in enumerate(product(nuc_index, react_index)):
                nuc = self.nuc_ind_map[i_nuc]
                rxn = self.rxn_ind_map[i_react]
                density = self.outer.number.get_atom_density(mat_id, nuc)
                self._results_cache[i_nuc, i_react] = self.outer.cross_sections[rxn][nuc] * density

            return self._results_cache

    def _get_helper_classes(self, helper_kwargs):
        """Get helper classes for calculating reation rates and fission yields"""

        fission_yield_opts = helper_kwargs['fission_yield_opts']

        rates = self.reaction_rates
        # Get classes to assit working with tallies
        nuc_ind_map = {ind: nuc for nuc, ind in rates.index_nuc.items()}
        rxn_ind_map = {ind: rxn for rxn, ind in rates.index_rx.items()}

        self._rate_helper = self.FluxTimesXSHelper(self.reaction_rates.n_nuc, self.reaction_rates.n_react, self)
        self._rate_helper.nuc_ind_map = nuc_ind_map
        self._rate_helper.rxn_ind_map = rxn_ind_map

        self._normalization_helper = SourceRateHelper()

        # Select and create fission yield helper
        fission_helper = ConstantFissionYieldHelper
        fission_yield_opts = (
            {} if fission_yield_opts is None else fission_yield_opts)
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

        # Use the flux spectra as a "source rate"
        rates = self._calculate_reaction_rates(self.flux_spectra)
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
                                print(
                                    "WARNING: nuclide ",
                                    nuc,
                                    " in material ",
                                    mat,
                                    " is negative (density = ",
                                    val,
                                    " at/barn-cm)")
                            number_i[mat, nuc] = 0.0

                # TODO Update densities on the Python side, otherwise the
                # summary.h5 file contains densities at the first time step

    def write_bos_data(self, step):
        """Document beginning of step data for a given step

        Called at the beginning of a depletion step and at
        the final point in the simulation.

        Parameters
        ----------
        step : int
            Current depletion step including restarts
        """
        # Since we aren't running a transport simulation, we simply pass
        pass


    @staticmethod
    def create_micro_xs_from_data_array(
            nuclides, reactions, data, units='barn'):
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

        FluxDepletionOperator._validate_micro_xs_inputs(nuclides, reactions, data)

        # Convert to cm^2
        if units == 'barn':
            data /= 1e24

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
            micro_xs /= 1e24

        return micro_xs

    # Convenience function for the micro_xs static methods
    @staticmethod
    def _validate_micro_xs_inputs(nuclides, reactions, data):
        check_iterable_type('nuclides', nuclides, str)
        check_iterable_type('reactions', reactions, str)
        check_type('data', data, np.ndarray, expected_iter_type=float)
        for reaction in reactions:
            check_value('reactions', reaction, valid_rxns)



