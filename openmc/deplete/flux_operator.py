"""Pure depletion operator

This module implements a pure depletion operator that uses user provided fluxes
and one-group cross sections.

"""


from warnings import warn

import numpy as np
import pandas as pd
from uncertainties import ufloat

from openmc.checkvalue import check_type, check_value, check_iterable_type
from openmc.mpi import comm
from .abc import TransportOperator, OperatorResult
from .atom_number import AtomNumber
from .chain import REACTIONS
from .reaction_rates import ReactionRates
from .helpers import ConstantFissionYieldHelper

valid_rxns = list(REACTIONS)
valid_rxns.append('fission')



class FluxDepletionOperator(TransportOperator):
    """Depletion operator that uses a user-provided flux spectrum and one-group
    cross sections to calculate reaction rates.

    Instances of this class can be used to perform depletion using one group
    cross sections and constant flux. Normally, a user needn't call methods of
    this class directly. Instead, an instance of this class is passed to an
    integrator class, such as :class:`openmc.deplete.CECMIntegrator`

    Parameters
    ----------
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

    def __init__(self,
                 volume,
                 nuclides,
                 micro_xs,
                 flux_spectra,
                 chain_file,
                 keff=None,
                 fission_q=None,
                 prev_results=None,
                 reduce_chain=False,
                 reduce_chain_level=None,
                 fission_yield_opts=None):
        super().__init__(chain_file, fission_q, 0.0, prev_results)
        self.round_number = False
        self.flux_spectra = flux_spectra
        self._init_nuclides = nuclides
        self._volume = volume

        # Reduce the chain to only those nuclides present
        if reduce_chain:
            init_nuc_names = set(nuclides.keys())
            self.chain = self.chain.reduce(init_nuc_names, reduce_chain_level)

        # Validate nuclides and micro-xs parameters
        check_type('nuclides', nuclides, dict, str)
        check_type('micro_xs', micro_xs, pd.DataFrame)

        self._micro_xs = micro_xs
        if keff is not None:
            check_type('keff', keff, tuple, float)
            keff = ufloat(keff)

        self._keff = keff

        self._all_nuclides = self._get_all_nuclides_in_simulation()

        # TODO: add support for loading previous results

        # Determine which nuclides have cross section data
        # This nuclides variables contains every nuclides
        # for which there is an entry in the micro_xs parameter
        self.nuclides_with_data = self._get_nuclides_with_data()

        # Select nuclides with data that are also in the chain
        self._burnable_nucs = [nuc.name for nuc in self.chain.nuclides
                               if nuc.name in self.nuclides_with_data]

        # Extract number densities from the geometry / previous depletion run
        self._extract_number(['0'],
                             {'0': volume},
                             self._all_nuclides,
                             self.prev_res)

        # Create reaction rates array
        self.reaction_rates = ReactionRates(
            ['0'], self._burnable_nucs, self.chain.reactions)

        # Select and create fission yield helper
        fission_helper = ConstantFissionYieldHelper
        fission_yield_opts = (
            {} if fission_yield_opts is None else fission_yield_opts)
        self._yield_helper = fission_helper.from_operator(
            self, **fission_yield_opts)

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

        # Update the number densities regardless of the source rate
        self.number.set_density(vec)
        self._update_materials()

        # Get all nuclides for which we will calculate reaction rates
        rxn_nuclides = self._get_reaction_nuclides()

        rates = self.reaction_rates
        rates.fill(0.0)

        # Form fast map
        nuc_ind = [rates.index_nuc[nuc] for nuc in rxn_nuclides]
        react_ind = [rates.index_rx[react] for react in self.chain.reactions]

        # Store fission yield dictionaries
        fission_yields = []

        # Create arrays to store fission Q values, reaction rates, and nuclide
        # numbers, zeroed out in material iteration
        number = np.empty(rates.n_nuc)

        fission_ind = rates.index_rx.get("fission")

        # Zero out reaction rates and nuclide numbers
        number.fill(0.0)

        # Get new number densities
        for nuc, i_nuc_results in zip(rxn_nuclides, nuc_ind):
            number[i_nuc_results] = self.number[0, nuc]

        # Calculate macroscopic cross sections and store them in rates array
        for nuc in rxn_nuclides:
            density = self.number.get_atom_density('0', nuc)
            for rxn in self.chain.reactions:
                rates.set(
                    '0',
                    nuc,
                    rxn,
                    self._micro_xs[rxn, rxn] * density)

        # Get reaction rate in reactions/sec
        rates *= self.flux_spectra

        # Compute fission yields for this material
        fission_yields.append(self._yield_helper.weighted_yields(0))

        # Divide by total number of atoms and store
        # the reason we do this is based on the mathematical equation;
        # in the equation, we multiply the depletion matrix by the nuclide
        # vector. Since what we want is the depletion matrix, we need to
        # divide the reaction rates by the number of atoms to get the right
        # units.
        mask = np.nonzero(number)
        results = rates[0]
        for col in range(results.shape[1]):
            results[mask, col] /= number[mask]
        rates[0] = results

        # Store new fission yields on the chain
        self.chain.fission_yields = fission_yields

        return OperatorResult(self._keff, rates)

    def initial_condition(self):
        """Performs final setup and returns initial condition.

        Returns
        -------
        list of numpy.ndarray
            Total density for initial conditions.
        """

        # Return number density vector
        return list(self.number.get_mat_slice(np.s_[:]))

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

    def get_results_info(self):
        """Returns volume list, cell lists, and nuc lists.

        Returns
        -------
        volume : dict of str to float
            Volumes corresponding to materials in burn_list
        nuc_list : list of str
            A list of all nuclide names. Used for sorting the simulation.
        burn_list : list of int
            A list of all cell IDs to be burned.  Used for sorting the
            simulation.
        full_burn_list : list of int
            All burnable materials in the geometry.
        """
        nuc_list = self.number.burnable_nuclides
        burn_list = ['0']

        volume = {'0': self._volume}

        # Combine volume dictionaries across processes
        volume_list = comm.allgather(volume)
        volume = {k: v for d in volume_list for k, v in d.items()}

        return volume, nuc_list, burn_list, burn_list

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

    def _get_reaction_nuclides(self):
        """Determine nuclides that should have reaction rates

        This method returns a list of all nuclides that have neutron data and
        are listed in the depletion chain. Technically, we should list nuclides
        that may not appear in the depletion chain because we still need to get
        the fission reaction rate for these nuclides in order to normalize
        power, but that is left as a future exercise.

        Returns
        -------
        list of str
            nuclides with reaction rates

        """
        # Create the set of all nuclides in the decay chain in materials marked
        # for burning in which the number density is greater than zero.
        nuc_set = set(self._all_nuclides).intersection(self.nuclides_with_data)
        for nuc in nuc_set:
            if not np.sum(self.number[:, nuc]) > 0.0:
                nuc_set.remove(nuc)

        # Communicate which nuclides have nonzeros to rank 0
        if comm.rank == 0:
            for i in range(1, comm.size):
                nuc_newset = comm.recv(source=i, tag=i)
                nuc_set |= nuc_newset
        else:
            comm.send(nuc_set, dest=0, tag=comm.rank)

        if comm.rank == 0:
            # Sort nuclides in the same order as self.number
            nuc_list = [nuc for nuc in self.number.nuclides
                        if nuc in nuc_set]
        else:
            nuc_list = None

        # Store list of tally nuclides on each process
        nuc_list = comm.bcast(nuc_list)
        return [nuc for nuc in nuc_list if nuc in self.chain]

    def _get_all_nuclides_in_simulation(self):
        """Determine nuclides that will show up in the depletion matrix.
        This is the union of the nuclides provided by the user and
        the nuclides present in the depletion chain.

        Returns
        -------
        all_nuclides : list of str
            Nuclides in order of how they'll appear in the simulation.

        """

        init_nuclides = sorted(self._init_nuclides.keys())

        # Construct a global nuclide dictionary, burned first
        all_nuclides = list(self.chain.nuclide_dict)
        for nuc in init_nuclides:
            if nuc not in all_nuclides:
                all_nuclides.append(nuc)

        return all_nuclides

    def _get_nuclides_with_data(self):
        """Finds nuclides with cross section data"""
        return set(self._micro_xs.index)


    def _extract_number(self, local_mats, volume, nuclides, prev_res=None):
        """Construct AtomNumber using geometry

        Parameters
        ----------
        local_mats : list of str
            Material IDs to be managed by this process
        volume : OrderedDict of str to float
            Volumes for the above materials in [cm^3]
        nuclides : list of str
            Nuclides to be used in the simulation.
        prev_res : Results, optional
            Results from a previous depletion calculation

        """
        self.number = AtomNumber(local_mats, nuclides, volume, len(self.chain))

        # Now extract and store the number densities
        # From the geometry if no previous depletion results
        if prev_res is None:
            for nuclide in nuclides:
                if nuclide in self._init_nuclides:
                    self.number.set_atom_density(
                        '0', nuclide, self._init_nuclides[nuclide])
                elif nuclide not in self._burnable_nucs:
                    self.number.set_atom_density('0', nuclide, 0)

        # Else from previous depletion results
        else:
            raise RuntimeError(
                "Loading from previous results not yet supported")
