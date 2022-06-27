"""Pure depletion operator

This module implements a pure depletion operator that uses user provided fluxes
and one-group cross sections.

"""


import copy
from collections import OrderedDict
import os
from warnings import warn

import numpy as np
from uncertainties import ufloat

import openmc
from openmc.checkvalue import check_value
from openmc.data import DataLibrary
from openmc.exceptions import DataError
from openmc.mpi import comm
from .abc import TransportOperator, OperatorResult
from .atom_number import AtomNumber
from .chain import _find_chain_file
from .reaction_rates import ReactionRates
from .results import Results
from .helpers import (
    FluxReactionRateHelper, ChainFissionHelper, SourceRateHelper,
    ConstantFissionYieldHelper)




class FluxSpectraDepletionOperator(TransportOperator):
    """Depletion operator that uses a user provided flux spectrum to
    calculate reaction rates. The flux provided must match the type of cross
    section library in use (contunuous flux for continuous energy library,
    multi-group flux for multi-group library)

    Instances of this class can be used to perform depletion using one group
    cross sections and constant flux. Normally, a user needn't call methods of
    this class directly. Instead, an instance of this class is passed to an
    integrator class, such as :class:`openmc.deplete.CECMIntegrator`

    Parameters
    ----------
    nuclides : dict of str to float
        Dictionary with nuclides names as keys and concentration as values.
    micro_xs : pandas.DataFrame
        DataFrame with nuclides names as index and microscopic cross section
        data in the columns.
    flux_spectra : ???
        Flux spectrum
    chain_file : str
        Path to the depletion chain XML file.
    fission_q : dict, optional
        Dictionary of nuclides and their fission Q values [eV]. If not given,
        values will be pulled from the ``chain_file``.
    dilute_initial : float, optional
        Initial atom density [atoms/cm^3] to add for nuclides that are zero
        in initial condition to ensure they exist in the decay chain.
        Only done for nuclides with reaction rates.
        Defaults to 1.0e3.
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
    dilute_initial : float
        Initial atom density [atoms/cm^3] to add for nuclides that are zero
        in initial condition to ensure they exist in the decay chain.
        Only done for nuclides with reaction rates.
    prev_res : Results or None
        Results from a previous depletion calculation. ``None`` if no
        results are to be used.
    """

    def __init__(self, nuclides, micro_xs, flux_spectra, chain_file, fission_q=None, dilute_initial=1.0e3,
                 prev_results=None, reduce_chain=False, reduce_chain_level=None):


        # TODO : validate nuclides and micro-xs parameters

        super.__init__(chain_file, fission_q, dilute_initial, prev_results)
        self.round_number = False

        self.flux_spectra = flux_spectra

        # Reduce the chain
        if reduce_chain:
            all_nuclides = set()
            for nuclide_name in nuclides.keys():
                all_isotopes.add(nuclide_name)
            self.chain = self.chain.reduce(all_nuclides, reduce_chain_level)

        # Clear out OpenMC, create task lists, distribute
        #self.burnable_mats, volume, nuclides = self.get_burnable_mats()

        # TODO: add support for loading previous results

        # TODO : Determine which nuclides have incident neutron data
        self.nuclides_with_data = self._get_nuclides_with_data()

        # Select nuclides with data that are also in the chain
        self._burnable_nucs = [nuc.name for nuc in self.chain.nuclides
                               if nuc.name in self.nuclides_with_data]

        # Extract number densities from the geometry / previous depletion run
        self._extract_number('0', volume, list(nuclides.keys()), self.prev_res)

        # Create reaction rates array
        self.reaction_rates = ReactionRates(
            '0', self._burnable_nucs, self.chain.reactions)

        # Initialize normalization helper
        if normalization_mode == "fission-q":
            self._normalization_helper = ChainFissionHelper()
        else:
            self._normalization_helper = SourceRateHelper()

        # Select and create fission yield helper
        fission_helper = ConstantFissionYieldHelper
        fission_yield_opts = (
            {} if fission_yield_opts is None else fission_yield_opts)
        self._yield_helper = fission_helper.from_operator(
            self, **fission_yield_opts)


    def __call__(self, vec, source_rate):
        """Runs a simulation.

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

        # Update nuclides data in preparation for transport solve
        nuclides = self._get_reaction_nuclides()
        self._yield_helper.update_tally_nuclides(nuclides)

        rates = self.reaction_rates
        rates.fill(0.0)

        # Get k and uncertainty
        # TODO : add functionality to get this from the transport solver
        keff = ...

        # Form fast map
        nuc_ind = [rates.index_nuc[nuc] for nuc in nuclides]
        react_ind = [rates.index_rx[react] for react in self.chain.reactions]

        # Keep track of energy produced from all reactions in eV per source
        # particle
        self._normalization_helper.reset()

        # Store fission yield dictionaries
        fission_yields = []

        # Create arrays to store fission Q values, reaction rates, and nuclide
        # numbers, zeroed out in material iteration
        number = np.empty(rates.n_nuc)

        fission_ind = rates.index_rx.get("fission")

        # Zero out reaction rates and nuclide numbers
        number.fill(0.0)

        # Get new number densities
        for nuc, i_nuc_results in zip(nuclides, nuc_ind):
            number[i_nuc_results] = self.number[0, nuc]

        # Store microscopic cross sections in rates array
        for nuc in nuclides:
            for rxn in self.chain.reactions:
                rates.set('0', nuc, rxn, self.micro_xs[rxn].loc[nuc])

        # Get reaction rate in reactions/sec
        rates *= self.flux_spectra

        # Compute fission yields for this material
        fission_yields.append(self._yield_helper.weighted_yields(0))

        # Accumulate energy from fission
        if fission_ind is not None:
            self._normalization_helper.update(rxn_rates[:, fission_ind])

        ## These don't seem relevant so I'm commenting them out for now
        ## will delete if they are indeed not useful
        # Divide by total number and store
        #rates[0] = self._rate_helper.divide_by_adens(number)

        # Scale reaction rates to obtain units of reactions/sec
        #rates *= self._normalization_helper.factor(source_rate)
        ##

        # Store new fission yields on the chain
        self.chain.fission_yields = fission_yields

        return OperatorResult(keff, rates)


    def initial_condition(self):
        """Performs final setup and returns initial condition.

        Returns
        -------
        list of numpy.ndarray
            Total density for initial conditions.
        """

        self._normalization_helper.prepare(
            self.chain.nuclides, self.reaction_rates.index_nuc)

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
        ...


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
                            # negative. CRAM does not guarantee positive values.
                            if val < -1.0e-21:
                                print("WARNING: nuclide ", nuc, " in material ", mat,
                                      " is negative (density = ", val, " at/barn-cm)")
                            number_i[mat, nuc] = 0.0


                #TODO Update densities on the Python side, otherwise the
                # summary.h5 file contains densities at the first time step

    def _get_reaction_nuclides(self):
        """Determine nuclides that should be tallied for reaction rates.

        This method returns a list of all nuclides that have neutron data and
        are listed in the depletion chain. Technically, we should tally nuclides
        that may not appear in the depletion chain because we still need to get
        the fission reaction rate for these nuclides in order to normalize
        power, but that is left as a future exercise.

        Returns
        -------
        list of str
            Tally nuclides

        """
        nuc_set = set()

        # Create the set of all nuclides in the decay chain in materials marked
        # for burning in which the number density is greater than zero.
        for nuc in self.number.nuclides:
            if nuc in self.nuclides_with_data:
                if np.sum(self.number[:, nuc]) > 0.0:
                    nuc_set.add(nuc)

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


