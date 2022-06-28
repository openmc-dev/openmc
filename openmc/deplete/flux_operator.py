"""Pure depletion operator

This module implements a pure depletion operator that uses user provided fluxes
and one-group cross sections.

"""


import copy
from collections import OrderedDict
import os
from warnings import warn

import numpy as np
import pandas as pd
from uncertainties import ufloat

import openmc
from openmc.checkvalue import check_type
from openmc.data import DataLibrary
from openmc.exceptions import DataError
from openmc.mpi import comm
from .abc import TransportOperator, OperatorResult
from .atom_number import AtomNumber
from .chain import _find_chain_file
from .reaction_rates import ReactionRates
from .results import Results
from .helpers import ConstantFissionYieldHelper




class FluxSpectraDepletionOperator(TransportOperator):
    """Depletion operator that uses a user provided flux spectrum and one-group
    cross sections calculate reaction rates.

    Instances of this class can be used to perform depletion using one group
    cross sections and constant flux. Normally, a user needn't call methods of
    this class directly. Instead, an instance of this class is passed to an
    integrator class, such as :class:`openmc.deplete.CECMIntegrator`

    Parameters
    ----------
    volume : float
        Volume of the material being depleted in cm^3
    nuclides : dict of str to float
        Dictionary with nuclides names as keys and concentration as values.
        Nuclide concentration is assumed to be in [at/cm^3]
    micro_xs : pandas.DataFrame
        DataFrame with nuclides names as index and microscopic cross section
        data in the columns.
    flux_spectra : float
        Flux spectrum [n cm^-2 s^-1]
    chain_file : str
        Path to the depletion chain XML file.
    keff : 2-tuple of float, optional
       keff eigenvalue and uncertainty from transport calculation.
       Defualt is None.
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
    number : openmc.deplete.AtomNumber
        Total number of atoms in simulation.
    nuclides_with_data : set of str
        A set listing all unique nuclides available from cross_sections.xml.
    chain : openmc.deplete.Chain
        The depletion chain information necessary to form matrices and tallies.
    reaction_rates : openmc.deplete.ReactionRates
        Reaction rates from the last operator step.
    heavy_metal : float
        Initial heavy metal inventory [g]
    prev_res : Results or None
        Results from a previous depletion calculation. ``None`` if no
        results are to be used.
    """

    def __init__(self, volume, nuclides, micro_xs, flux_spectra, chain_file, keff=None, fission_q=None, dilute_initial=1.0e3,
                 prev_results=None, reduce_chain=False, reduce_chain_level=None, fission_yield_opts=None):
        super().__init__(chain_file, fission_q, dilute_initial, prev_results)
        self.round_number = False
        self.flux_spectra = flux_spectra
        self._init_nuclides = nuclides
        self._volume = volume

        # Reduce the chain
        if reduce_chain:
            all_nuclides = set()
            for nuclide_name in nuclides.keys():
                all_isotopes.add(nuclide_name)
            self.chain = self.chain.reduce(all_nuclides, reduce_chain_level)

        # Validate nuclides and micro-xs parameters
        check_type('nuclides', nuclides, dict, str)
        check_type('micro_xs', micro_xs, pd.DataFrame)

        self._micro_xs = micro_xs
        self._keff = keff

        nuclides = self._get_all_depletion_nuclides(nuclides)

        # TODO: add support for loading previous results

        # Determine which nuclides have cross section data
        self.nuclides_with_data = self._get_nuclides_with_data()

        # Select nuclides with data that are also in the chain
        self._burnable_nucs = [nuc.name for nuc in self.chain.nuclides
                               if nuc.name in self.nuclides_with_data]

        # TODO : implement _extract_number
        # Extract number densities from the geometry / previous depletion run
        #!!! consider removing this and just using the nuclides dict...
        self._extract_number(['0'], {'0': volume}, nuclides, self.prev_res)

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
        #!!! See previous excamation point comment
        self.number.set_density(vec)
        ## TODO : make sure this function works w the current structure.
        self._update_materials()

        # Update nuclides data in preparation for transport solve

        ## TODO : make sure this function works w the current structure.
        # this nuclides varibale is all the nuclides that will show up via depletion
        # that have xs data
        nuclides = self._get_reaction_nuclides()
        self._yield_helper.update_tally_nuclides(nuclides)

        rates = self.reaction_rates
        rates.fill(0.0)

        # Form fast map
        nuc_ind = [rates.index_nuc[nuc] for nuc in nuclides]
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
        for nuc, i_nuc_results in zip(nuclides, nuc_ind):
            number[i_nuc_results] = self.number[0, nuc]

        # Calculate macroscopic cross sections and store them in rates array
        for nuc in nuclides:
            density = number.get_atom_density('0', nuc)
            for rxn in self.chain.reactions:
                rates.set('0', nuc, rxn, self._micro_xs[rxn].loc[nuc] * density)

        # Get reaction rate in reactions/sec
        rates *= self.flux_spectra

        # Compute fission yields for this material
        fission_yields.append(self._yield_helper.weighted_yields(0))

        # Accumulate energy from fission
        #if fission_ind is not None:
        #    self._normalization_helper.update(rxn_rates[:, fission_ind])

        # Divide by total number and store
        # the reason we do this is bc in the equation, we multiply the depletion matrix
        # by the nuclide vector. Since what we want is the depletion matrix, we need to
        # divide the reaction rates by the number of atoms to get the right units.
        mask = nonzero(number)
        results = rates[0]
        for col in range(results.shape[1]):
            results[mask, col] /= number[mask]
        rates[0] = results

        # Scale reaction rates to obtain units of reactions/sec
        #rates *= self._normalization_helper.factor(source_rate)
        ##

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

        #self._normalization_helper.prepare(
        #    self.chain.nuclides, self.reaction_rates.index_nuc)

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

    def _get_all_depletion_nuclides(self, nuclides):
        """Determine nuclides that will show up in simulation

        Parameters
        ----------
        nuclides : dict
            dict whose keys are nuclide symbols as strings,
            and whose values are nuclide concentrations in at/cm^3

        Returns
        -------
        nuclides : list of str
            Nuclides in order of how they'll appear in the simulation.

        """

        model_nuclides = set()

        #self.heavy_metal = 0.0

        #bu_mat = openmc.Material()
        #bu_mat.volume = volume
        for n in nuclides:
            model_nuclides.add(n)

        #self.heavy_metal = bu_mat.fissionable_mass

        # Sort the sets
        model_nuclides = sorted(model_nuclides)

        # Construct a global nuclide dictionary, burned first
        nuclides = list(self.chain.nuclide_dict)
        for nuc in model_nuclides:
            if nuc not in nuclides:
                nuclides.append(nuc)

        return nuclides

    def _get_nuclides_with_data(self):
        """Finds nuclides with cross section data"""
        nuclides = set()
        for name in self._micro_xs.index:
            nuclides.add(name)

        return nuclides

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

        if self.dilute_initial != 0.0:
            for nuc in self._burnable_nucs:
                self.number.set_atom_density(np.s_[:], nuc, self.dilute_initial)

        # Now extract and store the number densities
        # From the geometry if no previous depletion results
        if prev_res is None:
            for nuclide in nuclides:
                if nuclide in self._init_nuclides:
                    self.number.set_atom_density('0', nuclide, self._init_nuclides[nuclide])
                elif nuclide not in self._burnable_nucs:
                    self.number.set_atom_density('0', nuclide, 0)

        # Else from previous depletion results
        else:
           raise RuntimeError("Loading from previous results not yet supported")
