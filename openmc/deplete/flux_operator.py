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
from .operator import Operator, _distribute, _find_cross_sections
from .abc import TransportOperator, OperatorResult
from .atom_number import AtomNumber
from .chain import _find_chain_file
from .reaction_rates import ReactionRates
from .results import Results



class FluxSpectraDepletionOperator(TransportOperator, Operator):
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
    nuclides : pandas.DataFrame
        DataFrame contaning nuclide concentration as well as cross section data.
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

    def __init__(self, nuclides, flux_spectra, chain_file, fission_q=None, dilute_initial=1.0e3,
                 prev_results=None, reduce_chain=False, reduce_chain_level=None):

        # Determine cross sections / depletion chain
        # TODO : add support for mg cross sections to _find_cross_sections
        if chain_file is None:
            chain_file = _find_chain_file(cross_sections)

        TransportOperator.__init__(chain_file, fission_q, dilute_initial, prev_results)
        self.round_number = False

        self.flux_spectra = flux_spectra

        # Reduce the chain before we create more materials
        if reduce_chain:
            all_isotopes = set()
            for material in self.materials:
                if not material.depletable:
                    continue
                for name, _dens_percent, _dens_type in material.nuclides:
                    all_isotopes.add(name)
            self.chain = self.chain.reduce(all_isotopes, reduce_chain_level)

        # Clear out OpenMC, create task lists, distribute
        self.burnable_mats, volume, nuclides = Operator._get_burnable_mats()
        self.local_mats = _distribute(self.burnable_mats)

        # Generate map from local materials => material index
        self._mat_index_map = {
            lm: self.burnable_mats.index(lm) for lm in self.local_mats}

        # TODO: add support for loading previous results

        # Determine which nuclides have incident neutron data
        # TODO : add support for mg cross sections to Operator._get_nuclides_with_data
        self.nuclides_with_data = Operator._get_nuclides_with_data(mg_cross_sections)

        # Select nuclides with data that are also in the chain
        self._burnable_nucs = [nuc.name for nuc in self.chain.nuclides
                               if nuc.name in self.nuclides_with_data]

        # Extract number densities from the geometry / previous depletion run
        Operator._extract_number(self.local_mats, volume, nuclides, self.prev_res)

        # Create reaction rates array
        self.reaction_rates = ReactionRates(
            self.local_mats, self._burnable_nucs, self.chain.reactions)

        # TODO : set up rate, yield, and normalization helper objects


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

        # Update tally nuclides data in preparation for transport solve
        nuclides = Operator._get_tally_nuclides()
        self._rate_helper.nuclides = nuclides
        self._normalization_helper.nuclides = nuclides
        self._yield_helper.update_tally_nuclides(nuclides)



        rates = self.reaction_rates
        rates.fill(0.0)

        # Get k and uncertainty
        # TODO : add functionality to get this from the transport solver

        # Form fast map
        nuc_ind = [rates.index_nuc[nuc] for nuc in nuclides]
        react_ind = [rates.index_rx[react] for react in self.chain.reactions]

        # Create arrays to store fission Q values, reaction rates, and nuclide
        # numbers, zeroed out in material iteration
        number = np.empty(rates.n_nuc)

        # Extract results
        for i, mat in enumerate(self.local_mats):
            # Get tally index
            mat_index = self._mat_index_map[mat]

            # TODO : add machinery to multiply
            # each cross section by the corresponding flux

            ...


        # TODO: add flow control to determine if we need to scale
        # the reaction rates
        # Scale reaction rates to obtain units of reactions/sec
        rates *= self._normalization_helper.factor(source_rate)

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

        # TODO : implement these functions so they can store the
        # cross section data
        materials = [openmc.lib.materials[int(i)]
                     for i in self.burnable_mats]
        self._rate_helper.generate_tallies(materials, self.chain.reactions)
        self._normalization_helper.prepare(
            self.chain.nuclides, self.reaction_rates.index_nuc)
        # Tell fission yield helper what materials this process is
        # responsible for
        self._yield_helper.generate_tallies(
            materials, tuple(sorted(self._mat_index_map.values())))


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
