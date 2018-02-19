"""The OpenMC wrapper module.

This module implements the depletion -> OpenMC linkage.
"""

import copy
from collections import OrderedDict
from itertools import chain
import os
import random
import sys
import time
try:
    import lxml.etree as ET
    _have_lxml = True
except ImportError:
    import xml.etree.ElementTree as ET
    _have_lxml = False

import h5py
import numpy as np

import openmc
import openmc.capi
from openmc.data import JOULE_PER_EV
from . import comm
from .abc import Settings, Operator, OperatorResult
from .atom_number import AtomNumber
from .chain import Chain
from .reaction_rates import ReactionRates


def _distribute(items):
    min_size, extra = divmod(len(items), comm.size)
    j = 0
    for i in range(comm.size):
        chunk_size = min_size + int(i < extra)
        if comm.rank == i:
            return items[j:j + chunk_size]
        j += chunk_size


class OpenMCSettings(Settings):
    """Extends Settings to provide information OpenMC needs to run.

    Attributes
    ----------
    dt_vec : numpy.array
        Array of time steps to in units of [s]
    output_dir : pathlib.Path
        Path to output directory to save results.
    chain_file : str
        Path to the depletion chain XML file.  Defaults to the
        :envvar:`OPENMC_DEPLETE_CHAIN` environment variable if it exists.
    dilute_initial : float
        Initial atom density to add for nuclides that are zero in initial
        condition to ensure they exist in the decay chain.  Only done for
        nuclides with reaction rates. Defaults to 1.0e3.
    power : float
        Power of the reactor in [W]. For a 2D problem, the power can be given in
        W/cm as long as the "volume" assigned to a depletion material is
        actually an area in cm^2.
    round_number : bool
        Whether or not to round output to OpenMC to 8 digits.
        Useful in testing, as OpenMC is incredibly sensitive to exact values.
    settings : openmc.Settings
        Settings for OpenMC simulations

    """

    _depletion_attrs = {'dt_vec', '_output_dir', 'chain_file', 'dilute_initial',
                        'round_number', 'power'}

    def __init__(self):
        super().__init__()
        self.round_number = False

        # Avoid setattr to create OpenMC settings
        self.__dict__['settings'] = openmc.Settings()

    def __setattr__(self, name, value):
        if hasattr(self.__class__, name):
            # Use properties when appropriate
            prop = getattr(self.__class__, name)
            prop.fset(self, value)
        elif name in self._depletion_attrs:
            # For known attributes, store in dictionary
            self.__dict__[name] = value
        else:
            # otherwise, delegate to openmc.Settings
            setattr(self.__dict__['settings'], name, value)

    def __getattr__(self, name):
        if name in self._depletion_attrs:
            return self.__dict__[name]
        else:
            return getattr(self.__dict__['settings'], name)


class OpenMCOperator(Operator):
    """OpenMC transport operator

    Parameters
    ----------
    geometry : openmc.Geometry
        The OpenMC geometry object.
    settings : openmc.deplete.OpenMCSettings
        Settings object.

    Attributes
    ----------
    settings : OpenMCSettings
        Settings object. (From Operator)
    geometry : openmc.Geometry
        The OpenMC geometry object.
    number : openmc.deplete.AtomNumber
        Total number of atoms in simulation.
    nuclides_with_data : set of str
        A set listing all unique nuclides available from cross_sections.xml.
    chain : openmc.deplete.Chain
        The depletion chain information necessary to form matrices and tallies.
    reaction_rates : openmc.deplete.ReactionRates
        Reaction rates from the last operator step.
    burn_mat_to_ind : OrderedDict of str to int
        Dictionary mapping material ID (as a string) to an index in reaction_rates.
    burnable_mats : list of str
        All burnable material IDs

    """
    def __init__(self, geometry, settings):
        super().__init__(settings)
        self.geometry = geometry

        # Read depletion chain
        self.chain = Chain.from_xml(settings.chain_file)

        # Clear out OpenMC, create task lists, distribute
        openmc.reset_auto_ids()
        self.burnable_mats, volume, nuc_dict = self._get_burnable_mats()
        local_mats = _distribute(self.burnable_mats)

        # Determine which nuclides have incident neutron data
        self.nuclides_with_data = self._get_nuclides_with_data()
        self._burnable_nucs = [nuc for nuc in self.nuclides_with_data
                               if nuc in self.chain]

        # Extract number densities from the geometry
        self._extract_number(local_mats, volume, nuc_dict)

        # Create reaction rates array
        index_rx = {rx: i for i, rx in enumerate(self.chain.reactions)}
        self.reaction_rates = ReactionRates(
            self.burn_mat_to_ind, self._burnable_nucs, index_rx)

    def __call__(self, vec, print_out=True):
        """Runs a simulation.

        Parameters
        ----------
        vec : list of numpy.array
            Total atoms to be used in function.
        print_out : bool, optional
            Whether or not to print out time.

        Returns
        -------
        openmc.deplete.OperatorResult
            Eigenvalue and reaction rates resulting from transport operator

        """
        # Prevent OpenMC from complaining about re-creating tallies
        openmc.reset_auto_ids()

        # Update status
        self.number.set_density(vec)

        time_start = time.time()

        # Update material compositions and tally nuclides
        self._update_materials()
        self._tally.nuclides = self._get_tally_nuclides()

        # Run OpenMC
        openmc.capi.reset()
        openmc.capi.run()

        time_openmc = time.time()

        # Extract results
        op_result = self._unpack_tallies_and_normalize()

        if comm.rank == 0:
            time_unpack = time.time()

            if print_out:
                print("Time to openmc: ", time_openmc - time_start)
                print("Time to unpack: ", time_unpack - time_openmc)

        return copy.deepcopy(op_result)

    def _get_burnable_mats(self):
        """Determine depletable materials, volumes, and nuclids

        Returns
        -------
        burnable_mats : list of str
            List of burnable material IDs
        volume : OrderedDict of str to float
            Volume of each material in [cm^3]
        nuc_dict : OrderedDict of str to int
            Nuclides in order of how they'll appear in the simulation.

        """
        burnable_mats = set()
        model_nuclides = set()
        volume = OrderedDict()

        # Iterate once through the geometry to get dictionaries
        for mat in self.geometry.get_all_materials().values():
            for nuclide in mat.get_nuclides():
                model_nuclides.add(nuclide)
            if mat.depletable:
                burnable_mats.add(str(mat.id))
                if mat.volume is None:
                    raise RuntimeError("Volume not specified for depletable "
                                       "material with ID={}.".format(mat.id))
                volume[str(mat.id)] = mat.volume

        # Sort the sets
        burnable_mats = sorted(burnable_mats, key=int)
        model_nuclides = sorted(model_nuclides)

        # Construct a global nuclide dictionary, burned first
        nuc_dict = copy.deepcopy(self.chain.nuclide_dict)
        i = len(nuc_dict)
        for nuc in model_nuclides:
            if nuc not in nuc_dict:
                nuc_dict[nuc] = i
                i += 1

        return burnable_mats, volume, nuc_dict

    def _extract_number(self, local_mats, volume, nuc_dict):
        """Construct AtomNumber using geometry

        Parameters
        ----------
        local_mats : list of str
            Material IDs to be managed by this process
        volume : OrderedDict of str to float
            Volumes for the above materials in [cm^3]
        nuc_dict : OrderedDict of str to int
            Nuclides to be used in the simulation.

        """
        # Same with materials
        self.burn_mat_to_ind = OrderedDict()
        for i, mat in enumerate(local_mats):
            self.burn_mat_to_ind[mat] = i

        self.number = AtomNumber(self.burn_mat_to_ind, nuc_dict, volume,
                                 len(self.chain))

        if self.settings.dilute_initial != 0.0:
            for nuc in self._burnable_nucs:
                self.number.set_atom_density(np.s_[:], nuc,
                                             self.settings.dilute_initial)

        # Now extract the number densities and store
        for mat in self.geometry.get_all_materials().values():
            if str(mat.id) in self.burn_mat_to_ind:
                self._set_number_from_mat(mat)

    def _set_number_from_mat(self, mat):
        """Extracts material and number densities from openmc.Material

        Parameters
        ----------
        mat : openmc.Material
            The material to read from

        """
        mat_id = str(mat.id)
        mat_ind = self.number.mat_to_ind[mat_id]

        for nuclide, density in mat.get_nuclide_atom_densities().values():
            number = density * 1.0e24
            self.number.set_atom_density(mat_id, nuclide, number)

    def form_matrix(self, y, mat):
        """Forms the depletion matrix.

        Parameters
        ----------
        y : numpy.ndarray
            An array representing reaction rates for this cell.
        mat : int
            Material id.

        Returns
        -------
        scipy.sparse.csr_matrix
            Sparse matrix representing the depletion matrix.
        """

        return copy.deepcopy(self.chain.form_matrix(y[mat, :, :]))

    def initial_condition(self):
        """Performs final setup and returns initial condition.

        Returns
        -------
        list of numpy.array
            Total density for initial conditions.
        """

        # Create XML files
        if comm.rank == 0:
            self.geometry.export_to_xml()
            self.settings.settings.export_to_xml()
            self._generate_materials_xml()

        # Initialize OpenMC library
        comm.barrier()
        openmc.capi.init(comm)

        # Generate tallies in memory
        self._generate_tallies()

        # Return number density vector
        return list(self.number.get_mat_slice(np.s_[:]))

    def finalize(self):
        """Finalize a depletion simulation and release resources."""
        openmc.capi.finalize()

    def _update_materials(self):
        """Updates material compositions in OpenMC on all processes."""

        for rank in range(comm.size):
            number_i = comm.bcast(self.number, root=rank)

            for mat in number_i.mat_to_ind:
                nuclides = []
                densities = []
                for nuc in number_i.nuc_to_ind:
                    if nuc in self.nuclides_with_data:
                        val = 1.0e-24 * number_i.get_atom_density(mat, nuc)

                        # If nuclide is zero, do not add to the problem.
                        if val > 0.0:
                            if self.settings.round_number:
                                val_magnitude = np.floor(np.log10(val))
                                val_scaled = val / 10**val_magnitude
                                val_round = round(val_scaled, 8)

                                val = val_round * 10**val_magnitude

                            nuclides.append(nuc)
                            densities.append(val)
                        else:
                            # Only output warnings if values are significantly
                            # negative.  CRAM does not guarantee positive values.
                            if val < -1.0e-21:
                                print("WARNING: nuclide ", nuc, " in material ", mat,
                                      " is negative (density = ", val, " at/barn-cm)")
                            number_i[mat, nuc] = 0.0

                mat_internal = openmc.capi.materials[int(mat)]
                mat_internal.set_densities(nuclides, densities)

    def _generate_materials_xml(self):
        """Creates materials.xml from self.number.

        Due to uncertainty with how MPI interacts with OpenMC API, this
        constructs the XML manually.  The long term goal is to do this
        through direct memory writing.

        """
        materials = openmc.Materials(self.geometry.get_all_materials()
                                     .values())

        # Sort nuclides according to order in AtomNumber object
        nuclides = list(self.number.nuc_to_ind.keys())
        for mat in materials:
            mat._nuclides.sort(key=lambda x: nuclides.index(x[0]))

        materials.export_to_xml()

    def _get_tally_nuclides(self):
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

        # Create the set of all nuclides in the decay chain in cells marked for
        # burning in which the number density is greater than zero.
        for nuc in self.number.nuc_to_ind:
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
            nuc_list = [nuc for nuc in self.number.nuc_to_ind
                        if nuc in nuc_set]
        else:
            nuc_list = None

        # Store list of tally nuclides on each process
        nuc_list = comm.bcast(nuc_list)
        return [nuc for nuc in nuc_list if nuc in self.chain]

    def _generate_tallies(self):
        """Generates depletion tallies.

        Using information from the depletion chain as well as the nuclides
        currently in the problem, this function automatically generates a
        tally.xml for the simulation.

        """
        # Create tallies for depleting regions
        materials = [openmc.capi.materials[int(i)]
                     for i in self.burnable_mats]
        mat_filter = openmc.capi.MaterialFilter(materials)

        # Set up a tally that has a material filter covering each depletable
        # material and scores corresponding to all reactions that cause
        # transmutation. The nuclides for the tally are set later when eval() is
        # called.
        self._tally = openmc.capi.Tally()
        self._tally.scores = self.chain.reactions
        self._tally.filters = [mat_filter]

    def _unpack_tallies_and_normalize(self):
        """Unpack tallies from OpenMC and return an operator result

        This method uses OpenMC's C API bindings to determine the k-effective
        value and reaction rates from the simulation. The reaction rates are
        normalized by the user-specified power, summing the product of the
        fission reaction rate times the fission Q value for each material.

        Returns
        -------
        openmc.deplete.OperatorResult
            Eigenvalue and reaction rates resulting from transport operator

        """
        rates = self.reaction_rates
        rates[:, :, :] = 0.0

        k_combined = openmc.capi.keff()[0]

        # Extract tally bins
        materials = self.burnable_mats
        nuclides = self._tally.nuclides

        # Form fast map
        nuc_ind = [rates.index_nuc[nuc] for nuc in nuclides]
        react_ind = [rates.index_rx[react] for react in self.chain.reactions]

        # Compute fission power
        # TODO : improve this calculation

        # Keep track of energy produced from all reactions in eV per source
        # particle
        energy = 0.0

        # Create arrays to store fission Q values, reaction rates, and nuclide
        # numbers
        fission_Q = np.zeros(rates.n_nuc)
        rates_expanded = np.zeros((rates.n_nuc, rates.n_react))
        number = np.zeros(rates.n_nuc)

        fission_ind = rates.index_rx["fission"]

        for nuclide in self.chain.nuclides:
            if nuclide.name in rates.index_nuc:
                for rx in nuclide.reactions:
                    if rx.type == 'fission':
                        ind = rates.index_nuc[nuclide.name]
                        fission_Q[ind] = rx.Q
                        break

        # Extract results
        for i, mat in enumerate(self.burn_mat_to_ind):
            # Get tally index
            slab = materials.index(mat)

            # Get material results hyperslab
            results = self._tally.results[slab, :, 1]

            # Zero out reaction rates and nuclide numbers
            rates_expanded[:] = 0.0
            number[:] = 0.0

            # Expand into our memory layout
            j = 0
            for nuc, i_nuc_results in zip(nuclides, nuc_ind):
                number[i_nuc_results] = self.number[mat, nuc]
                for react in react_ind:
                    rates_expanded[i_nuc_results, react] = results[j]
                    j += 1

            # Accumulate energy from fission
            energy += np.dot(rates_expanded[:, fission_ind], fission_Q)

            # Divide by total number and store
            for i_nuc_results in nuc_ind:
                if number[i_nuc_results] != 0.0:
                    for react in react_ind:
                        rates_expanded[i_nuc_results, react] /= number[i_nuc_results]

            rates[i, :, :] = rates_expanded

        # Reduce energy produced from all processes
        energy = comm.allreduce(energy)

        # Determine power in eV/s
        power = self.settings.power / JOULE_PER_EV

        # Scale reaction rates to obtain units of reactions/sec
        rates[:, :, :] *= power / energy

        return OperatorResult(k_combined, rates)

    def _get_nuclides_with_data(self):
        """Loads a cross_sections.xml file to find participating nuclides.

        This allows for nuclides that are important in the decay chain but not
        important neutronically, or have no cross section data.
        """

        # Reads cross_sections.xml to create a dictionary containing
        # participating (burning and not just decaying) nuclides.

        try:
            filename = os.environ["OPENMC_CROSS_SECTIONS"]
        except KeyError:
            filename = None

        nuclides = set()

        try:
            tree = ET.parse(filename)
        except Exception:
            if filename is None:
                msg = "No cross_sections.xml specified in materials."
            else:
                msg = 'Cross section file "{}" is invalid.'.format(filename)
            raise IOError(msg)

        root = tree.getroot()
        for nuclide_node in root.findall('library'):
            mats = nuclide_node.get('materials')
            if not mats:
                continue
            for name in mats.split():
                # Make a burn list of the union of nuclides in cross_sections.xml
                # and nuclides in depletion chain.
                if name not in nuclides:
                    nuclides.add(name)

        return nuclides

    def get_results_info(self):
        """Returns volume list, cell lists, and nuc lists.

        Returns
        -------
        volume : dict of str float
            Volumes corresponding to materials in full_burn_dict
        nuc_list : list of str
            A list of all nuclide names. Used for sorting the simulation.
        burn_list : list of int
            A list of all cell IDs to be burned.  Used for sorting the simulation.
        full_burn_list : list
            List of all burnable material IDs

        """
        nuc_list = self.number.burn_nuc_list
        burn_list = list(self.burn_mat_to_ind)

        volume = {}
        for i, mat in enumerate(burn_list):
            volume[mat] = self.number.volume[i]

        # Combine volume dictionaries across processes
        volume_list = comm.allgather(volume)
        volume = {k: v for d in volume_list for k, v in d.items()}

        return volume, nuc_list, burn_list, self.burnable_mats
