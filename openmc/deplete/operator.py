"""OpenMC transport operator

This module implements a transport operator for OpenMC so that it can be used by
depletion integrators. The implementation makes use of the Python bindings to
OpenMC's C API so that reading tally results and updating material number
densities is all done in-memory instead of through the filesystem.

"""

import copy
from collections import OrderedDict
from itertools import chain
import os
import time
import xml.etree.ElementTree as ET

import h5py
import numpy as np

import openmc
import openmc.capi
from openmc.data import JOULE_PER_EV
from . import comm
from .abc import TransportOperator, OperatorResult
from .atom_number import AtomNumber
from .reaction_rates import ReactionRates


def _distribute(items):
    """Distribute items across MPI communicator

    Parameters
    ----------
    items : list
        List of items of distribute

    Returns
    -------
    list
        Items assigned to process that called

    """
    min_size, extra = divmod(len(items), comm.size)
    j = 0
    for i in range(comm.size):
        chunk_size = min_size + int(i < extra)
        if comm.rank == i:
            return items[j:j + chunk_size]
        j += chunk_size


class Operator(TransportOperator):
    """OpenMC transport operator for depletion.

    Instances of this class can be used to perform depletion using OpenMC as the
    transport operator. Normally, a user needn't call methods of this class
    directly. Instead, an instance of this class is passed to an integrator
    function, such as :func:`openmc.deplete.integrator.cecm`.

    Parameters
    ----------
    geometry : openmc.Geometry
        OpenMC geometry object
    settings : openmc.Settings
        OpenMC Settings object
    chain_file : str, optional
        Path to the depletion chain XML file.  Defaults to the
        :envvar:`OPENMC_DEPLETE_CHAIN` environment variable if it exists.

    Attributes
    ----------
    geometry : openmc.Geometry
        OpenMC geometry object
    settings : openmc.Settings
        OpenMC settings object
    dilute_initial : float
        Initial atom density to add for nuclides that are zero in initial
        condition to ensure they exist in the decay chain.  Only done for
        nuclides with reaction rates. Defaults to 1.0e3.
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
    local_mats : list of str
        All burnable material IDs being managed by a single process

    """
    def __init__(self, geometry, settings, chain_file=None):
        super().__init__(chain_file)
        self.round_number = False
        self.settings = settings
        self.geometry = geometry

        # Clear out OpenMC, create task lists, distribute
        openmc.reset_auto_ids()
        self.burnable_mats, volume, nuclides = self._get_burnable_mats()
        self.local_mats = _distribute(self.burnable_mats)

        # Determine which nuclides have incident neutron data
        self.nuclides_with_data = self._get_nuclides_with_data()
        self._burnable_nucs = [nuc for nuc in self.nuclides_with_data
                               if nuc in self.chain]

        # Extract number densities from the geometry
        self._extract_number(self.local_mats, volume, nuclides)

        # Create reaction rates array
        self.reaction_rates = ReactionRates(
            self.local_mats, self._burnable_nucs, self.chain.reactions)

    def __call__(self, vec, power, print_out=True):
        """Runs a simulation.

        Parameters
        ----------
        vec : list of numpy.ndarray
            Total atoms to be used in function.
        power : float
            Power of the reactor in [W]
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
        op_result = self._unpack_tallies_and_normalize(power)

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
        nuclides : list of str
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

        # Make sure there are burnable materials
        if not burnable_mats:
            raise RuntimeError(
                "No depletable materials were found in the model.")

        # Sort the sets
        burnable_mats = sorted(burnable_mats, key=int)
        model_nuclides = sorted(model_nuclides)

        # Construct a global nuclide dictionary, burned first
        nuclides = list(self.chain.nuclide_dict)
        for nuc in model_nuclides:
            if nuc not in nuclides:
                nuclides.append(nuc)

        return burnable_mats, volume, nuclides

    def _extract_number(self, local_mats, volume, nuclides):
        """Construct AtomNumber using geometry

        Parameters
        ----------
        local_mats : list of str
            Material IDs to be managed by this process
        volume : OrderedDict of str to float
            Volumes for the above materials in [cm^3]
        nuclides : list of str
            Nuclides to be used in the simulation.

        """
        self.number = AtomNumber(local_mats, nuclides, volume, len(self.chain))

        if self.dilute_initial != 0.0:
            for nuc in self._burnable_nucs:
                self.number.set_atom_density(np.s_[:], nuc, self.dilute_initial)

        # Now extract the number densities and store
        for mat in self.geometry.get_all_materials().values():
            if str(mat.id) in local_mats:
                self._set_number_from_mat(mat)

    def _set_number_from_mat(self, mat):
        """Extracts material and number densities from openmc.Material

        Parameters
        ----------
        mat : openmc.Material
            The material to read from

        """
        mat_id = str(mat.id)

        for nuclide, density in mat.get_nuclide_atom_densities().values():
            number = density * 1.0e24
            self.number.set_atom_density(mat_id, nuclide, number)

    def initial_condition(self):
        """Performs final setup and returns initial condition.

        Returns
        -------
        list of numpy.ndarray
            Total density for initial conditions.
        """

        # Create XML files
        if comm.rank == 0:
            self.geometry.export_to_xml()
            self.settings.export_to_xml()
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
        nuclides = list(self.number.nuclides)
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

    def _unpack_tallies_and_normalize(self, power):
        """Unpack tallies from OpenMC and return an operator result

        This method uses OpenMC's C API bindings to determine the k-effective
        value and reaction rates from the simulation. The reaction rates are
        normalized by the user-specified power, summing the product of the
        fission reaction rate times the fission Q value for each material.

        Parameters
        ----------
        power : float
            Power of the reactor in [W]

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
        for i, mat in enumerate(self.local_mats):
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
        power /= JOULE_PER_EV

        # Scale reaction rates to obtain units of reactions/sec
        rates *= power / energy

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
        """Returns volume list, material lists, and nuc lists.

        Returns
        -------
        volume : dict of str float
            Volumes corresponding to materials in full_burn_dict
        nuc_list : list of str
            A list of all nuclide names. Used for sorting the simulation.
        burn_list : list of int
            A list of all material IDs to be burned.  Used for sorting the simulation.
        full_burn_list : list
            List of all burnable material IDs

        """
        nuc_list = self.number.burnable_nuclides
        burn_list = self.local_mats

        volume = {}
        for i, mat in enumerate(burn_list):
            volume[mat] = self.number.volume[i]

        # Combine volume dictionaries across processes
        volume_list = comm.allgather(volume)
        volume = {k: v for d in volume_list for k, v in d.items()}

        return volume, nuc_list, burn_list, self.burnable_mats
