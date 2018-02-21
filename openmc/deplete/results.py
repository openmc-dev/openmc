"""The results module.

Contains results generation and saving capabilities.
"""

from collections import OrderedDict
import copy

import numpy as np
import h5py

from . import comm, have_mpi
from .reaction_rates import ReactionRates

_VERSION_RESULTS = (1, 0)


class Results(object):
    """Output of a depletion run

    Attributes
    ----------
    k : list of float
        Eigenvalue for each substep.
    time : list of float
        Time at beginning, end of step, in seconds.
    n_mat : int
        Number of mats.
    n_nuc : int
        Number of nuclides.
    rates : list of ReactionRates
        The reaction rates for each substep.
    volume : OrderedDict of int to float
        Dictionary mapping mat id to volume.
    mat_to_ind : OrderedDict of str to int
        A dictionary mapping mat ID as string to index.
    nuc_to_ind : OrderedDict of str to int
        A dictionary mapping nuclide name as string to index.
    mat_to_hdf5_ind : OrderedDict of str to int
        A dictionary mapping mat ID as string to global index.
    n_hdf5_mats : int
        Number of materials in entire geometry.
    n_stages : int
        Number of stages in simulation.
    data : numpy.ndarray
        Atom quantity, stored by stage, mat, then by nuclide.

    """
    def __init__(self):
        self.k = None
        self.time = None
        self.rates = None
        self.volume = None

        self.mat_to_ind = None
        self.nuc_to_ind = None
        self.mat_to_hdf5_ind = None

        self.data = None

    def __getitem__(self, pos):
        """Retrieves an item from results.

        Parameters
        ----------
        pos : tuple
            A three-length tuple containing a stage index, mat index and a nuc
            index.  All can be integers or slices.  The second two can be
            strings corresponding to their respective dictionary.

        Returns
        -------
        float
            The atoms for stage, mat, nuc

        """
        stage, mat, nuc = pos
        if isinstance(mat, str):
            mat = self.mat_to_ind[mat]
        if isinstance(nuc, str):
            nuc = self.nuc_to_ind[nuc]

        return self.data[stage, mat, nuc]

    def __setitem__(self, pos, val):
        """Sets an item from results.

        Parameters
        ----------
        pos : tuple
            A three-length tuple containing a stage index, mat index and a nuc
            index.  All can be integers or slices.  The second two can be
            strings corresponding to their respective dictionary.

        val : float
            The value to set data to.

        """
        stage, mat, nuc = pos
        if isinstance(mat, str):
            mat = self.mat_to_ind[mat]
        if isinstance(nuc, str):
            nuc = self.nuc_to_ind[nuc]

        self.data[stage, mat, nuc] = val

    @property
    def n_mat(self):
        return len(self.mat_to_ind)

    @property
    def n_nuc(self):
        return len(self.nuc_to_ind)

    @property
    def n_hdf5_mats(self):
        return len(self.mat_to_hdf5_ind)

    @property
    def n_stages(self):
        return self.data.shape[0]

    def allocate(self, volume, nuc_list, burn_list, full_burn_list, stages):
        """Allocates memory of Results.

        Parameters
        ----------
        volume : dict of str float
            Volumes corresponding to materials in full_burn_dict
        nuc_list : list of str
            A list of all nuclide names. Used for sorting the simulation.
        burn_list : list of int
            A list of all mat IDs to be burned.  Used for sorting the simulation.
        full_burn_list : list of str
            List of all burnable material IDs
        stages : int
            Number of stages in simulation.

        """
        self.volume = copy.deepcopy(volume)
        self.nuc_to_ind = {nuc: i for i, nuc in enumerate(nuc_list)}
        self.mat_to_ind = {mat: i for i, mat in enumerate(burn_list)}
        self.mat_to_hdf5_ind = {mat: i for i, mat in enumerate(full_burn_list)}

        # Create storage array
        self.data = np.zeros((stages, self.n_mat, self.n_nuc))

    def export_to_hdf5(self, filename, step):
        """Export results to an HDF5 file

        Parameters
        ----------
        filename : str
            The filename to write to
        step : int
            What step is this?

        """
        if have_mpi and h5py.get_config().mpi:
            kwargs = {'driver': 'mpio', 'comm': comm}
        else:
            kwargs = {}

        kwargs['mode'] = "w" if step == 0 else "a"

        with h5py.File(filename, **kwargs) as handle:
            self._to_hdf5(handle, step)

    def _write_hdf5_metadata(self, handle):
        """Writes result metadata in HDF5 file

        Parameters
        ----------
        handle : h5py.File or h5py.Group
            An hdf5 file or group type to store this in.

        """
        # Create and save the 5 dictionaries:
        # quantities
        #   self.mat_to_ind -> self.volume (TODO: support for changing volumes)
        #   self.nuc_to_ind
        # reactions
        #   self.rates[0].nuc_to_ind (can be different from above, above is superset)
        #   self.rates[0].react_to_ind
        # these are shared by every step of the simulation, and should be deduplicated.

        # Store concentration mat and nuclide dictionaries (along with volumes)

        handle.attrs['version'] = np.array(_VERSION_RESULTS)
        handle.attrs['filetype'] = np.string_('depletion results')

        mat_list = sorted(self.mat_to_hdf5_ind, key=int)
        nuc_list = sorted(self.nuc_to_ind)
        rxn_list = sorted(self.rates[0].index_rx)

        n_mats = self.n_hdf5_mats
        n_nuc_number = len(nuc_list)
        n_nuc_rxn = len(self.rates[0].index_nuc)
        n_rxn = len(rxn_list)
        n_stages = self.n_stages

        mat_group = handle.create_group("materials")

        for mat in mat_list:
            mat_single_group = mat_group.create_group(mat)
            mat_single_group.attrs["index"] = self.mat_to_hdf5_ind[mat]
            mat_single_group.attrs["volume"] = self.volume[mat]

        nuc_group = handle.create_group("nuclides")

        for nuc in nuc_list:
            nuc_single_group = nuc_group.create_group(nuc)
            nuc_single_group.attrs["atom number index"] = self.nuc_to_ind[nuc]
            if nuc in self.rates[0].index_nuc:
                nuc_single_group.attrs["reaction rate index"] = self.rates[0].index_nuc[nuc]

        rxn_group = handle.create_group("reactions")

        for rxn in rxn_list:
            rxn_single_group = rxn_group.create_group(rxn)
            rxn_single_group.attrs["index"] = self.rates[0].index_rx[rxn]

        # Construct array storage

        handle.create_dataset("number", (1, n_stages, n_mats, n_nuc_number),
                              maxshape=(None, n_stages, n_mats, n_nuc_number),
                              chunks=(1, 1, n_mats, n_nuc_number),
                              dtype='float64')

        handle.create_dataset("reaction rates", (1, n_stages, n_mats, n_nuc_rxn, n_rxn),
                              maxshape=(None, n_stages, n_mats, n_nuc_rxn, n_rxn),
                              chunks=(1, 1, n_mats, n_nuc_rxn, n_rxn),
                              dtype='float64')

        handle.create_dataset("eigenvalues", (1, n_stages),
                              maxshape=(None, n_stages), dtype='float64')

        handle.create_dataset("time", (1, 2), maxshape=(None, 2), dtype='float64')

    def _to_hdf5(self, handle, index):
        """Converts results object into an hdf5 object.

        Parameters
        ----------
        handle : h5py.File or h5py.Group
            An HDF5 file or group type to store this in.
        index : int
            What step is this?

        """
        if "/number" not in handle:
            comm.barrier()
            self._write_hdf5_metadata(handle)

        comm.barrier()

        # Grab handles
        number_dset = handle["/number"]
        rxn_dset = handle["/reaction rates"]
        eigenvalues_dset = handle["/eigenvalues"]
        time_dset = handle["/time"]

        # Get number of results stored
        number_shape = list(number_dset.shape)
        number_results = number_shape[0]

        new_shape = index + 1

        if number_results < new_shape:
            # Extend first dimension by 1
            number_shape[0] = new_shape
            number_dset.resize(number_shape)

            rxn_shape = list(rxn_dset.shape)
            rxn_shape[0] = new_shape
            rxn_dset.resize(rxn_shape)

            eigenvalues_shape = list(eigenvalues_dset.shape)
            eigenvalues_shape[0] = new_shape
            eigenvalues_dset.resize(eigenvalues_shape)

            time_shape = list(time_dset.shape)
            time_shape[0] = new_shape
            time_dset.resize(time_shape)

        # If nothing to write, just return
        if len(self.mat_to_ind) == 0:
            return

        # Add data
        # Note, for the last step, self.n_stages = 1, even if n_stages != 1.
        n_stages = self.n_stages
        inds = [self.mat_to_hdf5_ind[mat] for mat in self.mat_to_ind]
        low = min(inds)
        high = max(inds)
        for i in range(n_stages):
            number_dset[index, i, low:high+1, :] = self.data[i, :, :]
            rxn_dset[index, i, low:high+1, :, :] = self.rates[i][:, :, :]
            if comm.rank == 0:
                eigenvalues_dset[index, i] = self.k[i]
        if comm.rank == 0:
            time_dset[index, :] = self.time

    @classmethod
    def from_hdf5(cls, handle, step):
        """Loads results object from HDF5.

        Parameters
        ----------
        handle : h5py.File or h5py.Group
            An HDF5 file or group type to load from.
        step : int
            What step is this?

        """
        results = cls()

        # Grab handles
        number_dset = handle["/number"]
        eigenvalues_dset = handle["/eigenvalues"]
        time_dset = handle["/time"]

        results.data = number_dset[step, :, :, :]
        results.k = eigenvalues_dset[step, :]
        results.time = time_dset[step, :]

        # Reconstruct dictionaries
        results.volume = OrderedDict()
        results.mat_to_ind = OrderedDict()
        results.nuc_to_ind = OrderedDict()
        rxn_nuc_to_ind = OrderedDict()
        rxn_to_ind = OrderedDict()

        for mat, mat_handle in handle["/materials"].items():
            vol = mat_handle.attrs["volume"]
            ind = mat_handle.attrs["index"]

            results.volume[mat] = vol
            results.mat_to_ind[mat] = ind

        for nuc, nuc_handle in handle["/nuclides"].items():
            ind_atom = nuc_handle.attrs["atom number index"]
            results.nuc_to_ind[nuc] = ind_atom

            if "reaction rate index" in nuc_handle.attrs:
                rxn_nuc_to_ind[nuc] = nuc_handle.attrs["reaction rate index"]

        for rxn, rxn_handle in handle["/reactions"].items():
            rxn_to_ind[rxn] = rxn_handle.attrs["index"]

        results.rates = []
        # Reconstruct reactions
        for i in range(results.n_stages):
            rate = ReactionRates(results.mat_to_ind, rxn_nuc_to_ind, rxn_to_ind)

            rate[:] = handle["/reaction rates"][step, i, :, :, :]
            results.rates.append(rate)

        return results

    @staticmethod
    def save(op, x, op_results, t, step_ind):
        """Creates and writes depletion results to disk

        Parameters
        ----------
        op : openmc.deplete.TransportOperator
            The operator used to generate these results.
        x : list of list of numpy.array
            The prior x vectors.  Indexed [i][cell] using the above equation.
        op_results : list of openmc.deplete.OperatorResult
            Results of applying transport operator
        t : list of float
            Time indices.
        step_ind : int
            Step index.

        """
        # Get indexing terms
        vol_dict, nuc_list, burn_list, full_burn_list = op.get_results_info()

        # Create results
        stages = len(x)
        results = Results()
        results.allocate(vol_dict, nuc_list, burn_list, full_burn_list, stages)

        n_mat = len(burn_list)

        for i in range(stages):
            for mat_i in range(n_mat):
                results[i, mat_i, :] = x[i][mat_i][:]

        results.k = [r.k for r in op_results]
        results.rates = [r.rates for r in op_results]
        results.time = t

        results.export_to_hdf5("depletion_results.h5", step_ind)
