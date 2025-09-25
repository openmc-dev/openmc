"""The stepresult module.

Contains capabilities for generating and saving results of a single depletion
timestep.
"""

import copy
import warnings
from pathlib import Path

import h5py
import numpy as np

import openmc
from openmc.mpi import comm, MPI
from openmc.checkvalue import PathLike
from .reaction_rates import ReactionRates

VERSION_RESULTS = (1, 1)


__all__ = ["StepResult"]


class StepResult:
    """Result of a single depletion timestep

    .. versionchanged:: 0.13.1
        Name changed from ``Results`` to ``StepResult``

    Attributes
    ----------
    k : list of (float, float)
        Eigenvalue and uncertainty for each substep.
    time : list of float
        Time at beginning, end of step, in seconds.
    source_rate : float
        Source rate during timestep in [W] or [neutron/sec]
    n_mat : int
        Number of mats.
    n_nuc : int
        Number of nuclides.
    rates : list of ReactionRates
        The reaction rates for each substep.
    volume : dict of str to float
        Dictionary mapping mat id to volume.
    index_mat : dict of str to int
        A dictionary mapping mat ID as string to index.
    index_nuc : dict of str to int
        A dictionary mapping nuclide name as string to index.
    mat_to_hdf5_ind : dict of str to int
        A dictionary mapping mat ID as string to global index.
    n_hdf5_mats : int
        Number of materials in entire geometry.
    n_stages : int
        Number of stages in simulation.
    data : numpy.ndarray
        Atom quantity, stored by stage, mat, then by nuclide.
    proc_time : int
        Average time spent depleting a material across all
        materials and processes

    """
    def __init__(self):
        self.k = None
        self.time = None
        self.source_rate = None
        self.rates = None
        self.volume = None
        self.proc_time = None

        self.index_mat = None
        self.index_nuc = None
        self.mat_to_hdf5_ind = None

        self.data = None

    def __repr__(self):
        t = self.time[0]
        dt = self.time[1] - self.time[0]
        return f"<StepResult: t={t}, dt={dt}, source={self.source_rate}>"

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
        if isinstance(mat, openmc.Material):
            mat = str(mat.id)
        if isinstance(mat, str):
            mat = self.index_mat[mat]
        if isinstance(nuc, str):
            nuc = self.index_nuc[nuc]

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
            mat = self.index_mat[mat]
        if isinstance(nuc, str):
            nuc = self.index_nuc[nuc]

        self.data[stage, mat, nuc] = val

    @property
    def n_mat(self):
        return len(self.index_mat)

    @property
    def n_nuc(self):
        return len(self.index_nuc)

    @property
    def n_hdf5_mats(self):
        return len(self.mat_to_hdf5_ind)

    @property
    def n_stages(self):
        return self.data.shape[0]

    def allocate(self, volume, nuc_list, burn_list, full_burn_list, stages):
        """Allocate memory for depletion step data

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
        self.index_nuc = {nuc: i for i, nuc in enumerate(nuc_list)}
        self.index_mat = {mat: i for i, mat in enumerate(burn_list)}
        self.mat_to_hdf5_ind = {mat: i for i, mat in enumerate(full_burn_list)}

        # Create storage array
        self.data = np.zeros((stages, self.n_mat, self.n_nuc))

    def distribute(self, local_materials, ranges):
        """Create a new object containing data for distributed materials

        Parameters
        ----------
        local_materials : iterable of str
            Materials for this process
        ranges : iterable of int
            Slice-like object indicating indicies of ``local_materials``
            in the material dimension of :attr:`data` and each element
            in :attr:`rates`

        Returns
        -------
        StepResult
            New results object
        """
        new = StepResult()
        new.volume = {lm: self.volume[lm] for lm in local_materials}
        new.index_mat = {mat: idx for (idx, mat) in enumerate(local_materials)}

        # Direct transfer
        direct_attrs = ("time", "k", "source_rate", "index_nuc",
                        "mat_to_hdf5_ind", "proc_time")
        for attr in direct_attrs:
            setattr(new, attr, getattr(self, attr))
        # Get applicable slice of data
        new.data = self.data[:, ranges]
        new.rates = [r[ranges] for r in self.rates]
        return new

    def get_material(self, mat_id):
        """Return material object for given depleted composition

        .. versionadded:: 0.13.2

        Parameters
        ----------
        mat_id : str
            Material ID as a string

        Returns
        -------
        openmc.Material
            Equivalent material

        Raises
        ------
        KeyError
            If specified material ID is not found in the StepResult

        """
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', openmc.IDWarning)
            material = openmc.Material(material_id=int(mat_id))
        try:
            vol = self.volume[mat_id]
        except KeyError as e:
            raise KeyError(
                f'mat_id {mat_id} not found in StepResult. Available mat_id '
                f'values are {list(self.volume.keys())}'
            ) from e
        for nuc, _ in sorted(self.index_nuc.items(), key=lambda x: x[1]):
            atoms = self[0, mat_id, nuc]
            if atoms <= 0.0:
                continue
            atom_per_bcm = atoms / vol * 1e-24
            material.add_nuclide(nuc, atom_per_bcm)
        material.volume = vol
        return material

    def export_to_hdf5(self, filename, step, write_reaction_rates=True):
        """Export results to an HDF5 file

        Parameters
        ----------
        filename : str
            The filename to write to
        step : int
            What step is this?
        write_reaction_rates : bool, optional
            Whether to include reaction rate datasets in the results file.

        """
        # Write new file if first time step, else add to existing file
        kwargs = {'mode': "w" if step == 0 else "a"}

        if h5py.get_config().mpi and comm.size > 1:
            # Write results in parallel
            kwargs['driver'] = 'mpio'
            kwargs['comm'] = comm
            with h5py.File(filename, **kwargs) as handle:
                self._to_hdf5(handle, step, parallel=True,
                              write_reaction_rates=write_reaction_rates)
        else:
            # Gather results at root process
            all_results = comm.gather(self)

            # Only root process writes results
            if comm.rank == 0:
                with h5py.File(filename, **kwargs) as handle:
                    for res in all_results:
                        res._to_hdf5(handle, step, parallel=False,
                                     write_reaction_rates=write_reaction_rates)

    def _write_hdf5_metadata(self, handle, write_reaction_rates):
        """Writes result metadata in HDF5 file

        Parameters
        ----------
        handle : h5py.File or h5py.Group
            An hdf5 file or group type to store this in.
        write_reaction_rates : bool
            Whether reaction rate datasets are being written.

        """
        # Create and save the 5 dictionaries:
        # quantities
        #   self.index_mat -> self.volume (TODO: support for changing volumes)
        #   self.index_nuc
        # reactions
        #   self.rates[0].index_nuc (can be different from above, above is superset)
        #   self.rates[0].index_rx
        # these are shared by every step of the simulation, and should be deduplicated.

        # Store concentration mat and nuclide dictionaries (along with volumes)

        handle.attrs['version'] = np.array(VERSION_RESULTS)
        handle.attrs['filetype'] = np.bytes_('depletion results')

        mat_list = sorted(self.mat_to_hdf5_ind, key=int)
        nuc_list = sorted(self.index_nuc)

        include_rates = (
            write_reaction_rates
            and self.rates
            and len(self.rates) > 0
            and bool(self.rates[0].index_nuc)
            and bool(self.rates[0].index_rx)
        )
        rxn_list = sorted(self.rates[0].index_rx) if include_rates else []

        n_mats = self.n_hdf5_mats
        n_nuc_number = len(nuc_list)
        n_nuc_rxn = len(self.rates[0].index_nuc) if include_rates else 0
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
            nuc_single_group.attrs["atom number index"] = self.index_nuc[nuc]
            if include_rates and nuc in self.rates[0].index_nuc:
                nuc_single_group.attrs["reaction rate index"] = (
                    self.rates[0].index_nuc[nuc])

        if include_rates:
            rxn_group = handle.create_group("reactions")

            for rxn in rxn_list:
                rxn_single_group = rxn_group.create_group(rxn)
                rxn_single_group.attrs["index"] = (
                    self.rates[0].index_rx[rxn])

        # Construct array storage

        handle.create_dataset("number", (1, n_stages, n_mats, n_nuc_number),
                              maxshape=(None, n_stages, n_mats, n_nuc_number),
                              chunks=True,
                              dtype='float64')

        if include_rates and n_nuc_rxn > 0 and n_rxn > 0:
            handle.create_dataset(
                "reaction rates", (1, n_stages, n_mats, n_nuc_rxn, n_rxn),
                maxshape=(None, n_stages, n_mats, n_nuc_rxn, n_rxn),
                chunks=True, dtype='float64')

        handle.create_dataset("eigenvalues", (1, n_stages, 2),
                              maxshape=(None, n_stages, 2), dtype='float64')

        handle.create_dataset("time", (1, 2), maxshape=(None, 2), dtype='float64')

        handle.create_dataset("source_rate", (1, n_stages), maxshape=(None, n_stages),
                              dtype='float64')

        handle.create_dataset(
            "depletion time", (1,), maxshape=(None,),
            dtype="float64")

    def _to_hdf5(self, handle, index, parallel=False, write_reaction_rates=True):
        """Converts results object into an hdf5 object.

        Parameters
        ----------
        handle : h5py.File or h5py.Group
            An HDF5 file or group type to store this in.
        index : int
            What step is this?
        parallel : bool
            Being called with parallel HDF5?
        write_reaction_rates : bool, optional
            Whether reaction rate datasets are being written.

        """
        if "/number" not in handle:
            if parallel:
                comm.barrier()
            self._write_hdf5_metadata(handle, write_reaction_rates)

        if parallel:
            comm.barrier()

        # Grab handles
        number_dset = handle["/number"]
        has_reactions = ("reaction rates" in handle)
        if has_reactions:
            rxn_dset = handle["/reaction rates"]
        eigenvalues_dset = handle["/eigenvalues"]
        time_dset = handle["/time"]
        source_rate_dset = handle["/source_rate"]
        proc_time_dset = handle["/depletion time"]

        # Get number of results stored
        number_shape = list(number_dset.shape)
        number_results = number_shape[0]

        new_shape = index + 1

        if number_results < new_shape:
            # Extend first dimension by 1
            number_shape[0] = new_shape
            number_dset.resize(number_shape)

            if has_reactions:
                rxn_shape = list(rxn_dset.shape)
                rxn_shape[0] = new_shape
                rxn_dset.resize(rxn_shape)

            eigenvalues_shape = list(eigenvalues_dset.shape)
            eigenvalues_shape[0] = new_shape
            eigenvalues_dset.resize(eigenvalues_shape)

            time_shape = list(time_dset.shape)
            time_shape[0] = new_shape
            time_dset.resize(time_shape)

            source_rate_shape = list(source_rate_dset.shape)
            source_rate_shape[0] = new_shape
            source_rate_dset.resize(source_rate_shape)

            proc_shape = list(proc_time_dset.shape)
            proc_shape[0] = new_shape
            proc_time_dset.resize(proc_shape)

        # If nothing to write, just return
        if len(self.index_mat) == 0:
            return

        # Add data
        # Note, for the last step, self.n_stages = 1, even if n_stages != 1.
        n_stages = self.n_stages
        inds = [self.mat_to_hdf5_ind[mat] for mat in self.index_mat]
        low = min(inds)
        high = max(inds)
        for i in range(n_stages):
            number_dset[index, i, low:high+1] = self.data[i]
            if has_reactions:
                rxn_dset[index, i, low:high+1] = self.rates[i]
            if comm.rank == 0:
                eigenvalues_dset[index, i] = self.k[i]
        if comm.rank == 0:
            time_dset[index] = self.time
            source_rate_dset[index] = self.source_rate
            if self.proc_time is not None:
                proc_time_dset[index] = (
                    self.proc_time / (comm.size * self.n_hdf5_mats)
                )

    @classmethod
    def from_hdf5(cls, handle, step):
        """Loads results object from HDF5.

        Parameters
        ----------
        handle : h5py.File or h5py.Group
            An HDF5 file or group type to load from.
        step : int
            Index for depletion step
        """
        results = cls()

        # Grab handles
        number_dset = handle["/number"]
        eigenvalues_dset = handle["/eigenvalues"]
        time_dset = handle["/time"]
        if "source_rate" in handle:
            source_rate_dset = handle["/source_rate"]
        else:
            # Older versions used "power" instead of "source_rate"
            source_rate_dset = handle["/power"]

        results.data = number_dset[step, :, :, :]
        results.k = eigenvalues_dset[step, :]
        results.time = time_dset[step, :]
        results.source_rate = source_rate_dset[step, 0]

        if "depletion time" in handle:
            proc_time_dset = handle["/depletion time"]
            if step < proc_time_dset.shape[0]:
                results.proc_time = proc_time_dset[step]

        if results.proc_time is None:
            results.proc_time = np.array([np.nan])

        # Reconstruct dictionaries
        results.volume = {}
        results.index_mat = {}
        results.index_nuc = {}
        rxn_nuc_to_ind = {}
        rxn_to_ind = {}

        for mat, mat_handle in handle["/materials"].items():
            vol = mat_handle.attrs["volume"]
            ind = mat_handle.attrs["index"]

            results.volume[mat] = vol
            results.index_mat[mat] = ind

        for nuc, nuc_handle in handle["/nuclides"].items():
            ind_atom = nuc_handle.attrs["atom number index"]
            results.index_nuc[nuc] = ind_atom

            if "reaction rate index" in nuc_handle.attrs:
                rxn_nuc_to_ind[nuc] = nuc_handle.attrs["reaction rate index"]

        if "reactions" in handle:
            for rxn, rxn_handle in handle["/reactions"].items():
                rxn_to_ind[rxn] = rxn_handle.attrs["index"]

        results.rates = []
        # Reconstruct reactions
        for i in range(results.n_stages):
            rate = ReactionRates(results.index_mat, rxn_nuc_to_ind, rxn_to_ind, True)

            if "reaction rates" in handle:
                rate[:] = handle["/reaction rates"][step, i, :, :, :]
            results.rates.append(rate)

        return results

    @staticmethod
    def save(
        op,
        x,
        op_results,
        t,
        source_rate,
        step_ind,
        proc_time=None,
        write_reaction_rates: bool = False,
        path: PathLike = "depletion_results.h5"
    ):
        """Creates and writes depletion results to disk

        Parameters
        ----------
        op : openmc.deplete.abc.TransportOperator
            The operator used to generate these results.
        x : list of list of numpy.array
            The prior x vectors.  Indexed [i][cell] using the above equation.
        op_results : list of openmc.deplete.OperatorResult
            Results of applying transport operator
        t : list of float
            Time indices.
        source_rate : float
            Source rate during time step in [W] or [neutron/sec]
        step_ind : int
            Step index.
        proc_time : float or None
            Total process time spent depleting materials. This may
            be process-dependent and will be reduced across MPI
            processes.
        write_reaction_rates : bool, optional
            Whether reaction rates should be written to the results file.
        path : PathLike
            Path to file to write. Defaults to 'depletion_results.h5'.

            .. versionadded:: 0.14.0
        """
        # Get indexing terms
        vol_dict, nuc_list, burn_list, full_burn_list = op.get_results_info()

        stages = len(x)

        # Create results
        results = StepResult()
        results.allocate(vol_dict, nuc_list, burn_list, full_burn_list, stages)

        n_mat = len(burn_list)

        for i in range(stages):
            for mat_i in range(n_mat):
                results[i, mat_i, :] = x[i][mat_i]

        ks = []
        for r in op_results:
            if isinstance(r.k, type(None)):
                ks += [(None, None)]
            else:
                ks += [(r.k.nominal_value, r.k.std_dev)]
        results.k = ks
        results.rates = [r.rates for r in op_results]
        results.time = t
        results.source_rate = source_rate
        results.proc_time = proc_time
        if results.proc_time is not None:
            results.proc_time = comm.reduce(proc_time, op=MPI.SUM)

        if not Path(path).is_file():
            Path(path).parent.mkdir(parents=True, exist_ok=True)
        results.export_to_hdf5(path, step_ind, write_reaction_rates)

    def transfer_volumes(self, model):
        """Transfers volumes from depletion results to geometry

        Parameters
        ----------
        model : OpenMC model to be used in a depletion restart
            calculation

        """

        if not model.materials:
            materials = openmc.Materials(
                model.geometry.get_all_materials().values()
            )
        else:
            materials = model.materials

        for material in materials:
            if material.depletable:
                material.volume = self.volume[str(material.id)]
