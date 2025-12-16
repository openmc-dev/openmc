from datetime import datetime
from collections import namedtuple
import glob
import re
import os
import warnings

import h5py
import numpy as np
from pathlib import Path
from uncertainties import ufloat
from uncertainties.unumpy import uarray

import openmc
import openmc.checkvalue as cv

_VERSION_STATEPOINT = 18


KineticsParameters = namedtuple("KineticsParameters", ["generation_time", "beta_effective"])


class StatePoint:
    """State information on a simulation at a certain point in time (at the end
    of a given batch). Statepoints can be used to analyze tally results as well
    as restart a simulation.

    Parameters
    ----------
    filepath : str or Path
        Path to file to load
    autolink : bool, optional
        Whether to automatically link in metadata from a summary.h5 file and
        stochastic volume calculation results from volume_*.h5 files. Defaults
        to True.

    Attributes
    ----------
    cmfd_on : bool
        Indicate whether CMFD is active
    cmfd_balance : numpy.ndarray
        Residual neutron balance for each batch
    cmfd_dominance
        Dominance ratio for each batch
    cmfd_entropy : numpy.ndarray
        Shannon entropy of CMFD fission source for each batch
    cmfd_indices : numpy.ndarray
        Number of CMFD mesh cells and energy groups. The first three indices
        correspond to the x-, y-, and z- spatial directions and the fourth index
        is the number of energy groups.
    cmfd_srccmp : numpy.ndarray
        Root-mean-square difference between OpenMC and CMFD fission source for
        each batch
    cmfd_src : numpy.ndarray
        CMFD fission source distribution over all mesh cells and energy groups.
    current_batch : int
        Number of batches simulated
    date_and_time : datetime.datetime
        Date and time at which statepoint was written
    entropy : numpy.ndarray
        Shannon entropy of fission source at each batch
    filters : dict
        Dictionary whose keys are filter IDs and whose values are Filter
        objects
    generations_per_batch : int
        Number of fission generations per batch
    global_tallies : numpy.ndarray of compound datatype
        Global tallies for k-effective estimates and leakage. The compound
        datatype has fields 'name', 'sum', 'sum_sq', 'mean', and 'std_dev'.
    k_combined : uncertainties.UFloat
        Combined estimator for k-effective

        .. deprecated:: 0.13.1
    k_col_abs : float
        Cross-product of collision and absorption estimates of k-effective
    k_col_tra : float
        Cross-product of collision and tracklength estimates of k-effective
    k_abs_tra : float
        Cross-product of absorption and tracklength estimates of k-effective
    k_generation : numpy.ndarray
        Estimate of k-effective for each batch/generation
    keff : uncertainties.UFloat
        Combined estimator for k-effective

        .. versionadded:: 0.13.1
    meshes : dict
        Dictionary whose keys are mesh IDs and whose values are MeshBase objects
    n_batches : int
        Number of batches
    n_inactive : int
        Number of inactive batches
    n_particles : int
        Number of particles per generation
    n_realizations : int
        Number of tally realizations
    path : str
        Working directory for simulation
    photon_transport : bool
        Indicate whether photon transport is active
    run_mode : str
        Simulation run mode, e.g. 'eigenvalue'
    runtime : dict
        Dictionary whose keys are strings describing various runtime metrics
        and whose values are time values in seconds.
    seed : int
        Pseudorandom number generator seed
    stride : int
        Number of random numbers allocated for each particle history
    source : numpy.ndarray of compound datatype
        Array of source sites. The compound datatype has fields 'r', 'u',
        'E', 'wgt', 'delayed_group', 'surf_id', and 'particle', corresponding to
        the position, direction, energy, weight, delayed group, surface ID and
        particle type of the source site, respectively.
    source_present : bool
        Indicate whether source sites are present
    sparse : bool
        Whether or not the tallies uses SciPy's LIL sparse matrix format for
        compressed data storage
    tallies : dict
        Dictionary whose keys are tally IDs and whose values are Tally objects
    tallies_present : bool
        Indicate whether user-defined tallies are present
    tally_derivatives : dict
        Dictionary whose keys are tally derivative IDs and whose values are
        TallyDerivative objects
    version: tuple of Integral
        Version of OpenMC
    summary : None or openmc.Summary
        A summary object if the statepoint has been linked with a summary file

    """

    def __init__(self, filepath, autolink=True):
        filename = str(filepath)  # in case it's a Path
        self._f = h5py.File(filename, 'r')
        self._meshes = {}
        self._filters = {}
        self._tallies = {}
        self._derivs = {}

        # Check filetype and version
        cv.check_filetype_version(self._f, 'statepoint', _VERSION_STATEPOINT)

        # Set flags for what data has been read
        self._meshes_read = False
        self._filters_read = False
        self._tallies_read = False
        self._summary = None
        self._global_tallies = None
        self._sparse = False
        self._derivs_read = False

        # Automatically link in a summary file if one exists
        if autolink:
            path_summary = os.path.join(os.path.dirname(filename), 'summary.h5')
            if os.path.exists(path_summary):
                su = openmc.Summary(path_summary)
                self.link_with_summary(su)

            path_volume = os.path.join(os.path.dirname(filename), 'volume_*.h5')
            for path_i in glob.glob(path_volume):
                if re.search(r'volume_\d+\.h5', path_i):
                    vol = openmc.VolumeCalculation.from_hdf5(path_i)
                    self.add_volume_information(vol)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    @property
    def cmfd_on(self):
        return self._f.attrs['cmfd_on'] > 0

    @property
    def cmfd_balance(self):
        return self._f['cmfd/cmfd_balance'][()] if self.cmfd_on else None

    @property
    def cmfd_dominance(self):
        return self._f['cmfd/cmfd_dominance'][()] if self.cmfd_on else None

    @property
    def cmfd_entropy(self):
        return self._f['cmfd/cmfd_entropy'][()] if self.cmfd_on else None

    @property
    def cmfd_indices(self):
        return self._f['cmfd/indices'][()] if self.cmfd_on else None

    @property
    def cmfd_src(self):
        if self.cmfd_on:
            data = self._f['cmfd/cmfd_src'][()]
            return np.reshape(data, tuple(self.cmfd_indices), order='F')
        else:
            return None

    @property
    def cmfd_srccmp(self):
        return self._f['cmfd/cmfd_srccmp'][()] if self.cmfd_on else None

    @property
    def current_batch(self):
        return self._f['current_batch'][()]

    @property
    def date_and_time(self):
        s = self._f.attrs['date_and_time'].decode()
        return datetime.strptime(s, '%Y-%m-%d %H:%M:%S')

    @property
    def entropy(self):
        if self.run_mode == 'eigenvalue':
            return self._f['entropy'][()]
        else:
            return None

    @property
    def filters(self):
        if not self._filters_read:
            filters_group = self._f['tallies/filters']

            # Iterate over all Filters
            for group in filters_group.values():
                new_filter = openmc.Filter.from_hdf5(group, meshes=self.meshes)
                self._filters[new_filter.id] = new_filter

            self._filters_read = True

        return self._filters

    @property
    def generations_per_batch(self):
        if self.run_mode == 'eigenvalue':
            return self._f['generations_per_batch'][()]
        else:
            return None

    @property
    def global_tallies(self):
        if self._global_tallies is None:
            data = self._f['global_tallies'][()]
            gt = np.zeros(data.shape[0], dtype=[
                ('name', 'S14'), ('sum', 'f8'), ('sum_sq', 'f8'),
                ('mean', 'f8'), ('std_dev', 'f8')])
            gt['name'] = ['k-collision', 'k-absorption', 'k-tracklength',
                          'leakage']
            gt['sum'] = data[:,1]
            gt['sum_sq'] = data[:,2]

            # Calculate mean and sample standard deviation of mean
            n = self.n_realizations
            gt['mean'] = gt['sum']/n
            gt['std_dev'] = np.sqrt((gt['sum_sq']/n - gt['mean']**2)/(n - 1))

            self._global_tallies = gt

        return self._global_tallies

    @property
    def k_cmfd(self):
        if self.cmfd_on:
            return self._f['cmfd/k_cmfd'][()]
        else:
            return None

    @property
    def k_generation(self):
        if self.run_mode == 'eigenvalue':
            return self._f['k_generation'][()]
        else:
            return None

    @property
    def keff(self):
        if self.run_mode == 'eigenvalue':
            return ufloat(*self._f['k_combined'][()])
        else:
            return None

    @property
    def k_combined(self):
        warnings.warn(
            "The 'k_combined' property has been renamed to 'keff' and will be "
            "removed in a future version of OpenMC.", FutureWarning
        )
        return self.keff

    @property
    def k_col_abs(self):
        if self.run_mode == 'eigenvalue':
            return self._f['k_col_abs'][()]
        else:
            return None

    @property
    def k_col_tra(self):
        if self.run_mode == 'eigenvalue':
            return self._f['k_col_tra'][()]
        else:
            return None

    @property
    def k_abs_tra(self):
        if self.run_mode == 'eigenvalue':
            return self._f['k_abs_tra'][()]
        else:
            return None

    @property
    def meshes(self):
        if not self._meshes_read:
            mesh_group = self._f['tallies/meshes']

            # Iterate over all meshes
            for group in mesh_group.values():
                mesh = openmc.MeshBase.from_hdf5(group)
                self._meshes[mesh.id] = mesh

            self._meshes_read = True

        return self._meshes

    @property
    def n_batches(self):
        return self._f['n_batches'][()]

    @property
    def n_inactive(self):
        if self.run_mode == 'eigenvalue':
            return self._f['n_inactive'][()]
        else:
            return None

    @property
    def n_particles(self):
        return self._f['n_particles'][()]

    @property
    def n_realizations(self):
        return self._f['n_realizations'][()]

    @property
    def path(self):
        return self._f.attrs['path'].decode()

    @property
    def photon_transport(self):
        return self._f.attrs['photon_transport'] > 0

    @property
    def run_mode(self):
        return self._f['run_mode'][()].decode()

    @property
    def runtime(self):
        return {name: dataset[()]
                for name, dataset in self._f['runtime'].items()}

    @property
    def seed(self):
        return self._f['seed'][()]

    @property
    def stride(self):
        return self._f['stride'][()]

    @property
    def source(self):
        return self._f['source_bank'][()] if self.source_present else None

    @property
    def source_present(self):
        return self._f.attrs['source_present'] > 0

    @property
    def sparse(self):
        return self._sparse

    @sparse.setter
    def sparse(self, sparse):
        """Convert tally data from NumPy arrays to SciPy list of lists (LIL)
        sparse matrices, and vice versa.

        This property may be used to reduce the amount of data in memory during
        tally data processing. The tally data will be stored as SciPy LIL
        matrices internally within each Tally object. All tally data access
        properties and methods will return data as a dense NumPy array.

        """

        cv.check_type('sparse', sparse, bool)
        self._sparse = sparse

        # Update tally sparsities
        if self._tallies_read:
            for tally_id in self.tallies:
                self.tallies[tally_id].sparse = self.sparse

    @property
    def tallies(self):
        if self.tallies_present and not self._tallies_read:
            # Read the number of tallies
            tallies_group = self._f['tallies']
            n_tallies = tallies_group.attrs['n_tallies']

            # Read a list of the IDs for each Tally
            if n_tallies > 0:
                # Tally user-defined IDs
                tally_ids = tallies_group.attrs['ids']
            else:
                tally_ids = []

            # Ignore warnings about duplicate IDs
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', openmc.IDWarning)

                # Iterate over all tallies
                for tally_id in tally_ids:
                    group = tallies_group[f'tally {tally_id}']

                    # Check if tally is internal and therefore has no data
                    if group.attrs.get("internal"):
                        continue

                    # Create Tally object and assign basic properties
                    tally = openmc.Tally(tally_id)
                    tally._sp_filename = Path(self._f.filename)
                    tally.name = group['name'][()].decode() if 'name' in group else ''

                    # Check if tally has multiply_density attribute
                    if "multiply_density" in group.attrs:
                        tally.multiply_density = group.attrs["multiply_density"].item() > 0

                    # Check if tally has higher_moments attribute
                    if 'higher_moments' in group.attrs:
                        tally.higher_moments = bool(group.attrs['higher_moments'][()])

                    # Read the number of realizations
                    n_realizations = group['n_realizations'][()]

                    tally.estimator = group['estimator'][()].decode()
                    tally.num_realizations = n_realizations

                    # Read derivative information.
                    if 'derivative' in group:
                        deriv_id = group['derivative'][()]
                        tally.derivative = self.tally_derivatives[deriv_id]

                    # Read all filters
                    n_filters = group['n_filters'][()]
                    if n_filters > 0:
                        filter_ids = group['filters'][()]
                        filters_group = self._f['tallies/filters']
                        for filter_id in filter_ids:
                            filter_group = filters_group[f'filter {filter_id}']
                            new_filter = openmc.Filter.from_hdf5(
                                filter_group, meshes=self.meshes)
                            tally.filters.append(new_filter)

                    # Read nuclide bins
                    nuclide_names = group['nuclides'][()]

                    # Add all nuclides to the Tally
                    tally.nuclides = [name.decode().strip() for name in nuclide_names]

                    # Add the scores to the Tally
                    scores = group['score_bins'][()]
                    tally.scores = [score.decode() for score in scores]

                    # Add Tally to the global dictionary of all Tallies
                    tally.sparse = self.sparse
                    self._tallies[tally_id] = tally

            self._tallies_read = True

        return self._tallies

    @property
    def tallies_present(self):
        return self._f.attrs['tallies_present'] > 0

    @property
    def tally_derivatives(self):
        if not self._derivs_read:
            # Populate the dictionary if any derivatives are present.
            if 'derivatives' in self._f['tallies']:
                # Read the derivative ids.
                base = 'tallies/derivatives'
                deriv_ids = [int(k.split(' ')[1]) for k in self._f[base]]

                # Create each derivative object and add it to the dictionary.
                for d_id in deriv_ids:
                    group = self._f[f'tallies/derivatives/derivative {d_id}']
                    deriv = openmc.TallyDerivative(derivative_id=d_id)
                    deriv.variable = group['independent variable'][()].decode()
                    if deriv.variable == 'density':
                        deriv.material = group['material'][()]
                    elif deriv.variable == 'nuclide_density':
                        deriv.material = group['material'][()]
                        deriv.nuclide = group['nuclide'][()].decode()
                    elif deriv.variable == 'temperature':
                        deriv.material = group['material'][()]
                    self._derivs[d_id] = deriv

            self._derivs_read = True

        return self._derivs

    @property
    def version(self):
        return tuple(self._f.attrs['openmc_version'])

    @property
    def summary(self):
        return self._summary

    def close(self):
        """Close the statepoint HDF5 file and the corresponding
        summary HDF5 file if present.
        """
        self._f.close()
        if self._summary is not None:
            self._summary._f.close()

    def add_volume_information(self, volume_calc):
        """Add volume information to the geometry within the file

        Parameters
        ----------
        volume_calc : openmc.VolumeCalculation
            Results from a stochastic volume calculation

        """
        if self.summary is not None:
            self.summary.add_volume_information(volume_calc)

    def get_tally(self, scores=[], filters=[], nuclides=[],
                  name=None, id=None, estimator=None, exact_filters=False,
                  exact_nuclides=False, exact_scores=False,
                  multiply_density=None, derivative=None, filter_type=None):
        """Finds and returns a Tally object with certain properties.

        This routine searches the list of Tallies and returns the first Tally
        found which satisfies all of the input parameters.

        NOTE: If any of the "exact" parameters are False (default), the input
        parameters do not need to match the complete Tally specification and may
        only represent a subset of the Tally's properties. If an "exact"
        parameter is True then number of scores, filters, or nuclides in the
        parameters must precisely match those of any matching Tally.

        Parameters
        ----------
        scores : list, optional
            A list of one or more score strings (default is []).
        filters : list, optional
            A list of Filter objects (default is []).
        nuclides : list, optional
            A list of Nuclide objects (default is []).
        name : str, optional
            The name specified for the Tally (default is None).
        id : Integral, optional
            The id specified for the Tally (default is None).
        estimator: str, optional
            The type of estimator ('tracklength', 'analog'; default is None).
        exact_filters : bool
            If True, the number of filters in the parameters must be identical
            to those in the matching Tally. If False (default), the filters in
            the parameters may be a subset of those in the matching Tally.
        exact_nuclides : bool
            If True, the number of nuclides in the parameters must be identical
            to those in the matching Tally. If False (default), the nuclides in
            the parameters may be a subset of those in the matching Tally.
        exact_scores : bool
            If True, the number of scores in the parameters must be identical to
            those in the matching Tally. If False (default), the scores in the
            parameters may be a subset of those in the matching Tally. Default
            is None (no check).
        multiply_density : bool, optional
            If not None, the Tally must have the multiply_density attribute set
            to the same value as this parameter.
        derivative : openmc.TallyDerivative, optional
            TallyDerivative object to match.
        filter_type : type, optional
            If not None, the Tally must have at least one Filter that is an
            instance of this type. For example `openmc.MeshFilter`.

        Returns
        -------
        tally : openmc.Tally
            A tally matching the specified criteria

        Raises
        ------
        LookupError
            If a Tally meeting all of the input parameters cannot be found in
            the statepoint.

        """

        tally = None

        # Iterate over all tallies to find the appropriate one
        for test_tally in self.tallies.values():

            # Determine if Tally has queried name
            if name and name != test_tally.name:
                continue

            # Determine if Tally has queried id
            if id and id != test_tally.id:
                continue

            # Determine if Tally has queried estimator, only move on to next tally
            # if the estimator is both specified and the tally estimtor does not
            # match
            if estimator is not None and estimator != test_tally.estimator:
                continue

            # The number of filters, nuclides and scores must exactly match
            if exact_scores and len(scores) != test_tally.num_scores:
                continue
            if exact_nuclides and nuclides and len(nuclides) != test_tally.num_nuclides:
                continue
            if exact_nuclides and not nuclides and test_tally.nuclides != ['total']:
                continue
            if exact_filters and len(filters) != test_tally.num_filters:
                continue
            if derivative is not None and derivative != test_tally.derivative:
                continue
            if multiply_density is not None and multiply_density != test_tally.multiply_density:
                continue

            # Determine if Tally has the queried score(s)
            if scores:
                if not all(score in test_tally.scores for score in scores):
                    continue

            # Determine if Tally has the queried Filter(s)
            if filters:
                contains_filters = True

                # Iterate over the Filters requested by the user
                for outer_filter in filters:
                    contains_filters = False

                    # Test if requested filter is a subset of any of the test
                    # tally's filters and if so continue to next filter
                    for inner_filter in test_tally.filters:
                        if inner_filter.is_subset(outer_filter):
                            contains_filters = True
                            break

                    if not contains_filters:
                        break

                if not contains_filters:
                    continue

            if filter_type is not None:
                if not any(isinstance(f, filter_type) for f in test_tally.filters):
                    continue

            # Determine if Tally has the queried Nuclide(s)
            if nuclides:
                if not all(nuclide in test_tally.nuclides for nuclide in nuclides):
                    continue

            # If the current Tally met user's request, break loop and return it
            tally = test_tally
            break

        # If we did not find the Tally, return an error message
        if tally is None:
            raise LookupError('Unable to get Tally')

        return tally

    def get_mesh(
        self,
        id: int | None = None,
        name: str | None = None,
        mesh_type: openmc.MeshBase | None = None,
    ):
        """Return a Mesh object that matches all of the input parameters.

        Parameters
        ----------
        name : str, optional
            The name specified for the mesh (default is None).
        id : int, optional
            The id specified for the mesh (default is None).
        mesh_type : openmc.Mesh, optional
            The type of MeshBase, for example openmc.RegularMesh (default is None).

        Returns
        -------
        mesh : openmc.Mesh
            A mesh matching the specified criteria

        Raises
        ------
        LookupError
            If a mesh meeting all of the input parameters cannot be found in the
            statepoint.

        """

        mesh = None

        for test_mesh in self.meshes.values():
            if (
                (id and id != test_mesh.id)
                or (name and name != test_mesh.name)
                or (mesh_type and not isinstance(test_mesh, mesh_type))
            ):
                continue

            # If the current mesh met user's request, break loop and return it
            mesh = test_mesh
            break

        # If we did not find the mesh, return an error message
        if mesh is None:
            raise LookupError("Unable to get Mesh")

        return mesh

    def link_with_summary(self, summary):
        """Links Tallies and Filters with Summary model information.

        This routine retrieves model information (materials, geometry) from a
        Summary object populated with an HDF5 'summary.h5' file and inserts it
        into the Tally objects. This can be helpful when viewing and
        manipulating large scale Tally data. Note that it is necessary to link
        against a summary to populate the Tallies with any user-specified "name"
        XML tags.

        Parameters
        ----------
        summary : openmc.Summary
             A Summary object.

        Raises
        ------
        ValueError
            An error when the argument passed to the 'summary' parameter is not
            an openmc.Summary object.

        """

        if self.summary is not None:
            warnings.warn('A Summary object has already been linked.',
                          RuntimeWarning)
            return

        if not isinstance(summary, openmc.Summary):
            msg = f'Unable to link statepoint with "{summary}" which is not a' \
                  'Summary object'
            raise ValueError(msg)

        cells = summary.geometry.get_all_cells()

        for tally in self.tallies.values():
            tally.with_summary = True

            for tally_filter in tally.filters:
                if isinstance(tally_filter, (openmc.DistribcellFilter)):
                    cell_id = tally_filter.bins[0]
                    cell = cells[cell_id]
                    if not cell._paths:
                        summary.geometry.determine_paths()
                    tally_filter._paths = cell.paths

        self._summary = summary

    def get_kinetics_parameters(self) -> KineticsParameters:
        """Get kinetics parameters from IFP tallies.

        This method searches the tallies in the statepoint for the tallies
        required to compute kinetics parameters using the Iterated Fission
        Probability (IFP) method.

        Returns
        -------
        KineticsParameters
            A named tuple containing the generation time and effective delayed
            neutron fraction. If the necessary tallies for one or both
            parameters are not found, that parameter is returned as None.

        """

        denom_tally = None
        gen_time_tally = None
        beta_tally = None
        for tally in self.tallies.values():
            if 'ifp-denominator' in tally.scores:
                denom_tally = self.get_tally(scores=['ifp-denominator'])
            if 'ifp-time-numerator' in tally.scores:
                gen_time_tally = self.get_tally(scores=['ifp-time-numerator'])
            if 'ifp-beta-numerator' in tally.scores:
                beta_tally = self.get_tally(scores=['ifp-beta-numerator'])

        if denom_tally is None:
            return KineticsParameters(None, None)

        def get_ufloat(tally, score):
            return uarray(tally.get_values(scores=[score]),
                          tally.get_values(scores=[score], value='std_dev'))

        denom_values = get_ufloat(denom_tally, 'ifp-denominator')
        if gen_time_tally is None:
            generation_time = None
        else:
            gen_time_values = get_ufloat(gen_time_tally, 'ifp-time-numerator')
            gen_time_values /= denom_values*self.keff
            generation_time = gen_time_values.flatten()[0]

        if beta_tally is None:
            beta_effective = None
        else:
            beta_values = get_ufloat(beta_tally, 'ifp-beta-numerator')
            beta_values /= denom_values
            beta_effective = beta_values.flatten()
            if beta_effective.size == 1:
                beta_effective = beta_effective[0]

        return KineticsParameters(generation_time, beta_effective)
