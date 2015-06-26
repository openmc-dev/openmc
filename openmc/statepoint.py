import copy
import struct
import sys

import numpy as np
import scipy.stats

import openmc
from openmc.constants import *

if sys.version > '3':
    long = int


class SourceSite(object):
    """A single source site produced from fission.

    Attributes
    ----------
    weight : float
        Weight of the particle arising from the site
    xyz : list of float
        Cartesian coordinates of the site
    uvw : list of float
        Directional cosines for particles emerging from the site
    E : float
        Energy of the emerging particle in MeV

    """

    def __init__(self):
        self._weight = None
        self._xyz = None
        self._uvw = None
        self._E = None

    def __repr__(self):
        string = 'SourceSite\n'
        string += '{0: <16}{1}{2}\n'.format('\tweight', '=\t', self._weight)
        string += '{0: <16}{1}{2}\n'.format('\tE', '=\t', self._E)
        string += '{0: <16}{1}{2}\n'.format('\t(x,y,z)', '=\t', self._xyz)
        string += '{0: <16}{1}{2}\n'.format('\t(u,v,w)', '=\t', self._uvw)
        return string

    @property
    def weight(self):
        return self._weight

    @property
    def xyz(self):
        return self._xyz

    @property
    def uvw(self):
        return self._uvw

    @property
    def E(self):
        return self._E


class StatePoint(object):
    """State information on a simulation at a certain point in time (at the end of a
    given batch). Statepoints can be used to analyze tally results as well as
    restart a simulation.

    Attributes
    ----------
    k_combined : list
        Combined estimator for k-effective and its uncertainty
    n_particles : int
        Number of particles per generation
    n_batches : int
        Number of batches
    current_batch :
        Number of batches simulated
    results : bool
        Indicate whether tally results have been read
    source : ndarray of SourceSite
        Array of source sites
    with_summary : bool
        Indicate whether statepoint data has been linked against a summary file
    tallies : dict
        Dictionary whose keys are tally IDs and whose values are Tally objects
    tallies_present : bool
        Indicate whether user-defined tallies are present
    global_tallies : ndarray
        Global tallies and their uncertainties
    n_realizations : int
        Number of tally realizations

    """

    def __init__(self, filename):
        if filename.endswith('.h5'):
            import h5py
            self._f = h5py.File(filename, 'r')
            self._hdf5 = True
        else:
            self._f = open(filename, 'rb')
            self._hdf5 = False

        # Set flags for what data has been read
        self._results = False
        self._source = False
        self._with_summary = False

        # Read all metadata
        self._read_metadata()

        # Read information about tally meshes
        self._read_meshes()

        # Read tally metadata
        self._read_tallies()

    def close(self):
        self._f.close()

    @property
    def k_combined(self):
        return self._k_combined

    @property
    def n_particles(self):
        return self._n_particles

    @property
    def n_batches(self):
        return self._n_batches

    @property
    def current_batch(self):
        return self._current_batch

    @property
    def results(self):
        return self._results

    @property
    def source(self):
        return self._source

    @property
    def with_summary(self):
        return self._with_summary

    @property
    def tallies(self):
        return self._tallies

    @property
    def tallies_present(self):
        return self._tallies_present

    @property
    def global_tallies(self):
        return self._global_tallies

    @property
    def n_realizations(self):
        return self._n_realizations

    def _read_metadata(self):
        # Read filetype
        self._filetype = self._get_int(path='filetype')[0]

        # Read statepoint revision
        self._revision = self._get_int(path='revision')[0]
        if self._revision != 13:
            raise Exception('Statepoint Revision is not consistent.')

        # Read OpenMC version
        if self._hdf5:
            self._version = [self._get_int(path='version_major')[0],
                             self._get_int(path='version_minor')[0],
                             self._get_int(path='version_release')[0]]
        else:
            self._version = self._get_int(3)

        # Read date and time
        self._date_and_time = self._get_string(19, path='date_and_time')

        # Read path
        self._path = self._get_string(255, path='path').strip()

        # Read random number seed
        self._seed = self._get_long(path='seed')[0]

        # Read run information
        self._run_mode = self._get_int(path='run_mode')[0]
        self._n_particles = self._get_long(path='n_particles')[0]
        self._n_batches = self._get_int(path='n_batches')[0]

        # Read current batch
        self._current_batch = self._get_int(path='current_batch')[0]

        # Read whether or not the source site distribution is present
        self._source_present = self._get_int(path='source_present')[0]

        # Read criticality information
        if self._run_mode == 2:
            self._read_criticality()

    def _read_criticality(self):
        # Read criticality information
        if self._run_mode == 2:

            self._n_inactive = self._get_int(path='n_inactive')[0]
            self._gen_per_batch = self._get_int(path='gen_per_batch')[0]
            self._k_batch = self._get_double(
                 self._current_batch*self._gen_per_batch,
                 path='k_generation')
            self._entropy = self._get_double(
                 self._current_batch*self._gen_per_batch, path='entropy')

            self._k_col_abs = self._get_double(path='k_col_abs')[0]
            self._k_col_tra = self._get_double(path='k_col_tra')[0]
            self._k_abs_tra = self._get_double(path='k_abs_tra')[0]
            self._k_combined = self._get_double(2, path='k_combined')

            # Read CMFD information (if used)
            self._read_cmfd()

    def _read_cmfd(self):
        base = 'cmfd'

        # Read CMFD information
        self._cmfd_on = self._get_int(path='cmfd_on')[0]

        if self._cmfd_on == 1:

            self._cmfd_indices = self._get_int(4, path='{0}/indices'.format(base))
            self._k_cmfd = self._get_double(self._current_batch,
                 path='{0}/k_cmfd'.format(base))
            self._cmfd_src = self._get_double_array(np.product(self._cmfd_indices),
                 path='{0}/cmfd_src'.format(base))
            self._cmfd_src = np.reshape(self._cmfd_src, tuple(self._cmfd_indices),
                  order='F')
            self._cmfd_entropy = self._get_double(self._current_batch,
                  path='{0}/cmfd_entropy'.format(base))
            self._cmfd_balance = self._get_double(self._current_batch,
                  path='{0}/cmfd_balance'.format(base))
            self._cmfd_dominance = self._get_double(self._current_batch,
                  path='{0}/cmfd_dominance'.format(base))
            self._cmfd_srccmp = self._get_double(self._current_batch,
                  path='{0}/cmfd_srccmp'.format(base))

    def _read_meshes(self):
        # Initialize dictionaries for the Meshes
        # Keys     - Mesh IDs
        # Values - Mesh objects
        self._meshes = {}

        # Read the number of Meshes
        self._n_meshes = self._get_int(path='tallies/meshes/n_meshes')[0]

        # Read a list of the IDs for each Mesh
        if self._n_meshes > 0:

            # OpenMC Mesh IDs (redefined internally from user definitions)
            self._mesh_ids = self._get_int(self._n_meshes,
                 path='tallies/meshes/ids')

            # User-defined Mesh IDs
            self._mesh_keys = self._get_int(self._n_meshes,
                 path='tallies/meshes/keys')

        else:
            self._mesh_keys = []
            self._mesh_ids = []

        # Build dictionary of Meshes
        base = 'tallies/meshes/mesh '

        # Iterate over all Meshes
        for mesh_key in self._mesh_keys:

            # Read the user-specified Mesh ID and type
            mesh_id = self._get_int(path='{0}{1}/id'.format(base, mesh_key))[0]
            mesh_type = self._get_int(path='{0}{1}/type'.format(base, mesh_key))[0]

            # Get the Mesh dimension
            n_dimension = self._get_int(
                 path='{0}{1}/n_dimension'.format(base, mesh_key))[0]

            # Read the mesh dimensions, lower-left coordinates,
            # upper-right coordinates, and width of each mesh cell
            dimension = self._get_int(
                 n_dimension, path='{0}{1}/dimension'.format(base, mesh_key))
            lower_left = self._get_double(
                 n_dimension, path='{0}{1}/lower_left'.format(base, mesh_key))
            upper_right = self._get_double(
                 n_dimension, path='{0}{1}/upper_right'.format(base, mesh_key))
            width = self._get_double(
                 n_dimension, path='{0}{1}/width'.format(base, mesh_key))

            # Create the Mesh and assign properties to it
            mesh = openmc.Mesh(mesh_id)

            mesh.dimension = dimension
            mesh.width = width
            mesh.lower_left = lower_left
            mesh.upper_right = upper_right

            #FIXME: Set the mesh type to 'rectangular' by default
            mesh.type = 'rectangular'

            # Add mesh to the global dictionary of all Meshes
            self._meshes[mesh_id] = mesh

    def _read_tallies(self):
        # Initialize dictionaries for the Tallies
        # Keys     - Tally IDs
        # Values   - Tally objects
        self._tallies = {}

        # Read the number of tallies
        self._n_tallies = self._get_int(path='/tallies/n_tallies')[0]

        # Read a list of the IDs for each Tally
        if self._n_tallies > 0:

            # OpenMC Tally IDs (redefined internally from user definitions)
            self._tally_ids = self._get_int(
                 self._n_tallies, path='tallies/ids')

            # User-defined Tally IDs
            self._tally_keys = self._get_int(
                 self._n_tallies, path='tallies/keys')

        else:
            self._tally_keys = []
            self._tally_ids = []

        base = 'tallies/tally '

        # Iterate over all Tallies
        for tally_key in self._tally_keys:

            # Read integer Tally estimator type code (analog or tracklength)
            estimator_type = self._get_int(
                 path='{0}{1}/estimator'.format(base, tally_key))[0]

            # Read the Tally size specifications
            n_realizations = self._get_int(
                 path='{0}{1}/n_realizations'.format(base, tally_key))[0]

            # Create Tally object and assign basic properties
            tally = openmc.Tally(tally_key)
            tally.estimator = ESTIMATOR_TYPES[estimator_type]
            tally.num_realizations = n_realizations

            # Read the number of Filters
            n_filters = self._get_int(
                 path='{0}{1}/n_filters'.format(base, tally_key))[0]

            subbase = '{0}{1}/filter '.format(base, tally_key)

            # Initialize all Filters
            for j in range(1, n_filters+1):

                # Read the integer Filter type code
                filter_type = self._get_int(
                     path='{0}{1}/type'.format(subbase, j))[0]

                # Read the Filter offset
                offset = self._get_int(
                     path='{0}{1}/offset'.format(subbase, j))[0]

                n_bins = self._get_int(
                     path='{0}{1}/n_bins'.format(subbase, j))[0]

                if n_bins <= 0:
                    msg = 'Unable to create Filter {0} for Tally ID={2} ' \
                          'since no bins were specified'.format(j, tally_key)
                    raise ValueError(msg)

                # Read the bin values
                if FILTER_TYPES[filter_type] in ['energy', 'energyout']:
                    bins = self._get_double(
                         n_bins+1, path='{0}{1}/bins'.format(subbase, j))

                elif FILTER_TYPES[filter_type] in ['mesh', 'distribcell']:
                    bins = self._get_int(
                         path='{0}{1}/bins'.format(subbase, j))[0]

                else:
                    bins = self._get_int(
                         n_bins, path='{0}{1}/bins'.format(subbase, j))

                # Create Filter object
                filter = openmc.Filter(FILTER_TYPES[filter_type], bins)
                filter.offset = offset
                filter.num_bins = n_bins

                if FILTER_TYPES[filter_type] == 'mesh':
                    key = self._mesh_keys[self._mesh_ids.index(bins)]
                    filter.mesh = self._meshes[key]

                # Add Filter to the Tally
                tally.add_filter(filter)

            # Read Nuclide bins
            n_nuclides = self._get_int(
                 path='{0}{1}/n_nuclides'.format(base, tally_key))[0]

            nuclide_zaids = self._get_int(
                 n_nuclides, path='{0}{1}/nuclides'.format(base, tally_key))

            # Add all Nuclides to the Tally
            for nuclide_zaid in nuclide_zaids:
                tally.add_nuclide(nuclide_zaid)

            # Read score bins
            n_score_bins = self._get_int(
                 path='{0}{1}/n_score_bins'.format(base, tally_key))[0]

            tally.num_score_bins = n_score_bins

            scores = [SCORE_TYPES[j] for j in self._get_int(
                 n_score_bins, path='{0}{1}/score_bins'.format(base, tally_key))]
            n_user_scores = self._get_int(
                 path='{0}{1}/n_user_score_bins'.format(base, tally_key))[0]

            # Compute and set the filter strides
            for i in range(n_filters):
                filter = tally.filters[i]
                filter.stride = n_score_bins * n_nuclides

                for j in range(i+1, n_filters):
                    filter.stride *= tally.filters[j].num_bins

            # Read scattering moment order strings (e.g., P3, Y-1,2, etc.)
            moments = []
            subbase = '{0}{1}/moments/'.format(base, tally_key)

            # Extract the moment order string for each score
            for k in range(len(scores)):
                moment = self._get_string(8,
                     path='{0}order{1}'.format(subbase, k+1))
                moment = moment.lstrip('[\'')
                moment = moment.rstrip('\']')

                # Remove extra whitespace
                moment.replace(" ", "")
                moments.append(moment)

            # Add the scores to the Tally
            for j, score in enumerate(scores):
                # If this is a scattering moment, insert the scattering order
                if '-n' in score:
                    score = score.replace('-n', '-' + str(moments[j]))
                elif '-pn' in score:
                    score = score.replace('-pn', '-' + str(moments[j]))
                elif '-yn' in score:
                    score = score.replace('-yn', '-' + str(moments[j]))

                tally.add_score(score)

            # Add Tally to the global dictionary of all Tallies
            self.tallies[tally_key] = tally

    def read_results(self):
        """Read tally results and store them in the ``tallies`` attribute. No results
        are read when the statepoint is instantiated.

        """

        # Number of realizations for global Tallies
        self._n_realizations = self._get_int(path='n_realizations')[0]

        # Read global Tallies
        n_global_tallies = self._get_int(path='n_global_tallies')[0]

        if self._hdf5:
            data = self._f['global_tallies'].value
            self._global_tallies = np.column_stack((data['sum'], data['sum_sq']))

        else:
            self._global_tallies = np.array(self._get_double(2*n_global_tallies))
            self._global_tallies.shape = (n_global_tallies, 2)

        # Flag indicating if Tallies are present
        self._tallies_present = self._get_int(path='tallies/tallies_present')[0]

        base = 'tallies/tally '

        # Read Tally results
        if self._tallies_present:

            # Iterate over and extract the results for all Tallies
            for tally_key in self._tally_keys:

                # Get this Tally
                tally = self._tallies[tally_key]

                # Compute the total number of bins for this Tally
                num_tot_bins = tally.num_bins

                # Extract Tally data from the file
                if self._hdf5:
                    data = self._f['{0}{1}/results'.format(base, tally_key)].value
                    sum = data['sum']
                    sum_sq = data['sum_sq']

                else:
                    results = np.array(self._get_double(2*num_tot_bins))
                    sum = results[0::2]
                    sum_sq = results[1::2]

                # Define a routine to convert 0 to 1
                def nonzero(val):
                    return 1 if not val else val

                # Reshape the results arrays
                new_shape = (nonzero(tally.num_filter_bins),
                             nonzero(tally.num_nuclides),
                             nonzero(tally.num_score_bins))

                sum = np.reshape(sum, new_shape)
                sum_sq = np.reshape(sum_sq, new_shape)

                # Set the data for this Tally
                tally.sum = sum
                tally.sum_sq = sum_sq

        # Indicate that Tally results have been read
        self._results = True

    def read_source(self):
        """Read and store source sites from the statepoint file. By default, source
        sites are not loaded upon initialization.

        """

        # Check whether Tally results have been read
        if not self._results:
            self.read_results()

        # Check if source bank is in statepoint
        if not self._source_present:
            print('Unable to read source since it is not in statepoint file')
            return

        # Initialize a NumPy array for the source sites
        self._source = np.empty(self._n_particles, dtype=SourceSite)

        # For HDF5 state points, copy entire bank
        if self._hdf5:
            source_sites = self._f['source_bank'].value

        # Initialize SourceSite object for each particle
        for i in range(self._n_particles):
            # Initialize new source site
            site = SourceSite()

            # Read position, angle, and energy
            if self._hdf5:
                site._weight, site._xyz, site._uvw, site._E = source_sites[i]
            else:
                site._weight = self._get_double()[0]
                site._xyz = self._get_double(3)
                site._uvw = self._get_double(3)
                site._E = self._get_double()[0]

            # Store the source site in the NumPy array
            self._source[i] = site

    def compute_ci(self, confidence=0.95):
        """Computes confidence intervals for each Tally bin.

        This method is equivalent to calling compute_stdev(...) when the
        confidence is known as opposed to its corresponding t value.

        Parameters
        ----------
        confidence : float, optional
            Confidence level. Defaults to 0.95.

        """

        # Determine significance level and percentile for two-sided CI
        alpha = 1 - confidence
        percentile = 1 - alpha/2

        # Calculate t-value
        t_value = scipy.stats.t.ppf(percentile, self._n_realizations - 1)
        self.compute_stdev(t_value)

    def compute_stdev(self, t_value=1.0):
        """Computes the sample mean and the standard deviation of the mean
        for each Tally bin.

        Parameters
        ----------
        t_value : float, optional
            Student's t-value applied to the uncertainty. Defaults to 1.0,
            meaning the reported value is the sample standard deviation.

        """

        # Determine number of realizations
        n = self._n_realizations

        # Calculate the standard deviation for each global tally
        for i in range(len(self._global_tallies)):

            # Get sum and sum of squares
            s, s2 = self._global_tallies[i]

            # Calculate sample mean and replace value
            s /= n
            self._global_tallies[i, 0] = s

            # Calculate standard deviation
            if s != 0.0:
                self._global_tallies[i, 1] = t_value * np.sqrt((s2 / n - s**2) / (n-1))

        # Calculate sample mean and standard deviation for user-defined Tallies
        for tally_id, tally in self.tallies.items():
            tally.compute_std_dev(t_value)

    def get_tally(self, scores=[], filters=[], nuclides=[],
                  name=None, id=None, estimator=None):
        """Finds and returns a Tally object with certain properties.

        This routine searches the list of Tallies and returns the first Tally
        found it finds which satisfies all of the input parameters.
        NOTE: The input parameters do not need to match the complete Tally
        specification and may only represent a subset of the Tally's properties.

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
        id : int, optional
            The id specified for the Tally (default is None).
        estimator: str, optional
            The type of estimator ('tracklength', 'analog'; default is None).

        Returns
        -------
        tally : Tally
            A tally matching the specified criteria

        Raises
        ------
        LookupError
            If a Tally meeting all of the input parameters cannot be found in
            the statepoint.

        """

        tally = None

        # Iterate over all tallies to find the appropriate one
        for tally_id, test_tally in self.tallies.items():

            # Determine if Tally has queried name
            if name and name != test_tally.name:
                continue

            # Determine if Tally has queried id
            if id and id != test_tally.id:
                continue

            # Determine if Tally has queried estimator
            if estimator and not estimator == test_tally.estimator:
                continue

            # Determine if Tally has the queried score(s)
            if scores:
                contains_scores = True

                # Iterate over the scores requested by the user
                for score in scores:
                    if score not in test_tally.scores:
                        contains_scores = False
                        break

                if not contains_scores:
                    continue

            # Determine if Tally has the queried Filter(s)
            if filters:
                contains_filters = True

                # Iterate over the Filters requested by the user
                for filter in filters:
                    if filter not in test_tally.filters:
                        contains_filters = False
                        break

                if not contains_filters:
                    continue

            # Determine if Tally has the queried Nuclide(s)
            if nuclides:
                contains_nuclides = True

                # Iterate over the Nuclides requested by the user
                for nuclide in nuclides:
                    if nuclide not in test_tally.nuclides:
                        contains_nuclides = False
                        break

                if not contains_nuclides:
                    continue

            # If the current Tally met user's request, break loop and return it
            tally = test_tally
            break

        # If we did not find the Tally, return an error message
        if tally is None:
            raise LookupError('Unable to get Tally')

        return tally

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
        summary : Summary
             A Summary object.

        Raises
        ------
        ValueError
            An error when the argument passed to the 'summary' parameter is not
            an openmc.Summary object.

        """

        if not isinstance(summary, openmc.summary.Summary):
            msg = 'Unable to link statepoint with {0} which ' \
                  'is not a Summary object'.format(summary)
            raise ValueError(msg)

        for tally_id, tally in self.tallies.items():
            # Get the Tally name from the summary file
            tally.name = summary.tallies[tally_id].name
            tally.with_summary = True

            nuclide_zaids = copy.deepcopy(tally.nuclides)

            for nuclide_zaid in nuclide_zaids:
                tally.remove_nuclide(nuclide_zaid)
                if nuclide_zaid == -1:
                    tally.add_nuclide(openmc.Nuclide('total'))
                else:
                    tally.add_nuclide(summary.nuclides[nuclide_zaid])

            for filter in tally.filters:
                if filter.type == 'surface':
                    surface_ids = []
                    for bin in filter.bins:
                        surface_ids.append(summary.surfaces[bin].id)
                    filter.bins = surface_ids

                if filter.type in ['cell', 'distribcell']:
                    distribcell_ids = []
                    for bin in filter.bins:
                        distribcell_ids.append(summary.cells[bin].id)
                    filter.bins = distribcell_ids

                if filter.type == 'universe':
                    universe_ids = []
                    for bin in filter.bins:
                        universe_ids.append(summary.universes[bin].id)
                    filter.bins = universe_ids

                if filter.type == 'material':
                    material_ids = []
                    for bin in filter.bins:
                        material_ids.append(summary.materials[bin].id)
                    filter.bins = material_ids

        self._with_summary = True

    def _get_data(self, n, typeCode, size):
        return list(struct.unpack('={0}{1}'.format(n, typeCode),
                    self._f.read(n*size)))

    def _get_int(self, n=1, path=None):
        if self._hdf5:
            return [int(v) for v in self._f[path].value]
        else:
            return [int(v) for v in self._get_data(n, 'i', 4)]

    def _get_long(self, n=1, path=None):
        if self._hdf5:
            return [long(v) for v in self._f[path].value]
        else:
            return [long(v) for v in self._get_data(n, 'q', 8)]

    def _get_float(self, n=1, path=None):
        if self._hdf5:
            return [float(v) for v in self._f[path].value]
        else:
            return [float(v) for v in self._get_data(n, 'f', 4)]

    def _get_double(self, n=1, path=None):
        if self._hdf5:
            return [float(v) for v in self._f[path].value]
        else:
            return [float(v) for v in self._get_data(n, 'd', 8)]

    def _get_double_array(self, n=1, path=None):
        if self._hdf5:
            return self._f[path].value
        else:
            return self._get_data(n, 'd', 8)

    def _get_string(self, n=1, path=None):
        if self._hdf5:
            return str(self._f[path].value)
        else:
            return str(self._get_data(n, 's', 1)[0])
