#!/usr/bin/env python2

import struct
from collections import OrderedDict

import numpy as np
import scipy.stats

filter_types = {1: 'universe', 2: 'material', 3: 'cell', 4: 'cellborn',
                5: 'surface', 6: 'mesh', 7: 'energyin', 8: 'energyout'}

score_types = {-1: 'flux',
                -2: 'total',
                -3: 'scatter',
                -4: 'nu-scatter',
                -5: 'scatter-n',
                -6: 'scatter-pn',
                -7: 'transport',
                -8: 'n1n',
                -9: 'absorption',
                -10: 'fission',
                -11: 'nu-fission',
                -12: 'kappa-fission',
                -13: 'current',
                -14: 'events',
                1: '(n,total)',
                2: '(n,elastic)',
                4: '(n,level)',
                11: '(n,2nd)',
                16: '(n,2n)',
                17: '(n,3n)',
                18: '(n,fission)',
                19: '(n,f)',
                20: '(n,nf)',
                21: '(n,2nf)',
                22: '(n,na)',
                23: '(n,n3a)',
                24: '(n,2na)',
                25: '(n,3na)',
                28: '(n,np)',
                29: '(n,n2a)',
                30: '(n,2n2a)',
                32: '(n,nd)',
                33: '(n,nt)',
                34: '(n,nHe-3)',
                35: '(n,nd2a)',
                36: '(n,nt2a)',
                37: '(n,4n)',
                38: '(n,3nf)',
                41: '(n,2np)',
                42: '(n,3np)',
                44: '(n,n2p)',
                45: '(n,npa)',
                91: '(n,nc)',
                101: '(n,disappear)',
                102: '(n,gamma)',
                103: '(n,p)',
                104: '(n,d)',
                105: '(n,t)',
                106: '(n,3He)',
                107: '(n,a)',
                108: '(n,2a)',
                109: '(n,3a)',
                111: '(n,2p)',
                112: '(n,pa)',
                113: '(n,t2a)',
                114: '(n,d2a)',
                115: '(n,pd)',
                116: '(n,pt)',
                117: '(n,da)',
                201: '(n,Xn)',
                202: '(n,Xgamma)',
                203: '(n,Xp)',
                204: '(n,Xd)',
                205: '(n,Xt)',
                206: '(n,X3He)',
                207: '(n,Xa)',
                444: '(damage)',
                649: '(n,pc)',
                699: '(n,dc)',
                749: '(n,tc)',
                799: '(n,3Hec)',
                849: '(n,tc)'}
score_types.update({MT: '(n,n' + str(MT-50) + ')' for MT in range(51,91)})
score_types.update({MT: '(n,p' + str(MT-600) + ')' for MT in range(600,649)})
score_types.update({MT: '(n,d' + str(MT-650) + ')' for MT in range(650,699)})
score_types.update({MT: '(n,t' + str(MT-700) + ')' for MT in range(700,749)})
score_types.update({MT: '(n,3He' + str(MT-750) + ')' for MT in range(750,649)})
score_types.update({MT: '(n,a' + str(MT-800) + ')' for MT in range(800,849)})

class Mesh(object):
    def __init__(self):
        pass

    def __repr__(self):
        if hasattr(self, "dimension"):
            return "<Mesh: {0}>".format(tuple(self.dimension))
        else:
            return "<Mesh>"

class Filter(object):
    def __init__(self):
        self.type = 0
        self.bins = []

    def __repr__(self):
        return "<Filter: {0}>".format(self.type)

class Tally(object):
    def __init__(self):
        self.filters = OrderedDict()


class SourceSite(object):
    def __init__(self):
        self.weight = None
        self.xyz = None
        self.uvw = None
        self.E = None

    def __repr__(self):
        return "<SourceSite: xyz={0} at E={1}>".format(self.xyz, self.E)


class StatePoint(object):
    def __init__(self, filename):
        if filename.endswith('.h5'):
            import h5py
            self._f = h5py.File(filename, 'r')
            self._hdf5 = True
        else:
            self._f = open(filename, 'rb')
            self._hdf5 = False

        # Set flags for what data  was read
        self._metadata = False
        self._results = False
        self._source = False

        # Initialize arrays for meshes and tallies
        self.meshes = []
        self.tallies = []
        self.source = []

        # Read all metadata
        self._read_metadata()

    def _read_metadata(self):
        # Read filetype
        self.filetype = self._get_int(path='filetype')[0]

        # Read statepoint revision
        self.revision = self._get_int(path='revision')[0]
        if self.revision != 11:
          raise Exception('Statepoint Revision is not consistent.')

        # Read OpenMC version
        if self._hdf5:
            self.version = [self._get_int(path='version_major')[0],
                            self._get_int(path='version_minor')[0],
                            self._get_int(path='version_release')[0]]
        else:
            self.version = self._get_int(3)

        # Read date and time
        self.date_and_time = self._get_string(19, path='date_and_time')

        # Read path
        self.path = self._get_string(255, path='path').strip()

        # Read random number seed
        self.seed = self._get_long(path='seed')[0]

        # Read run information
        self.run_mode = self._get_int(path='run_mode')[0]
        self.n_particles = self._get_long(path='n_particles')[0]
        self.n_batches = self._get_int(path='n_batches')[0]

        # Read current batch
        self.current_batch = self._get_int(path='current_batch')[0]

        # Read criticality information
        if self.run_mode == 2:
            self.n_inactive = self._get_int(path='n_inactive')[0]
            self.gen_per_batch = self._get_int(path='gen_per_batch')[0]
            self.k_batch = self._get_double(
                self.current_batch*self.gen_per_batch, path='k_generation')
            self.entropy = self._get_double(
                self.current_batch*self.gen_per_batch, path='entropy')
            self.k_col_abs = self._get_double(path='k_col_abs')[0]
            self.k_col_tra = self._get_double(path='k_col_tra')[0]
            self.k_abs_tra = self._get_double(path='k_abs_tra')[0]
            self.k_combined = self._get_double(2, path='k_combined')

            # Read CMFD information
            cmfd_present = self._get_int(path='cmfd_on')[0]
            if cmfd_present == 1:
                self.cmfd_indices = self._get_int(4, path='cmfd/indices')
                self.k_cmfd = self._get_double(self.current_batch,
                              path='cmfd/k_cmfd')
                self.cmfd_src = self._get_double_array(np.product(self.cmfd_indices),
                                path='cmfd/cmfd_src')
                self.cmfd_src = np.reshape(self.cmfd_src,
                                tuple(self.cmfd_indices), order='F')
                self.cmfd_entropy = self._get_double(self.current_batch,
                                    path='cmfd/cmfd_entropy')
                self.cmfd_balance = self._get_double(self.current_batch,
                                    path='cmfd/cmfd_balance')
                self.cmfd_dominance = self._get_double(self.current_batch,
                                      path='cmfd/cmfd_dominance')
                self.cmfd_srccmp = self._get_double(self.current_batch,
                                   path='cmfd/cmfd_srccmp')

        # Read number of meshes
        n_meshes = self._get_int(path='tallies/n_meshes')[0]

        # Read meshes
        for i in range(n_meshes):
            m = Mesh()
            self.meshes.append(m)

            base = 'tallies/mesh' + str(i+1) + '/'

            # Read id, mesh type, and number of dimensions
            m.id = self._get_int(path=base+'id')[0]
            m.type = self._get_int(path=base+'type')[0]
            n = self._get_int(path=base+'n_dimension')[0]

            # Read mesh size, lower-left coordinates, upper-right coordinates,
            # and width of each mesh cell
            m.dimension = self._get_int(n, path=base+'dimension')
            m.lower_left = self._get_double(n, path=base+'lower_left')
            m.upper_right = self._get_double(n, path=base+'upper_right')
            m.width = self._get_double(n, path=base+'width')

        # Read number of tallies
        n_tallies = self._get_int(path='tallies/n_tallies')[0]

        for i in range(n_tallies):
            # Create Tally object and add to list of tallies
            t = Tally()
            self.tallies.append(t)

            base = 'tallies/tally' + str(i+1) + '/'

            # Read id and number of realizations
            t.id = self._get_int(path=base+'id')[0]
            t.n_realizations = self._get_int(path=base+'n_realizations')[0]

            # Read sizes of tallies
            t.total_score_bins = self._get_int(path=base+'total_score_bins')[0]
            t.total_filter_bins = self._get_int(path=base+'total_filter_bins')[0]

            # Read number of filters
            n_filters = self._get_int(path=base+'n_filters')[0]

            for j in range(n_filters):
                # Create Filter object
                f = Filter()

                base = 'tallies/tally{0}/filter{1}/'.format(i+1, j+1)

                # Get type of filter
                f.type = filter_types[self._get_int(path=base+'type')[0]]

                # Add to filter dictionary
                t.filters[f.type] = f

                # Determine how many bins are in this filter
                f.length = self._get_int(path=base+'n_bins')[0]
                assert f.length > 0
                if f.type == 'energyin' or f.type == 'energyout':
                    f.bins = self._get_double(f.length + 1, path=base+'bins')
                elif f.type == 'mesh':
                    f.bins = self._get_int(path=base+'bins')
                else:
                    f.bins = self._get_int(f.length, path=base+'bins')

            base = 'tallies/tally' + str(i+1) + '/'

            # Read nuclide bins
            n_nuclides = self._get_int(path=base+'n_nuclide_bins')[0]
            t.n_nuclides = n_nuclides
            t.nuclides = self._get_int(n_nuclides, path=base+'nuclide_bins')

            # Read score bins and scattering order
            t.n_scores = self._get_int(path=base+'n_score_bins')[0]
            t.scores = [score_types[j] for j in self._get_int(
                    t.n_scores, path=base+'score_bins')]
            t.scatt_order = self._get_int(t.n_scores, path=base+'scatt_order')

            # Read number of user score bins
            t.n_user_scores = self._get_int(path=base+'n_user_score_bins')[0]

            # Set up stride
            stride = 1
            for f in t.filters.values()[::-1]:
                f.stride = stride
                stride *= f.length

        # Source bank present
        source_present = self._get_int(path='source_present')[0]
        if source_present == 1:
            self.source_present = True       
        else:
            self.source_present = False

        # Set flag indicating metadata has already been read
        self._metadata = True

    def read_results(self):
        # Check whether metadata has been read
        if not self._metadata:
            self._read_metadata()

        # Number of realizations for global tallies
        self.n_realizations = self._get_int(path='n_realizations')[0]

        # Read global tallies
        n_global_tallies = self._get_int(path='n_global_tallies')[0]
        if self._hdf5:
            data = self._f['global_tallies'].value
            self.global_tallies = np.column_stack((data['sum'], data['sum_sq']))
        else:
            self.global_tallies = np.array(self._get_double(2*n_global_tallies))
            self.global_tallies.shape = (n_global_tallies, 2)

        # Flag indicating if tallies are present
        tallies_present = self._get_int(path='tallies/tallies_present')[0]

        # Read tally results
        if tallies_present:
            for i, t in enumerate(self.tallies):
                n = t.total_score_bins * t.total_filter_bins
                if self._hdf5:
                    path = 'tallies/tally{0}/results'.format(i+1)
                    data = self._f[path].value
                    t.results = np.column_stack((data['sum'], data['sum_sq']))
                    t.results.shape = (t.total_filter_bins, t.total_score_bins, 2)
                else:
                    t.results = np.array(self._get_double(2*n))
                    t.results.shape = (t.total_filter_bins, t.total_score_bins, 2)

        # Indicate that tally results have been read
        self._results = True

    def read_source(self):
        # Check whether tally results have been read
        if not self._results:
            self.read_results()

        # Check if source bank is in statepoint
        if not self.source_present:
            print('Source not in statepoint file.')
            return

        # For HDF5 state points, copy entire bank
        if self._hdf5:
            source_sites = self._f['source_bank'].value

        for i in range(self.n_particles):
            s = SourceSite()
            self.source.append(s)

            # Read position, angle, and energy
            if self._hdf5:
                s.weight, s.xyz, s.uvw, s.E = source_sites[i]
            else:
                s.weight = self._get_double()[0]
                s.xyz = self._get_double(3)
                s.uvw = self._get_double(3)
                s.E = self._get_double()[0]

    def generate_ci(self, confidence=0.95):
        """Calculates confidence intervals for each tally bin."""

        # Determine number of realizations
        n = self.n_realizations

        # Determine significance level and percentile for two-sided CI
        alpha = 1 - confidence
        percentile = 1 - alpha/2

        # Calculate t-value
        t_value = scipy.stats.t.ppf(percentile, n - 1)
        self.generate_stdev(t_value)

    def generate_stdev(self, t_value=1.0):
        """
        Calculates the sample mean and standard deviation of the mean for each
        tally bin.
        """

        # Determine number of realizations
        n = self.n_realizations

        # Global tallies
        for i in range(len(self.global_tallies)):
            # Get sum and sum of squares
            s, s2 = self.global_tallies[i]

            # Calculate sample mean and replace value
            s /= n
            self.global_tallies[i,0] = s

            # Calculate standard deviation
            if s != 0.0:
                self.global_tallies[i,1] = t_value*np.sqrt((s2/n - s*s)/(n-1))

        # Regular tallies
        for t in self.tallies:
            for i in range(t.results.shape[0]):
                for j in range(t.results.shape[1]):
                    # Get sum and sum of squares
                    s, s2 = t.results[i,j]

                    # Calculate sample mean and replace value
                    s /= n
                    t.results[i,j,0] = s

                    # Calculate standard deviation
                    if s != 0.0:
                        t.results[i,j,1] = t_value*np.sqrt((s2/n - s*s)/(n-1))

    def get_value(self, tally_index, spec_list, score_index):
        """Returns a tally score given a list of filters to satisfy.

        Parameters
        ----------
        tally_index : int
            Index for tally in StatePoint.tallies list

        spec_list : list
            A list of tuples where the first value in each tuple is the filter
            type, e.g. 'cell', and the second value is the desired index. If the
            first value in the tuple is 'mesh', the second value should be a
            tuple with three integers specifying the mesh indices.

            Example: [('cell', 1), ('mesh', (14,17,20)), ('energyin', 2)]

        score_index : int
            Index corresponding to score for tally, i.e. the second index in
            Tally.results[:,:,:].

        """

        # Get Tally object given the index
        t = self.tallies[tally_index]

        # Initialize index for filter in Tally.results[:,:,:]
        filter_index = 0

        # Loop over specified filters in spec_list
        for f_type, f_index in spec_list:

            # Treat mesh filter separately
            if f_type == 'mesh':
                # Get index in StatePoint.meshes
                mesh_index = t.filters['mesh'].bins[0] - 1

                # Get dimensions of corresponding mesh
                nx, ny, nz = self.meshes[mesh_index].dimension

                # Convert (x,y,z) to a single bin -- this is similar to
                # subroutine mesh_indices_to_bin in openmc/src/mesh.F90.
                value = ((f_index[0] - 1)*ny*nz +
                         (f_index[1] - 1)*nz +
                         (f_index[2] - 1))
                filter_index += value*t.filters[f_type].stride
            else:
                filter_index += f_index*t.filters[f_type].stride

        # Return the desired result from Tally.results. This could be the sum and
        # sum of squares, or it could be mean and stdev if self.generate_stdev()
        # has been called already.
        return t.results[filter_index, score_index]

    def extract_results(self, tally_id, score_str):
        """Returns a tally results dictionary given a tally_id and score string.

           Parameters
           ----------
           tally_id : int
               Index for the tally in StatePoint.tallies list

           score_str : string
               Corresponds to the string entered for a score in tallies.xml.
               For a flux score extraction it would be 'score'

        """

        # get tally
        try:
            tally = self.tallies[tally_id-1]
        except:
            print 'Tally does not exist'
            return

        # get the score index if it is present
        try:
            idx = tally.scores.index(score_str)
        except ValueError:
            print 'Score does not exist'
            print tally.scores
            return

        # create numpy array for mean and 95% CI
        n_bins = len(tally.results)
        n_filters = len(tally.filters)
        n_scores = len(tally.scores)
        meanv = np.zeros(n_bins)
        unctv = np.zeros(n_bins)
        filters = np.zeros((n_bins,n_filters))
        filtmax = np.zeros(n_filters+1)
        meshmax = np.zeros(4)
        filtmax[0] = 1
        meshmax[0] = 1

        # get number of realizations
        n = tally.n_realizations

        # get t-value
        t_value = scipy.stats.t.ppf(0.975, n - 1)

        # calculate mean
        meanv = tally.results[:,idx,0]
        meanv = meanv / n

        # calculate 95% two-sided CI
        unctv = tally.results[:,idx,1]
        unctv = t_value*np.sqrt((unctv/n - meanv*meanv)/(n-1))/meanv

        # create output dictionary
        data = {'mean':meanv,'CI95':unctv}

        # get bounds of filter bins
        for akey in tally.filters.keys():
            idx = tally.filters.keys().index(akey)
            filtmax[n_filters - idx] = tally.filters[akey].length

        # compute bin info
        for i in range(n_filters):

            # compute indices for filter combination
            filters[:,n_filters - i - 1] = np.floor((np.arange(n_bins) %
                   np.prod(filtmax[0:i+2]))/(np.prod(filtmax[0:i+1]))) + 1

            # append in dictionary bin with filter
            data.update({tally.filters.keys()[n_filters - i - 1]:
                         filters[:,n_filters - i - 1]})

            # check for mesh
            if tally.filters.keys()[n_filters - i - 1] == 'mesh':
                dims = list(self.meshes[tally.filters['mesh'].bins[0] - 1].dimension)
                dims.reverse()
                dims = np.asarray(dims)
                if score_str == 'current':
                    dims += 1
                meshmax[1:4] = dims
                mesh_bins = np.zeros((n_bins,3))
                mesh_bins[:,2] = np.floor(((filters[:,n_filters - i - 1] - 1) %
                            np.prod(meshmax[0:2]))/(np.prod(meshmax[0:1]))) + 1
                mesh_bins[:,1] = np.floor(((filters[:,n_filters - i - 1] - 1) %
                            np.prod(meshmax[0:3]))/(np.prod(meshmax[0:2]))) + 1
                mesh_bins[:,0] = np.floor(((filters[:,n_filters - i - 1] - 1) %
                            np.prod(meshmax[0:4]))/(np.prod(meshmax[0:3]))) + 1
                data.update({'mesh':zip(mesh_bins[:,0],mesh_bins[:,1],
                            mesh_bins[:,2])})
            i += 1

        # add in maximum bin filters and order
        b = tally.filters.keys()
        b.reverse()
        filtmax = list(filtmax[1:])
        try:
            idx = b.index('mesh')
            filtmax[idx] = np.max(mesh_bins[:,2])
            filtmax.insert(idx,np.max(mesh_bins[:,1]))
            filtmax.insert(idx,np.max(mesh_bins[:,0]))

        except ValueError:
            pass
        data.update({'bin_order':b,'bin_max':filtmax})

        return data

    def _get_data(self, n, typeCode, size):
        return list(struct.unpack('={0}{1}'.format(n,typeCode),
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
