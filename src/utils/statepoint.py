#!/usr/bin/env python

import struct
from math import sqrt

import numpy as np
import scipy.stats

filter_types = {1: 'universe', 2: 'material', 3: 'cell', 4: 'cellborn',
                5: 'surface', 6: 'mesh', 7: 'energyin', 8: 'energyout'}

score_types = {-1: 'flux', -2: 'total', -3: 'scatter', -4: 'nu-scatter', 
               -5: 'scatter-1', -6: 'scatter-2', -7: 'scatter-3', 
               -8: 'transport', -9: 'diffusion', -10: 'n1n', -11: 'n2n',
               -12: 'n3n', -13: 'n4n', -14: 'absorption', -15: 'fission',
                -16: 'nu-fission', -17: 'current'}

class BinaryFile(object):
    def __init__(self, filename):
        self._f = open(filename, 'rb')

    def _get_data(self, n, typeCode, size):
        return list(struct.unpack('={0}{1}'.format(n,typeCode),
                                  self._f.read(n*size)))
    
    def _get_int(self, n=1):
        return self._get_data(n, 'i', 4)

    def _get_long(self, n=1):
        return self._get_data(n, 'q', 8)

    def _get_float(self, n=1):
        return self._get_data(n, 'f', 4)

    def _get_double(self, n=1):
        return self._get_data(n, 'd', 8)

    def _get_string(self, length, n=1):
        data = self._get_data(length*n, 's', 1)[0]
        return [data[i*length:(i+1)*length] for i in range(n)]

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
        self.filters = []


class SourceSite(object):
    def __init__(self):
        self.weight = None
        self.xyz = None
        self.uvw = None
        self.E = None

    def __repr__(self):
        return "<SourceSite: xyz={0} at E={1}>".format(self.xyz, self.E)


class StatePoint(BinaryFile):
    def __init__(self, filename):
        super(StatePoint, self).__init__(filename)

        # Set flags for what data  was read
        self._metadata = False
        self._values = False
        self._source = False

        # Initialize arrays for meshes and tallies
        self.meshes = []
        self.tallies = []
        self.source = []

        # Read all metadata
        self._read_metadata()

    def _read_metadata(self):
        # Read statepoint revision
        self.revision = self._get_int()[0]

        # Read OpenMC version
        self.version = self._get_int(3)

        # Read date and time
        self.date_and_time = self._get_string(19)[0]

        # Read random number seed
        self.seed = self._get_long()[0]

        # Read run information
        self.run_mode = self._get_int()[0]
        self.n_particles = self._get_long()[0]
        self.n_batches = self._get_int()[0]

        # Read current batch
        self.current_batch = self._get_int()[0]

        # Read criticality information
        if self.run_mode == 2:
            self.n_inactive, self.gen_per_batch = self._get_int(2)
            self.k_batch = self._get_double(self.current_batch)
            self.entropy = self._get_double(self.current_batch)

        # Read global tallies
        n_global_tallies = self._get_int()[0]
        self.global_tallies = np.array(self._get_double(2*n_global_tallies))
        self.global_tallies.shape = (n_global_tallies, 2)

        # Read number of meshes
        n_meshes = self._get_int()[0]

        # Read meshes
        for i in range(n_meshes):
            m = Mesh()
            self.meshes.append(m)

            # Read type of mesh and number of dimensions
            m.type = self._get_int()[0]
            n = self._get_int()[0]

            # Read mesh size, lower-left coordinates, upper-right coordinates,
            # and width of each mesh cell
            m.dimension = self._get_int(n)
            m.lower_left = self._get_double(n)
            m.upper_right = self._get_double(n)
            m.width = self._get_double(n)

        # Read number of tallies
        n_tallies = self._get_int()[0]

        for i in range(n_tallies):
            # Create Tally object and add to list of tallies
            t = Tally()
            self.tallies.append(t)

            # Read sizes of tallies
            t.n_score_bins = self._get_int()[0]
            t.n_filter_bins = self._get_int()[0]

            # Read number of filters
            n_filters = self._get_int()[0]

            for j in range(n_filters):
                # Create Filter object and add to list of filters
                f = Filter()
                t.filters.append(f)

                # Get type of filter
                f.type = filter_types[self._get_int()[0]]

                # Determine how many bins are in this filter
                f.length = self._get_int()[0]
                assert f.length > 0
                if f.type == 'energyin' or f.type == 'energyout':
                    f.bins = self._get_double(f.length + 1)
                elif f.type == 'mesh':
                    f.bins = self._get_int()
                else:
                    f.bins = self._get_int(f.length)
            
            # Read nuclide bins
            n_nuclides = self._get_int()[0]
            t.nuclides = self._get_int(n_nuclides)

            # Read score bins
            n_scores = self._get_int()[0]
            t.scores = [score_types[j] for j in self._get_int(n_scores)]

        # Set flag indicating metadata has already been read
        self._metadata = True

    def read_values(self):
        # Check whether metadata has been read
        if not self._metadata:
            self._read_metadata()

        # Flag indicating if tallies are present
        tallies_present = self._get_int()[0]

        # Read tally results
        if tallies_present:
            for t in self.tallies:
                n = t.n_score_bins * t.n_filter_bins
                t.values = np.array(self._get_double(2*n))
                t.values.shape = (t.n_filter_bins, t.n_score_bins, 2)

        # Indicate that tally values have been read
        self._values = True

    def read_source(self):
        # Check whether tally values have been read
        if not self._values:
            self.read_values()

        for i in range(self.n_particles):
            s = SourceSite()
            self.source.append(s)

            # Read position, angle, and energy
            s.weight = self._get_double()[0]
            s.xyz = self._get_double(3)
            s.uvw = self._get_double(3)
            s.E = self._get_double()[0]

    def generate_ci(self, confidence=0.95):
        """Calculates confidence intervals for each tally bin."""

        # Determine number of realizations
        if self.run_mode == 2:
            n = self.current_batch - self.n_inactive
        else:
            n = self.current_batch

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
        if self.run_mode == 2:
            n = self.current_batch - self.n_inactive
        else:
            n = self.current_batch

        for t in self.tallies:
            for i in range(t.values.shape[0]):
                for j in range(t.values.shape[1]):
                    # Get sum and sum of squares
                    s, s2 = t.values[i,j]
                    
                    # Calculate sample mean and replace value
                    s /= n
                    t.values[i,j,0] = s

                    # Calculate standard deviation
                    if s != 0.0:
                        t.values[i,j,1] = t_value*sqrt((s2/n - s*s)/(n-1))
