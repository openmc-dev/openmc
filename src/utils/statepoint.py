#!/usr/bin/env python

import struct

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
    
    def get_int(self, n=1):
        return self._get_data(n, 'i', 4)

    def get_long(self, n=1):
        return self._get_data(n, 'q', 8)

    def get_float(self, n=1):
        return self._get_data(n, 'f', 4)

    def get_double(self, n=1):
        return self._get_data(n, 'd', 8)

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

class StatePoint(BinaryFile):
    def __init__(self, filename):
        super(StatePoint, self).__init__(filename)

        # Initialize arrays for meshes and tallies
        self.meshes = []
        self.tallies = []

        # Read all metadata
        self._read_metadata()

    def _read_metadata(self):
        # Read statepoint revision
        self.revision = self.get_int()[0]

        # Read OpenMC version
        self.version = self.get_int(3)

        # Read random number seed
        self.seed = self.get_long()[0]

        # Read run information
        self.run_mode = self.get_int()[0]
        self.n_particles = self.get_long()[0]
        self.n_batches, self.n_inactive, self.gen_per_batch = self.get_int(3)

        # Read current batch
        self.current_batch = self.get_int()[0]

        # Read batch keff and entropy
        keff = self.get_double(self.current_batch)
        entropy = self.get_double(self.current_batch)

        # Read global tallies
        self.n_global_tallies = self.get_int()[0]
        self.global_tallies = self.get_double(2*self.n_global_tallies)

        # Read number of meshes
        n_meshes = self.get_int()[0]

        # Read meshes
        for i in range(n_meshes):
            m = Mesh()
            self.meshes.append(m)

            # Read type of mesh and number of dimensions
            m.type = self.get_int()[0]
            n = self.get_int()[0]

            # Read mesh size, lower-left coordinates, upper-right coordinates,
            # and width of each mesh cell
            m.dimension = self.get_int(n)
            m.lower_left = self.get_double(n)
            m.upper_right = self.get_double(n)
            m.width = self.get_double(n)


        # Read number of tallies
        n_tallies = self.get_int()[0]

        for i in range(n_tallies):
            # Create Tally object and add to list of tallies
            t = Tally()
            self.tallies.append(t)

            # Read sizes of tallies
            t.n_score_bins = self.get_int()[0]
            t.n_filter_bins = self.get_int()[0]

            # Read number of filters
            n_filters = self.get_int()[0]

            for j in range(n_filters):
                # Create Filter object and add to list of filters
                f = Filter()
                t.filters.append(f)

                # Get type of filter
                f.type = filter_types[self.get_int()[0]]

                # Determine how many bins are in this filter
                f.length = self.get_int()[0]
                assert f.length > 0
                if f.type == 'energyin' or f.type == 'energyout':
                    f.bins = self.get_double(f.length + 1)
                elif f.type == 'mesh':
                    f.bins = self.get_int()
                else:
                    f.bins = self.get_int(f.length)
            
            # Read nuclide bins
            n_nuclides = self.get_int()[0]
            t.nuclides = self.get_int(n_nuclides)

            # Read score bins
            n_scores = self.get_int()[0]
            t.scores = [score_types[j] for j in self.get_int(n_scores)]
