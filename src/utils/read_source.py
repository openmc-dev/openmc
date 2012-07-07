#!/usr/bin/env python

import struct

class SourceFile(object):
    def __init__(self, filename):
        # Open source file for reading
        self.f = open(filename, 'r')

        # Read number of source sites
        self.n_sites = struct.unpack('q', self.f.read(8))[0]

        # Create list to store source sites
        self.sites = []

        for i in range(self.n_sites):
            # Read position, angle, and energy
            (uid,) = struct.unpack('q', self.f.read(8))
            (weight,) = struct.unpack('d', self.f.read(8))
            xyz = struct.unpack('3d', self.f.read(24))
            uvw = struct.unpack('3d', self.f.read(24))
            (E,) = struct.unpack('d', self.f.read(8))

            # Create source site and append to list
            self.sites.append(SourceSite(weight, xyz, uvw, E))

class SourceSite(object):
    def __init__(self, weight, xyz, uvw, E):
        self.weight = weight
        self.xyz = xyz
        self.uvw = uvw
        self.E = E

    def __repr__(self):
        return "<SourceSite: xyz={0} at E={1}>".format(self.xyz, self.E)
