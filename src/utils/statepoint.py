#!/usr/bin/env python2

from sys import exit
import struct
from collections import OrderedDict

import numpy as np
import scipy.stats

filter_types = {1: 'universe', 2: 'material', 3: 'cell',
                4: 'cellborn', 5: 'surface', 6: 'mesh',
                7: 'energyin', 8: 'energyout', 9: 'distribcell'}

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
        self.offset = -1

    def __repr__(self):
        return "<Filter: {0}>".format(self.type)


class Tally(object):

    def __init__(self):

        self.filters = OrderedDict()
        self.id = None
	self.label = None
	self.estimator = None
        self.n_realizations = None
        self.total_score_bins = None
        self.total_filter_bins = None

    def getScoreIndex(self, score):
        return self.scores.index(score)


class SourceSite(object):
    def __init__(self):
        self.weight = None
        self.xyz = None
        self.uvw = None
        self.E = None

    def __repr__(self):
        return "<SourceSite: xyz={0} at E={1}>".format(self.xyz, self.E)

class Geometry_Data(object):
    def __init__(self):
        self.univ = []
        self.cell = []
        self.lat  = []
        self.surf = []
        self.mat  = []
        self.n_cells = 0
        self.n_lattices = 0
        self.n_universes = 0
        self.n_surfaces = 0
        self.n_materials = 0

    def _print_all(self):

        print 'Geometry Information:'
        print 'n_cells:',self.n_cells
        print 'n_universes:',self.n_universes
        print 'n_lattices:',self.n_lattices
        print 'n_surfaces:',self.n_surfaces
        print 'n_materials:',self.n_materials
	print ''
	print 'Universes:'
	print ''
        for i in self.univ:
          print 'Universe:',i.ID
          print '--Contains Cells:',i.cells

	print ''
	print ''	
	print 'Cells:'
	print ''
        for i in self.cell:
          print 'Cell:',i.ID
          print '--UserID:',i.userID
          print '--Filltype:',i.filltype
          print '--Fill:'
          print i.fill
          print '--Offset:',i.offset
          print '--Material:',i.material

	print ''
	print ''	
	print 'Lattices:'
	print ''
        for i in self.lat:
          print 'Lattice:',i.ID
          print '--Fill:',i.fill
          if type(i.offset) == int:
            pass
          else:
            print '--Offset Dimensions:',i.offset.shape
          print '--Offset:',np.squeeze(i.offset)
          print '--Dimenisions:',i.dim

	print ''
	print ''	
	print 'Surfaces:'
	print ''
        for i in self.surf:
          print 'Surface:',i.ID
          print '--Type:',i.s_type
          print '--Coefficients:'
          print i.coeffs
	  print '--# Positive Neighbors:',i.n_pos
          print '--Positive Neighbors:'
	  print i.neighbor_pos
	  print '--# Negative Neighbors:',i.n_neg
          print '--Negative Neighbors:'
	  print i.neighbor_neg
          print 'Boundary Condition:',i.bc


	print ''
	print ''	
	print 'Materials:'
	print ''
        for i in self.mat:
          print 'Material:',i.ID
          print '--# Nuclides:',i.n_nuclide
          print '--Nuclides:'
	  print i.nuclide
          print '--Density:',i.density
          print '--Atom Densities:'
	  print i.atom_density
          print '--# Sab Tables:',i.n_sab
          print '--Sab Nuclides:'
	  print i.i_sab_nuclides
          print '--Tables:'
	  print i.i_sab_tables

    def _get_offset(self, path, filter_offset):
        """
        Returns the corresponding location in the results array for a given
        path and filter number.

        Parameters
        ----------
        path : list 
            A list of IDs that form the path to the target. It should begin with 0
            for the base universe, and should cover every universe, cell, and lattice passed through.
            For the case of the lattice, a tuple should be provided to indicate which coordinates
            in the lattice should be entered. This should be in the form: (lat_id, i_x, i_y, i_z)

        filter_offset : int
            An integer that specifies which offset map the filter is using

        """
        prev = -1
        prevtype = ''
        # The variable returned at completion
        offset = 0
        # Verify path length at least 2. 1 for base universe, 1 so it could end in a cell
        if len(path) < 2:
          error = "ERROR: len(path) < 2. It cannot meet all requirements."
          exit(error)
        # Iterate over path
        for i in path:

          if prev == -1:
            # We should get the base universe first
            prev = i
            prevtype = 'U'
            if prev != 0:
              error = "ERROR: First element in path was not 0."
              exit(error)

          else:
            if isinstance(i,tuple):
              if len(i) != 4:
                error = "ERROR: Tuple " + str(i) + " does not have exactly 4 elements."
                exit(error)

              # WILL'S CHANGE
              if i[0] != self.cell[self.cellKeys.index(prev)].fill:
                error = "ERROR: Previous cell did not contain lattice " + str(i[0]) + "."
                exit(error)

              # Find the lattice in the lattice list
              for j in self.lat:
                if j.ID == i[0]:
                  if (i[1] > j.dim[0] or i[1] < 1 or 
                      i[2] > j.dim[1] or i[2] < 1 or
                      i[3] > j.dim[2] or i[3] < 1):
                    error = "Bad lattice index specified for lattice "+str(i[0])+"."
                    exit(error)
                  offset += j.offset[i[3]-1,i[2]-1,i[1]-1,filter_offset-1]
                  prev = i[0]
                  prevtype = 'L'
                  break
              else:
                # Couldn't find the lattice - Not good 
                    error = "ERROR: Could not find lattice data for ID " + str(i[0]) + "."
                    exit(error)             
            
            else:
              # Expect universe or cell or lattice (as final target)
              if prevtype == 'U':
                # Last construct was a universe
                # Need to set j to be the OpenMC ID of the cell, as that is what the universes store
                for c in self.cell:
                  if i == c.userID:
                    jTrue = c.ID
                l = 0
                for u in self.univ:
                  if u.ID == prev:                    
                    break
                  l += 1
                if jTrue not in self.univ[l].cells:
                    error = "ERROR: Requested cell " + str(i) + " from universe ",str(prev)," was not found."
                    exit(error)
                # Get key index
                prev = i
                prevtype = 'C'
              elif prevtype == 'C':                
                # Last construct was a cell 
                # Check if the previous cell was a normal cell
                # Throw Exception if it was, there's no deeper cells
                for c in self.cell:
#                  if (prev == c.userID and c.filltype == "normal").any():
                  if prev == c.userID and c.filltype == 'normal':
                    error = "ERROR: Cell " + str(c.userID)  + " is normal cell, cannot contain any lower levels."
                    exit(error)
                # Else check if fill or lattice
                # If fill, just go to that universe
                for c in self.cell:
                  if prev == c.userID:
                    if c.filltype == "fill_cell":
                      if c.fill != i[0]:
                        error = "ERROR: Fill Cell " + str(c.userID) + " does not contain universe " + str(i[0])  + "."
                        exit(error)
                      prevtype = 'U'
                      prev = i[0]      
                      offset += c.offset[filter_offset-1]                   
                    elif c.filltype == "lattice":
                      prev = i[0]
                      prevtype = 'L'    
              elif prevtype == 'L':
                # Last construct was a lattice 
                # Double check that this universe is contained by that lattice
                for j in self.lat:
                  if j.ID == prev:
                    if not (j.fill == i).any():
                      error = "ERROR: Lattice " + str(prev) + " does not contain universe " + str(i) + "."
                      exit(error)
                # Found a lattice that matched and contains this universe
                prevtype = 'U'
                prev = i  
        return offset


class Universe(object):

    def __init__(self, ID, cells):

        self.ID = ID
        self.cells = cells


class Cell(object):

    def __init__(self, ID, userID, filltype):

        self.ID = ID
        self.userID = userID
        self.filltype = filltype
        self.offset = []
        self.material = -1
        self.fill = -1


    def _set_material(self, material):
        self.material = material


    def _set_fill(self, fill):
        self.fill = fill


    def _set_offset(self, offset):
        self.offset = offset


class Lattice(object):

    def __init__(self, ID, fill, offset, dim):

        self.ID = ID
        self.fill = fill
        self.offset = offset
        self.dim = dim
        if len(self.dim) == 2:
            self.dim.append(1)


class Surface(object):

    def __init__(self, ID, s_type, coeffs, n_pos, neighbor_pos, n_neg, neighbor_neg, bc):

        self.ID = ID
        self.s_type = s_type
        self.coeffs = coeffs
        self.n_pos = n_pos
        self.neighbor_pos = neighbor_pos
        self.n_neg = n_neg
        self.neighbor_neg = neighbor_neg
        self.bc = bc


class Material(object):

    def __init__(self, ID, n_nuclide, nuclide, density, atom_density, n_sab, 
                 i_sab_nuclides, i_sab_tables):

        self.ID = ID
        self.n_nuclide = n_nuclide
        self.nuclide = nuclide
        self.density = density
        self.atom_density = atom_density
        self.n_sab = n_sab
        self.i_sab_nuclides = i_sab_nuclides
        self.i_sab_tables = i_sab_tables


class StatePoint(object):

    def __init__(self, filename):

        self.tallies_present = False

        if filename.endswith('.h5'):
            import h5py
            self._f = h5py.File(filename, 'r')
            self._hdf5 = True
        else:
            self._f = open(filename, 'rb')
            self._hdf5 = False

        # Tally ID -> Index Dictionary
        self.tallyID = dict()

        # Set flags for what data was read
        self._metadata = False
        self._results = False
        self._source = False

        # Initialize arrays for meshes and tallies
        self.meshes = []
        self.tallies = []
        self.source = []

        # Initialize geometry data
        self.geom = Geometry_Data()

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

	#for i in xrange(300):
	#  print i,':',self._get_int(path='geometry')

        # Read geometry information
        self.geom.n_cells = self._get_int(path='geometry/n_cells')[0]
        self.geom.n_universes = self._get_int(path='geometry/n_universes')[0]
        self.geom.n_lattices = self._get_int(path='geometry/n_lattices')[0]
        self.geom.n_surfaces = self._get_int(path='geometry/n_surfaces')[0]
        self.geom.n_materials = self._get_int(path='geometry/n_materials')[0]

        if self.geom.n_lattices > 0:
          latticeList = self._get_int(self.geom.n_lattices,path='geometry/lattice_ids')

        univList = self._get_int(self.geom.n_universes,path='geometry/universe_ids')

        surfList = self._get_int(self.geom.n_surfaces,path='geometry/surface_ids')

        matList = self._get_int(self.geom.n_materials,path='geometry/material_ids')

        # User Inputs  
        self.geom.cellKeys = self._get_int(self.geom.n_cells,path='geometry/cell_keys')
        # OpenMC IDs
        self.geom.cellList = self._get_int(self.geom.n_cells,path='geometry/cell_ids')

        # Build list of cells
        base = 'geometry/cells/cell '
        for i in range(self.geom.n_cells):
          filltypeInt = self._get_int(path=base + str(self.geom.cellKeys[i])+'/fill_type')[0]
          if filltypeInt == 1:
            filltype = "normal"
          elif filltypeInt == 2:
            filltype = "fill cell"
          else:
            filltype = "lattice"
          self.geom.cell.append(Cell(self.geom.cellList[i],self.geom.cellKeys[i],filltype))
          if filltype == "normal":
            material = self._get_int(path=base + str(self.geom.cellKeys[i])+'/material')[0]
            self.geom.cell[-1]._set_material(material)
          else:
            if filltype == "lattice":
              fill = self._get_int(path=base + str(self.geom.cellKeys[i])+'/lattice')[0]
            elif filltype == "fill cell":
              fill = self._get_int(path=base + str(self.geom.cellKeys[i])+'/fill')[0]
              maps = self._get_int(path=base + str(self.geom.cellKeys[i])+'/maps')[0]
              if maps > 0:
                offset = self._get_int(path=base + str(self.geom.cellKeys[i])+'/offset')
                self.geom.cell[-1]._set_offset(offset)
            self.geom.cell[-1]._set_fill(fill)

        # Build list of universes
        base = 'geometry/universes/universe '
        for i in range(self.geom.n_universes):
          n = self._get_int(path=base + str(univList[i])+'/cells')[0]
          if n != 0:
            self.geom.univ.append(Universe(univList[i],self._get_int(n, path=base + str(univList[i])+'/cells')))

        # Build list of lattices
        base = 'geometry/lattices/lattice '
        for i in range(self.geom.n_lattices):
          struct = self._get_int(path=base + str(latticeList[i])+'/type')
          dim = self._get_int(3,path=base + str(latticeList[i])+'/dimension')
          maps = self._get_int(path=base + str(latticeList[i])+'/maps')[0]
          offset_size = self._get_int(path=base + str(latticeList[i])+'/offset_size')[0]
          offset = 0
          if offset_size > 0:
            if self._hdf5:
              path = base + str(latticeList[i])+'/offset'
              data = self._f[path].value
              offset = data    
              path = base + str(latticeList[i])+'/universes'
              data = self._f[path].value
              fill = data
            else:
              offset = np.array(self._get_int(dim[0]*dim[1]*dim[2]*maps))
              offset.shape = (dim[2], dim[1], dim[0],maps)     
              fill = np.array(self._get_int(dim[0]*dim[1]*dim[2]))
              fill.shape = (dim[0], dim[1], dim[2]) 
          else:
            if self._hdf5:    
              path = base + str(latticeList[i])+'/universes'
              data = self._f[path].value
              fill = data
            else:
              fill = np.array(self._get_int(dim[0]*dim[1]*dim[2]))
              fill.shape = (dim[0], dim[1], dim[2]) 
          self.geom.lat.append(Lattice(latticeList[i],fill, offset, dim))

	# Build list of surfaces
	base = 'geometry/surfaces/surface '
	for i in range(self.geom.n_surfaces):
          s_type = self._get_int(path=base + str(surfList[i])+'/type')[0]
	  bc = self._get_int(path=base + str(surfList[i])+'/bc')[0]
	  n_coeff = self._get_int(path=base + str(surfList[i])+'/n_coeffs')[0]
          coeffs = self._get_double(n_coeff,path=base + str(surfList[i])+'/coeffs')
	  n_pos = self._get_int(path=base + str(surfList[i])+'/n_neighbor_pos')[0]
	  if n_pos > 0:
            pos_n = self._get_int(n_pos,path=base + str(surfList[i])+'/neighbor_pos')
	  else:
	    pos_n = None
	  n_neg = self._get_int(path=base + str(surfList[i])+'/n_neighbor_neg')[0]
	  if n_neg > 0:
            neg_n = self._get_int(n_neg,path=base + str(surfList[i])+'/neighbor_neg')
	  else:
	    neg_n = None
          self.geom.surf.append(Surface(surfList[i], s_type, coeffs, n_pos, pos_n, n_neg, neg_n, bc))

	# Build list of materials
	base = 'geometry/materials/material '
	for i in range(self.geom.n_materials):
          n_nuclide = self._get_int(path=base + str(matList[i])+'/n_nuclides')[0]
          nuclide = self._get_int(n_nuclide,path=base + str(matList[i])+'/nuclide')
	  dens = self._get_double(path=base + str(matList[i])+'/density')[0]
	  a_dens = self._get_double(n_nuclide,path=base + str(matList[i])+'/atom_density')
	  n_sab = self._get_int(path=base + str(matList[i])+'/n_sab')[0]
	  sab_n = self._get_int(n_sab,path=base + str(matList[i])+'/i_sab_nuclides')
	  sab_t = self._get_int(n_sab,path=base + str(matList[i])+'/i_sab_tables')
	
          self.geom.mat.append(Material(matList[i], n_nuclide, nuclide, dens, a_dens, n_sab, sab_n, sab_t))

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

	#print 'PRINTING'
        for i in range(n_tallies):
            # Create Tally object and add to list of tallies
            t = Tally()
            self.tallies.append(t)

            base = 'tallies/tally' + str(i+1) + '/'

            # Read id and number of realizations
            t.id = self._get_int(path=base+'id')[0]
            labelsize = self._get_int(path=base+'id')[0]
	    t.label = self._get_string(labelsize,path=base+'id')
	    t.estimator = self._get_int(path=base+'id')[0]
            t.n_realizations = self._get_int(path=base+'n_realizations')[0]

            # Read sizes of tallies
            t.total_score_bins = self._get_int(path=base+'total_score_bins')[0]
            print 't.total_score_bins:',t.total_score_bins
            t.total_filter_bins = self._get_int(path=base+'total_filter_bins')[0]
            print 't.total_filter_bins:',t.total_filter_bins

            # Add tally to dictionary
            self.tallyID[t.id] = i
#            self.tallyID[i] = t.id

            # Read number of filters
            n_filters = self._get_int(path=base+'n_filters')[0]
            for j in range(n_filters):
                # Create Filter object
                f = Filter()

                base = 'tallies/tally{0}/filter{1}/'.format(i+1, j+1)

                # Get type of filter
                f.type = filter_types[self._get_int(path=base+'type')[0]]

                # Get offset of filter
                f.offset = self._get_int(path=base+'offset')[0]

                # Add to filter dictionary
                t.filters[f.type] = f

                # Determine how many bins are in this filter
                f.length = self._get_int(path=base+'n_bins')[0]
                assert f.length > 0
                if f.type == 'energyin' or f.type == 'energyout':
                    f.bins = self._get_double(f.length + 1, path=base+'bins')
                elif f.type == 'mesh' or f.type == 'distribcell':
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
        self.tallies_present = self._get_int(path='tallies/tallies_present')[0]

        # Read tally results
        if self.tallies_present:
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

    def get_value(self, tally_ID, spec_list, score_index):
        """Returns a tally score given a list of filters to satisfy.

        Parameters
        ----------
        tally_ID : int
            The ID of the tally as specified in the openMC input file

        spec_list : list
            A list of tuples where the first value in each tuple is the filter
            type, e.g. 'cell', and the second value is the desired index. If the
            first value in the tuple is 'mesh', the second value should be a
            tuple with three integers specifying the mesh indices.

            Example: [('cell', 1), ('mesh', (14,17,20)), ('energyin', 2)]
            Example: [('distribcell', path)] or [('distribcell', 3)]

        score_index : int
            Index corresponding to score for tally, i.e. the second index in
            Tally.results[:,:,:].

        """

        # Get Tally object given the index
        t = self.tallies[self.tallyID[tally_ID]]

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
            elif f_type == 'distribcell':
                if type(f_index) != int:
                    filter_offset = t.filters['distribcell'].offset
                    value = self.geom._get_offset(f_index, filter_offset)
                    filter_index += value*t.filters[f_type].stride          
                else:
                    filter_index += f_index*t.filters[f_type].stride
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
               The ID of the tally as specified in the openMC input file

           score_str : string
               Corresponds to the string entered for a score in tallies.xml.
               For a flux score extraction it would be 'score'

        """

        # get tally
        try:
            tally = self.tallies[self.tallyID[tally_id]]
        except:
            print 'Tally does not exist'
            return

        # get the score index if it is present
        try:
            idx = tally.scores.index(score_str)
        except ValueError:
            print 'Score does not exist'
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
