#!/usr/bin/env python

import struct
import numpy as np
import scipy.stats
import openmc
from openmc.constants import *



class SourceSite(object):

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



class StatePoint(object):

  def __init__(self, filename):

    openmc.reset_auto_ids()

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

    # Read all metadata
    self._read_metadata()

    # Read information needed to reconstruct the Materials and Geometry
    self._read_geometry()

    # Read information about tally meshes
    self._read_meshes()

    # Read tally metadata
    self._read_tallies()


  def _read_metadata(self):

    # Read filetype
    self._filetype = self._get_int(path='metadata/filetype')[0]

    # Read statepoint revision
    self._revision = self._get_int(path='metadata/revision')[0]
    if self._revision != 12:
      raise Exception('Statepoint Revision is not consistent.')

    # Read OpenMC version
    if self._hdf5:
      self._version = [self._get_int(path='metadata/version_major')[0],
                       self._get_int(path='metadata/version_minor')[0],
                       self._get_int(path='metadata/version_release')[0]]
    else:
      self._version = self._get_int(3)

    # Read date and time
    self._date_and_time = self._get_string(19, path='metadata/date_and_time')

    # Read path
    self._path = self._get_string(255, path='metadata/path').strip()

    # Read random number seed
    self._seed = self._get_long(path='metadata/seed')[0]

    # Read run information
    self._run_mode = self._get_int(path='metadata/run_mode')[0]
    self._n_particles = self._get_long(path='metadata/n_particles')[0]
    self._n_batches = self._get_int(path='metadata/n_batches')[0]

    # Read current batch
    self._current_batch = self._get_int(path='metadata/current_batch')[0]

    # Read whether or not the source site distribution is present
    self._source_present = self._get_int(path='metadata/source_present')[0]

    # Read criticality information
    if self._run_mode == 2:
      self._read_criticality()


  def _read_criticality(self):

    # Read criticality information
    if self._run_mode == 2:

      self._n_inactive = self._get_int(path='criticality/n_inactive')[0]
      self._gen_per_batch = self._get_int(path='criticality/gen_per_batch')[0]
      self._k_batch = self._get_double(
           self._current_batch*self._gen_per_batch,
           path='criticality/k_generation')
      self._entropy = self._get_double(
           self._current_batch*self._gen_per_batch, path='criticality/entropy')

      self._k_col_abs = self._get_double(path='criticality/k_col_abs')[0]
      self._k_col_tra = self._get_double(path='criticality/k_col_tra')[0]
      self._k_abs_tra = self._get_double(path='criticality/k_abs_tra')[0]
      self._k_combined = self._get_double(2, path='criticality/k_combined')

      # Read CMFD information (if used)
      self._read_cmfd()


  def _read_cmfd(self):

    base = 'criticality/cmfd'

    # Read CMFD information
    self._cmfd_on = self._get_int(path='{0}/cmfd_on'.format(base))[0]

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


  def _read_geometry(self):

    # Read in and initialize the Materials and Geometry
    self._read_nuclides()
    self._read_materials()
    self._read_surfaces()
    self._read_cells()
    self._read_universes()
    self._read_lattices()
    self._finalize_geometry()



  def _read_nuclides(self):

    # Initialize dictionary for each Nuclide
    # Keys   - Nuclide IDs
    # Values - Nuclide objects
    self._nuclides = dict()

    self._n_nuclides = self._get_int(path='model/nuclides/n_nuclides')[0]

    # OpenMC Nuclide IDs (defined internally)
    self._nuclide_ids = self._get_int(
         self._n_nuclides, path='model/nuclides/ids')

    # User-defined Nuclide string names
    self._nuclide_keys = dict()
    for i in range(self._n_nuclides):
      nuclide_name = self._get_string(11,
           path='model/nuclides/keys/{0}'.format(i+1))

      # Remove leading and trailing characters from string name
      nuclide_name = nuclide_name.lstrip('[\'')
      nuclide_name = nuclide_name.rstrip('\']')

      self._nuclide_keys[i+1] = nuclide_name

    # Iterate over and initalize all Nuclides
    for i, nuclide_id in enumerate(self._nuclide_ids):

      # Read the Nuclide's name (e.g., 'H-1' or 'U-235')
      name = self._nuclide_keys[i+1].split('.')[0]

      # Read the Nuclide's cross-section identifier (e.g., '70c')
      xs = self._nuclide_keys[i+1].split('.')[1]

      # Initialize this Nuclide and add it to the global dictionary of Nuclides
      self._nuclides[nuclide_id] = openmc.Nuclide(name=name, xs=xs)



  def _read_materials(self):

    # Initialize dictionary for each Material
    # Keys   - Material keys
    # Values - Material objects
    self._materials = dict()

    self._n_materials = self._get_int(path='model/materials/n_materials')[0]

    # OpenMC Material IDs (redefined internally from user definitions)
    self._material_ids = self._get_int(
         self._n_materials, path='model/materials/ids')

    # User-defined Material IDs
    self._material_keys = self._get_int(
         self._n_materials, path='model/materials/keys')

    # Build dictionary of Materials
    base = 'model/materials/material '

    # Iterate over all Materials
    for material_key in self._material_keys:

      # Read the number of Nuclides
      n_nuclides = self._get_int(
           path='{0}{1}/n_nuclides'.format(base, material_key))[0]

      # Read the Nuclide IDs
      nuclide_ids = self._get_int(
           n_nuclides, path='{0}{1}/nuclide_ids'.format(base, material_key))

      # Read the Material's density
      density = self._get_double(
           path='{0}{1}/density'.format(base, material_key))[0]

      # Read the atom/b-cm densities for each Nuclide
      nuclide_densities = self._get_double(
           n_nuclides, path='{0}{1}/nuclide_densities'.format(base, material_key))

      # Read the number of nuclides with S(a,b) tables
      n_sab = self._get_int(path='{0}{1}/n_sab'.format(base, material_key))[0]

      sab_names = list()
      sab_xs = list()

      # Read the names of the S(a,b) tables for this Material
      for i in range(1, n_sab+1):

        sab_table = self._get_string(12,
             path='{0}{1}/sab_tables/{2}'.format(base, material_key, i))

        # Read the cross-section identifiers for each S(a,b) table
        sab_names.append(sab_table.split('.')[0])
        sab_xs.append(sab_table.split('.')[1])

      # Create the Material
      material = openmc.Material(material_id=material_key)

      # Set the Material's density to g/cm3 - this is what is used in OpenMC
      material.setDensity(density=density, units='g/cm3')

      # Add all Nuclides to the Material
      for i, nuclide_id in enumerate(nuclide_ids):
        nuclide = self._nuclides[nuclide_id]
        density = nuclide_densities[i]
        material.addNuclide(nuclide, percent=density, percent_type='ao')

      # Add S(a,b) table(s?) to the Material
      for i in range(n_sab):
        name = sab_names[i]
        xs = sab_xs[i]
        material.setSAlphaBeta(name, xs)

      # Add the Material to the global dictionary of all Materials
      self._materials[material_key] = material


  def _read_surfaces(self):

    # Initialize dictionary for each Surface
    # Keys   - Surface keys
    # Values - Surfacee objects
    self._surfaces = dict()

    self._n_surfaces = self._get_int(path='model/surfaces/n_surfaces')[0]

    # OpenMC Surface IDs (redefined internally from user definitions)
    self._surface_ids = self._get_int(
         self._n_surfaces, path='model/surfaces/ids')

    # User-defined Surface IDs
    self._surface_keys = self._get_int(
         self._n_surfaces, path='model/surfaces/keys')

    #Build dictionary of Surfaces
    base = 'model/surfaces/surface '

    # Iterate over all Surfaces
    for surface_key in self._surface_keys:

      # Read the Surface type and boundary condition type codes
      type = self._get_int(path='{0}{1}/type'.format(base, surface_key))[0]
      bc = self._get_int(path='{0}{1}/bc'.format(base, surface_key))[0]

      # Read the Surface coefficients
      n_coeffs = self._get_int(path='{0}{1}/n_coeffs'.format(base, surface_key))[0]
      coeffs = self._get_double(
           n_coeffs, path='{0}{1}/coeffs'.format(base, surface_key))

      # Create the Surface based on it's type

      if SURFACE_TYPES[type] == 'x-plane':
        x0 = coeffs[0]
        surface = openmc.XPlane(surface_key, BC_TYPES[bc], x0)

      elif SURFACE_TYPES[type] == 'y-plane':
        y0 = coeffs[0]
        surface = openmc.YPlane(surface_key, BC_TYPES[bc], y0)

      elif SURFACE_TYPES[type] == 'z-plane':
        z0 = coeffs[0]
        surface = openmc.ZPlane(surface_key, BC_TYPES[bc], z0)

      elif SURFACE_TYPES[type] == 'plane':
        A = coeffs[0]
        B = coeffs[1]
        C = coeffs[2]
        D = coeffs[3]
        surface = openmc.Plane(surface_key, BC_TYPES[bc], A, B, C, D)

      elif SURFACE_TYPES[type] == 'x-cylinder':
        y0 = coeffs[0]
        z0 = coeffs[1]
        R = coeffs[2]
        surface = openmc.XCylinder(surface_key, BC_TYPES[bc], y0, z0, R)

      elif SURFACE_TYPES[type] == 'y-cylinder':
        x0 = coeffs[0]
        z0 = coeffs[1]
        R = coeffs[2]
        surface = openmc.YCylinder(surface_key, BC_TYPES[bc], x0, z0, R)

      elif SURFACE_TYPES[type] == 'z-cylinder':
        x0 = coeffs[0]
        y0 = coeffs[1]
        R = coeffs[2]
        surface = openmc.ZCylinder(surface_key, BC_TYPES[bc], x0, y0, R)

      elif SURFACE_TYPES[type] == 'sphere':
        x0 = coeffs[0]
        y0 = coeffs[1]
        z0 = coeffs[2]
        R = coeffs[3]
        surface = openmc.Sphere(surface_key, BC_TYPES[bc], x0, y0, z0, R)

      elif SURFACE_TYPES[type] in ['x-cone', 'y-cone', 'z-cone']:
        x0 = coeffs[0]
        y0 = coeffs[1]
        z0 = coeffs[2]
        R2 = coeffs[3]

        if SURFACE_TYPES[type] == 'x-cone':
          surface = openmc.XCone(surface_key, BC_TYPES[bc], x0, y0, z0, R2)
        if SURFACE_TYPES[type] == 'y-cone':
          surface = openmc.YCone(surface_key, BC_TYPES[bc], x0, y0, z0, R2)
        if SURFACE_TYPES[type] == 'z-cone':
          surface = openmc.ZCone(surface_key, BC_TYPES[bc], x0, y0, z0, R2)

      # Add Surface to global dictionary of all Surfaces
      self._surfaces[surface_key] = surface


  def _read_cells(self):

    # Initialize dictionary for each Cell
    # Keys   - Cell keys
    # Values - Cell objects
    self._cells = dict()

    self._n_cells = self._get_int(path='model/cells/n_cells')[0]

    # OpenMC Cell IDs (redefined internally from user definitions)
    self._cell_ids = self._get_int(self._n_cells, path='model/cells/ids')

    # User-defined Cell IDs
    self._cell_keys = self._get_int(self._n_cells, path='model/cells/keys')

    # Build dictionary of Cells
    base = 'model/cells/cell '

    # This is used later when the Universes and Lattices have been created
    # Keys    - Cell ID (user-defined)
    # Values  - (fill type string, fill ID) tuple
    self._cell_fills = dict()

    # Iterate over all Cells
    for cell_key in self._cell_keys:

      # Create this Cell
      cell = openmc.Cell(cell_id=cell_key)

      # Read the integer fill type code
      fill_type = self._get_int(
           path='{0}{1}/fill_type'.format(base, cell_key))[0]

      if FILL_TYPES[fill_type] == 'normal':
        fill_id = self._get_int(
             path='{0}{1}/material'.format(base, cell_key))[0]

      elif FILL_TYPES[fill_type] == 'fill':
        fill_id = self._get_int(path='{0}{1}/fill'.format(base, cell_key))[0]
        maps = self._get_int(path='{0}{1}/maps'.format(base, cell_key))[0]

        if maps > 0:
          cell.setOffset(self._get_int(
               maps, path='{0}{1}/offset'.format(base, cell_key)))

      else:
        fill_id = self._get_int(path='{0}{1}/lattice'.format(base, cell_key))[0]

      # Store the Cell fill information for use after Universe/Lattice creation
      self._cell_fills[cell_key] = (FILL_TYPES[fill_type], fill_id)

      # Read the number of Surfaces and Surface IDs
      n_surfaces = self._get_int(
        path='{0}{1}/n_surfaces'.format(base, cell_key))[0]

      # If the Cell contains Surfaces, add them to the Cell
      if n_surfaces > 0:

        surfaces = self._get_int(
          n_surfaces, path='{0}{1}/surface_ids'.format(base, cell_key))

        # Iterate over all Surfaces and add them to the Cell
        for surface_key in surfaces:

          # Determine Surface halfspace
          halfspace = np.sign(surface_key)

          surface_key = np.abs(surface_key)

          # Get the user-defined Sell ID and the corresponding Sell
          surface_index = self._surface_ids.index(surface_key)
          surface_id = self._surface_keys[surface_index]
          surface = self._surfaces[surface_id]

          # Add this Surface, halsfpace to the Cell
          cell.addSurface(surface, halfspace)

      # Add the Cell to the global dictionary of all Cells
      self._cells[cell_key] = cell



  def _read_universes(self):

    # Initialize dictionary for each Universe
    # Keys   - Universe keys
    # Values - Universe objects
    self._universes = dict()

    self._n_universes = self._get_int(path='model/universes/n_universes')[0]

    # OpenMC Universe IDs (redefined internally from user definitions)
    self._universe_ids = self._get_int(
         self._n_universes, path='model/universes/ids')

    # User-defined Universe IDs
    self._universe_keys = self._get_int(
         self._n_universes, path='model/universes/keys')

    # Build dictionary of Universes
    base = 'model/universes/universe '

    # Iterate over all Universes
    for universe_key in self._universe_keys:

      # Create this Universe
      universe = openmc.Universe(universe_id=universe_key)

      # Read the number of Cells in this Universe
      n_cells = self._get_int(
           path='{0}{1}/n_cells'.format(base, universe_key))[0]

      # If the Universe contains Cells, add them to the Universe
      if n_cells > 0:

        # Read Cell keys (internal OpenMC IDs)
        cell_keys = self._get_int(
             n_cells, path='{0}{1}/cell_ids'.format(base, universe_key))

        # Add each Cell to the Universe
        for cell_key in cell_keys:

          # Get the user-defined Cell ID and the corresponding Cell
          cell_index = self._cell_ids.index(cell_key)
          cell_id = self._cell_keys[cell_index]
          cell = self._cells[cell_id]

          # Add the Cell to the Universe
          universe.addCell(cell)

      # Add the Universe to the global list of Universes
      self._universes[universe_key] = universe


  def _read_lattices(self):

    # Initialize lattices for each Lattice
    # Keys   - Lattice keys
    # Values - Lattice objects
    self._lattices = dict()

    self._n_lattices = self._get_int(path='model/lattices/n_lattices')[0]

    # Read a list of the IDs for each object in the geometry
    if self._n_lattices > 0:

      # OpenMC Lattice IDs (redefined internally from user definitions)
      self._lattice_ids = self._get_int(
        self._n_lattices, path='model/lattices/ids')

      # User-defined Lattice IDs
      self._lattice_keys = self._get_int(
        self._n_lattices, path='model/lattices/keys')

    else:
      self._lattice_keys = dict()
      self._lattice_ids = dict()

    # Build dictionary of Lattices
    base = 'model/lattices/lattice '

    # Iterate over all Lattices
    for lattice_key in self._lattice_keys:

      # Read the type of Lattice (e.g., rectangular or hexagonal)
      type = self._get_int(path='{0}{1}/type'.format(base, lattice_key))[0]

      # Read the Lattice dimensions
      dimension = self._get_int(
           3, path='{0}{1}/dimension'.format(base, lattice_key))

      # Read the Lattice cell widths
      width = self._get_double(
           3, path='{0}{1}/width'.format(base, lattice_key))

      # Read the lower left coordinates
      lower_left = self._get_double(
           3, path='{0}{1}/lower_left'.format(base, lattice_key))

      # Read the ID of the Material outside the lattice (-1 for void)
      outside = self._get_int(
           path='{0}{1}/outside'.format(base, lattice_key))[0]

      # Read the offset(s)
      offset_size = self._get_int(
           path='{0}{1}/offset_size'.format(base, lattice_key))[0]

      # Read the lattice maps
      maps = self._get_int(path='{0}{1}/maps'.format(base, lattice_key))[0]

      if offset_size > 0:

        path = '{0}{1}/offsets'.format(base, lattice_key)

        if self._hdf5:
          offsets = self._f[path].value
          offsets = np.swapaxes(offsets, 0, 1)
          offsets = np.swapaxes(offsets, 1, 2)
        else:
          offsets = np.array(self._get_int(offset_size, path=path))
          offsets.shape = (dimension[2], dimension[1], dimension[0], maps)

      # Read the Universes filling the Lattice
      path = '{0}{1}/universes'.format(base, lattice_key)

      # Import the Universe IDs and convert to 3D array in (x,y,z) order
      if self._hdf5:
        universe_ids = self._f[path][...]
        universe_ids = np.swapaxes(universe_ids, 0, 1)
        universe_ids = np.swapaxes(universe_ids, 1, 2)
      else:
        num_values = dimension[0] * dimension[1] * dimension[2]
        universe_ids = np.array(self._get_int(num_values, path=path))
        universe_ids.shape = dimension

      # Create the Lattice
      lattice = openmc.Lattice(lattice_id=lattice_key, type=LATTICE_TYPES[type])
      lattice.setDimension(tuple(dimension))
      lattice.setLowerLeft(lower_left)
      lattice.setWidth(width)

      # Build array of Universe pointers for the Lattice
      universes = np.ndarray(tuple(universe_ids.shape),
                             dtype=openmc.Universe)

      for x in range(universe_ids.shape[0]):
        for y in range(universe_ids.shape[1]):
          for z in range(universe_ids.shape[2]):
            universes[x,y,z] = self._universes[universe_ids[x,y,z]]

      lattice.setUniverses(universes)

      # If the Material specified outside the Lattice is not void (-1)
      if outside != -1:
        lattice.setOutside(self._materials[outside])

      if offset_size > 0:
        lattice.setOffsets(offsets)

      # Add the Lattice to the global dictionary of all Lattices
      self._lattices[lattice_key] = lattice


  def _finalize_geometry(self):

    # Initialize Geometry object
    self._geometry = openmc.Geometry()

    # Iterate over all Cells and add fill Materials, Universes and Lattices
    for cell_key in self._cell_fills.keys():

      # Determine the fill type ('normal', 'universe', or 'lattice') and ID
      fill_type = self._cell_fills[cell_key][0]
      fill_id = self._cell_fills[cell_key][1]

      # Retrieve the object corresponding to the fill type and ID
      if fill_type == 'normal':
        if fill_id > 0:
          fill = self._materials[fill_id]
        else:
          fill = 'void'
      elif fill_type == 'fill':
        fill = self._universes[fill_id]
      else:
        fill = self._lattices[fill_id]

      # Set the fill for the Cell
      self._cells[cell_key].setFill(fill)

    # Set the root universe for the Geometry
    root_universe = self._universes[0]
    self._geometry.setRootUniverse(root_universe)


  def _read_meshes(self):

    # Initialize dictionaries for the Meshes
    # Keys   - Mesh IDs
    # Values - Mesh objects
    self._meshes = dict()

    # Read the number of Meshes
    self._n_meshes = self._get_int(path='meshes/n_meshes')[0]

    # Read a list of the IDs for each Mesh
    if self._n_meshes > 0:

      # OpenMC Mesh IDs (redefined internally from user definitions)
      self._mesh_ids = self._get_int(
        self._n_meshes, path='meshes/ids')

      # User-defined Mesh IDs
      self._mesh_keys = self._get_int(
        self._n_meshes, path='meshes/keys')

    else:
      self._mesh_keys = list()
      self._mesh_ids = list()

    # Build dictionary of Meshes
    base = 'meshes/mesh '

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

      mesh.setDimension(dimension)
      mesh.setWidth(width)
      mesh.setLowerLeft(lower_left)
      mesh.setUpperRight(upper_right)

      #FIXME: Set the mesh type to 'rectangular' by default
      mesh.setType('rectangular')

      # Add mesh to the global dictionary of all Meshes
      self._meshes[mesh_id] = mesh



  def _read_tallies(self):

    # Initialize dictionaries for the Tallies
    # Keys   - Tally IDs
    # Values - Tally objects
    self._tallies = dict()

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
      self._tally_keys = list()
      self._tally_ids = list()

    base = 'tallies/tally '

    # Iterate over all Tallies
    for tally_key in self._tally_keys:

      # Read user-specified Tally label (if specified)
      label_size = self._get_int(
           path='{0}{1}/label_size'.format(base, tally_key))[0]

      if label_size > 0:
        label = self._get_string(
          label_size, path='{0}{1}/label'.format(base, tally_key))

      # Remove leading and trailing characters from string label
      label = label.lstrip('[\'')
      label = label.rstrip('\']')

      # Read integer Tally estimator type code (analog or tracklength)
      estimator_type = self._get_int(
           path='{0}{1}/estimator'.format(base, tally_key))[0]

      # Read the Tally size specifications
      n_realizations = self._get_int(
           path='{0}{1}/n_realizations'.format(base, tally_key))[0]

      # Create Tally object and assign basic properties
      tally = openmc.Tally(tally_key, label)
      tally.setEstimator(ESTIMATOR_TYPES[estimator_type])
      tally.setNumRealizations(n_realizations)

      # Read the number of Filters
      n_filters = self._get_int(
           path='{0}{1}/n_filters'.format(base, tally_key))[0]

      subbase = '{0}{1}/filter '.format(base, tally_key)

      # Initialize the stride
      stride = 1

      # Initialize all Filters
      for j in range(1, n_filters+1):

        # Read the integer Filter type code
        filter_type = self._get_int(
             path='{0}{1}/type'.format(subbase, j))[0]

        # Read the Filter offset
        offset = self._get_int(path='{0}{1}/offset'.format(subbase, j))[0]

        n_bins = self._get_int(path='{0}{1}/n_bins'.format(subbase, j))[0]

        if n_bins <= 0:
          msg = 'Unable to create Filter {0} for Tally ID={2} since ' \
                'no bins were specified'.format(j, tally_key)
          raise ValueError(msg)

        # Read the bin values
        if FILTER_TYPES[filter_type] in ['energy', 'energyout']:
          bins = self._get_double(
            n_bins+1, path='{0}{1}/bins'.format(subbase, j))

        elif FILTER_TYPES[filter_type] == 'mesh':
          bins = self._get_int(path='{0}{1}/bins'.format(subbase, j))[0]

        elif FILTER_TYPES[filter_type] == 'surface':

          bins = self._get_int(n_bins, path='{0}{1}/bins'.format(subbase, j))

          # Get the user-defined Surface IDs from the internal OpenMC
          # Surface IDs for the bin from the Surface key mapping
          for k, bin in enumerate(bins):
            surface_index = self._surface_ids.index(bin)
            surface_id = self._surface_keys[surface_index]
            bins[k] = surface_id

        elif FILTER_TYPES[filter_type] == 'distribcell':

          bins = self._get_int(path='{0}{1}/bins'.format(subbase, j))[0]

          # Get the user-defined Cell ID from the internal OpenMC
          # Cell ID for the bin from the cell key mapping
          cell_index = self._cell_ids.index(bins)
          cell_id = self._cell_keys[cell_index]
          bins = cell_id

        elif FILTER_TYPES[filter_type] in ['cell', 'cellborn']:

          bins = self._get_int(n_bins, path='{0}{1}/bins'.format(subbase, j))

          # Get the user-defined Cell IDs from the internal OpenMC
          # Cell IDs for the bin from the cell key mapping
          for k, bin in enumerate(bins):
            cell_index = self._cell_ids.index(bin)
            cell_id = self._cell_keys[cell_index]
            bins[k] = cell_id

        elif FILTER_TYPES[filter_type] == 'universe':

          bins = self._get_int(
            n_bins, path='{0}{1}/bins'.format(subbase, j))

          # Get the user-defined Universe IDs from the internal OpenMC
          # Surface IDs for the bin from the Surface key mapping
          for k, bin in enumerate(bins):
            universe_index = self._universe_ids.index(bin)
            universe_id = self._universe_keys[universe_index]
            bins[k] = universe_id

        elif FILTER_TYPES[filter_type] == 'material':

          bins = self._get_int(n_bins, path='{0}{1}/bins'.format(subbase, j))

          # Get the user-defined Material IDs from the internal OpenMC
          # Surface IDs for the bin from the Surface key mapping
          for k, bin in enumerate(bins):
            material_index = self._material_ids.index(bin)
            material_id = self._material_keys[material_index]
            bins[k] = material_id

        # Create Filter object
        filter = openmc.Filter(FILTER_TYPES[filter_type], bins)
        filter.setOffset(offset)
        filter.setStride(stride)
        filter.setNumBins(n_bins)

        if FILTER_TYPES[filter_type] == 'mesh':
          filter.setMesh(self._meshes[bins])

        # Add Filter to the Tally
        tally.addFilter(filter)

        # Update the stride for the next Filter
        stride *= n_bins

      # Read Nuclide bins
      n_nuclides = self._get_int(
           path='{0}{1}/n_nuclides'.format(base, tally_key))[0]

      nuclide_ids = self._get_int(
           n_nuclides, path='{0}{1}/nuclides'.format(base, tally_key))

      # Add all Nuclides to the Tally
      for nuclide_id in nuclide_ids:
        if nuclide_id == -1:
          tally.addNuclide(openmc.Nuclide('total'))
        else:
          tally.addNuclide(self._nuclides[nuclide_id])

      # Read score bins
      n_score_bins = self._get_int(
           path='{0}{1}/n_score_bins'.format(base, tally_key))[0]

      tally.setNumScoreBins(n_score_bins)

      scores = [SCORE_TYPES[j] for j in self._get_int(
                n_score_bins, path='{0}{1}/score_bins'.format(base, tally_key))]

      # Read the scattering moment order for all scores
      scatt_order = self._get_int(
           n_score_bins, path='{0}{1}/moment_order'.format(base, tally_key))

      # Add the scores to the Tally
      for j, score in enumerate(scores):

        # If this is a scattering moment, insert the scattering order
        if 'scatter-n' in score or 'scatter-yn' in score or 'scatter-pn' in score:
          score = score.replace('n', str(scatt_order[j]))

        tally.addScore(score)

      # Add Tally to the global dictionary of all Tallies
      self._tallies[tally_key] = tally



  def read_results(self):

     # Number of realizations for global Tallies
    self._n_realizations = self._get_int(path='n_realizations')[0]

    #TODO: Make global tallies use the general Tally class
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
        num_tot_bins = tally.getNumBins()

        # Extract Tally data from the file
        if self._hdf5:
          data = self._f['{0}{1}/results'.format(base, tally_key)].value
          sum = data['sum']
          sum_sq = data['sum_sq']
#          sum = np.column_stack(data['sum'])
#          sum_sq = np.column_stack(data['sum_sq'])

        else:
          results = np.array(self._get_double(2*num_tot_bins))
          sum = results[0::2]
          sum_sq = results[1::2]

        # Define a routine to convert 0 to 1
        nonzero = lambda val: 1 if not val else val

        # Reshape the results arrays
        new_shape = (nonzero(tally.getNumFilterBins()),
                     nonzero(tally.getNumNuclides()),
                     nonzero(tally.getNumScoreBins()))

        sum = np.reshape(sum, new_shape)
        sum_sq = np.reshape(sum_sq, new_shape)

        # Set the data for this Tally
        tally.setResults(sum=sum, sum_sq=sum_sq)

    # Indicate that Tally results have been read
    self._results = True


  def read_source(self):

    # Check whether Tally results have been read
    if not self._results:
      self.read_results()

    # Check if source bank is in statepoint
    if not self._source_present:
      print('Unable to read the source since it is not in the statepoint file')
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
    """Computes confidence intervals for each Tally bin."""

    # Determine significance level and percentile for two-sided CI
    alpha = 1 - confidence
    percentile = 1 - alpha/2

    # Calculate t-value
    t_value = scipy.stats.t.ppf(percentile, self._n_realizations - 1)
    self.compute_stdev(t_value)


  def compute_stdev(self, t_value=1.0):
    """
    Computes the sample mean and the standard deviation of the mean
    for each Tally bin.
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


    # Calculate the sample mean and standard deviation for user-defined Tallies
    for tally_id, tally in self._tallies.items():
      tally.computeStdDev(t_value)


  def get_tally(self, score, filters, nuclides,
                label='', estimator='tracklength'):
    """Finds and returns a Tally object with certain properties.

    Parameters
    ----------
    score : str
        The score string

    filters : list
        A list of Filter objects

    nuclides : list
        A list of Nuclide objects

    label : str
        The label specified for the Tally (default is '')

    estimator: str
        The type of estimator ('tracklength' (default) or 'analog')
    """

    # Loop over the domain-to-tallies mapping to find the Tally
    tally = None

    # Iterate over all tallies to find the appropriate one
    for tally_id, test_tally in self._tallies.items():

      # Determine if the queried Tally label is the same as this Tally
      if not label == test_tally._label:
        continue

      # Determine if the queried Tally estimator is the same as this Tally
      if not estimator == test_tally._estimator:
        continue

      # Determine if the queried Tally scores are the same as this Tally
      if not score in test_tally._scores:
        continue

      # Determine if the queried Tally filters is the same length as this Tally
      if len(filters) != len(test_tally._filters):
        continue

      # Determine if the queried Tally filters are the same as this Tally
      contains_filters = True

      # Iterate over the filters requested by the user
      for filter in filters:
        if not filter in test_tally._filters:
          contains_filters = False
          break

      # Determine if the queried Nuclide is in this Tally
      contains_nuclides = True

      # Iterate over the Nuclides requested by the user
      for nuclide in nuclides:
        if not nuclide in test_tally._nuclides:
          contains_nuclides = False
          break

      # If the Tally contained all Filters and Nuclides, return the Tally
      if contains_filters and contains_nuclides:
        tally = test_tally
        break

    # If we did not find the Tally, return an error message
    if tally is None:
      raise LookupError('Unable to get Tally')

    return tally


  def get_tally_id(self, score, filters, label='', estimator='tracklength'):
    """Retrieve the Tally ID for a given list of filters and score(s).

    Parameters
    ----------
    score : str
        The score string

    filters : list
        A list of Filter objects

    label : str
        The label specified for the Tally (default is '')

    estimator: str
        The type of estimator ('tracklength' (default) or 'analog')
    """

    tally = self.get_tally(score, filters, label, estimator)
    return tally._id



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
