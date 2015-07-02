import numpy as np

import openmc


class Summary(object):
    """Information summarizing the geometry, materials, and tallies used in a
    simulation.

    """

    def __init__(self, filename):
        # A user may not have h5py, but they can still use the rest of the
        # Python API so we'll only try to import h5py if the user actually inits
        # a Summary object.
        import h5py

        openmc.reset_auto_ids()

        if not filename.endswith(('.h5', '.hdf5')):
            msg = 'Unable to open "{0}" which is not an HDF5 summary file'
            raise ValueError(msg)

        self._f = h5py.File(filename, 'r')
        self.openmc_geometry = None
        self.opencg_geometry = None

        self._read_metadata()
        self._read_geometry()
        self._read_tallies()

    def _read_metadata(self):
        # Read OpenMC version
        self.version = [self._f['version_major'][0],
                         self._f['version_minor'][0],
                         self._f['version_release'][0]]
        # Read date and time
        self.date_and_time = self._f['date_and_time'][...]

        self.n_batches = self._f['n_batches'][0]
        self.n_particles = self._f['n_particles'][0]
        self.n_active = self._f['n_active'][0]
        self.n_inactive = self._f['n_inactive'][0]
        self.gen_per_batch = self._f['gen_per_batch'][0]
        self.n_procs = self._f['n_procs'][0]

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
        self.n_nuclides = self._f['nuclides/n_nuclides'][0]

        # Initialize dictionary for each Nuclide
        # Keys   - Nuclide ZAIDs
        # Values - Nuclide objects
        self.nuclides = {}

        for key in self._f['nuclides'].keys():
            if key == 'n_nuclides':
                continue

            index = self._f['nuclides'][key]['index'][0]
            alias = self._f['nuclides'][key]['alias'][0]
            zaid = self._f['nuclides'][key]['zaid'][0]

            # Read the Nuclide's name (e.g., 'H-1' or 'U-235')
            name = alias.split('.')[0]

            # Read the Nuclide's cross-section identifier (e.g., '70c')
            xs = alias.split('.')[1]

            # Initialize this Nuclide and add to global dictionary of Nuclides
            if 'nat' in name:
                self.nuclides[zaid] = openmc.Element(name=name, xs=xs)
            else:
                self.nuclides[zaid] = openmc.Nuclide(name=name, xs=xs)
                self.nuclides[zaid].zaid = zaid

    def _read_materials(self):
        self.n_materials = self._f['materials/n_materials'][0]

        # Initialize dictionary for each Material
        # Keys     - Material keys
        # Values   - Material objects
        self.materials = {}

        for key in self._f['materials'].keys():
            if key == 'n_materials':
                continue

            material_id = int(key.lstrip('material '))
            index = self._f['materials'][key]['index'][0]
            name = self._f['materials'][key]['name'][0]
            density = self._f['materials'][key]['atom_density'][0]
            nuc_densities = self._f['materials'][key]['nuclide_densities'][...]
            nuclides = self._f['materials'][key]['nuclides'][...]
            n_sab = self._f['materials'][key]['n_sab'][0]

            sab_names = []
            sab_xs = []

            # Read the names of the S(a,b) tables for this Material
            for i in range(1, n_sab+1):
                sab_table = self._f['materials'][key]['sab_tables'][str(i)][0]

                # Read the cross-section identifiers for each S(a,b) table
                sab_names.append(sab_table.split('.')[0])
                sab_xs.append(sab_table.split('.')[1])

            # Create the Material
            material = openmc.Material(material_id=material_id, name=name)

            # Set the Material's density to g/cm3 - this is what is used in OpenMC
            material.set_density(density=density, units='g/cm3')

            # Add all Nuclides to the Material
            for i, zaid in enumerate(nuclides):
                nuclide = self.get_nuclide_by_zaid(zaid)
                density = nuc_densities[i]

                if isinstance(nuclide, openmc.Nuclide):
                    material.add_nuclide(nuclide, percent=density, percent_type='ao')
                elif isinstance(nuclide, openmc.Element):
                    material.add_element(nuclide, percent=density, percent_type='ao')

            # Add S(a,b) table(s?) to the Material
            for i in range(n_sab):
                name = sab_names[i]
                xs = sab_xs[i]
                material.add_s_alpha_beta(name, xs)

            # Add the Material to the global dictionary of all Materials
            self.materials[index] = material

    def _read_surfaces(self):
        self.n_surfaces = self._f['geometry/n_surfaces'][0]

        # Initialize dictionary for each Surface
        # Keys     - Surface keys
        # Values   - Surfacee objects
        self.surfaces = {}

        for key in self._f['geometry/surfaces'].keys():
            if key == 'n_surfaces':
                continue

            surface_id = int(key.lstrip('surface '))
            index = self._f['geometry/surfaces'][key]['index'][0]
            name = self._f['geometry/surfaces'][key]['name'][0]
            surf_type = self._f['geometry/surfaces'][key]['type'][...][0]
            bc = self._f['geometry/surfaces'][key]['boundary_condition'][...][0]
            coeffs = self._f['geometry/surfaces'][key]['coefficients'][...]

            # Create the Surface based on its type
            if surf_type == 'X Plane':
                x0 = coeffs[0]
                surface = openmc.XPlane(surface_id, bc, x0, name)

            elif surf_type == 'Y Plane':
                y0 = coeffs[0]
                surface = openmc.YPlane(surface_id, bc, y0, name)

            elif surf_type == 'Z Plane':
                z0 = coeffs[0]
                surface = openmc.ZPlane(surface_id, bc, z0, name)

            elif surf_type == 'Plane':
                A = coeffs[0]
                B = coeffs[1]
                C = coeffs[2]
                D = coeffs[3]
                surface = openmc.Plane(surface_id, bc, A, B, C, D, name)

            elif surf_type == 'X Cylinder':
                y0 = coeffs[0]
                z0 = coeffs[1]
                R = coeffs[2]
                surface = openmc.XCylinder(surface_id, bc, y0, z0, R, name)

            elif surf_type == 'Y Cylinder':
                x0 = coeffs[0]
                z0 = coeffs[1]
                R = coeffs[2]
                surface = openmc.YCylinder(surface_id, bc, x0, z0, R, name)

            elif surf_type == 'Z Cylinder':
                x0 = coeffs[0]
                y0 = coeffs[1]
                R = coeffs[2]
                surface = openmc.ZCylinder(surface_id, bc, x0, y0, R, name)

            elif surf_type == 'Sphere':
                x0 = coeffs[0]
                y0 = coeffs[1]
                z0 = coeffs[2]
                R = coeffs[3]
                surface = openmc.Sphere(surface_id, bc, x0, y0, z0, R, name)

            elif surf_type in ['X Cone', 'Y Cone', 'Z Cone']:
                x0 = coeffs[0]
                y0 = coeffs[1]
                z0 = coeffs[2]
                R2 = coeffs[3]

                if surf_type == 'X Cone':
                    surface = openmc.XCone(surface_id, bc, x0, y0, z0, R2, name)
                if surf_type == 'Y Cone':
                    surface = openmc.YCone(surface_id, bc, x0, y0, z0, R2, name)
                if surf_type == 'Z Cone':
                    surface = openmc.ZCone(surface_id, bc, x0, y0, z0, R2, name)

            # Add Surface to global dictionary of all Surfaces
            self.surfaces[index] = surface

    def _read_cells(self):
        self.n_cells = self._f['geometry/n_cells'][0]

        # Initialize dictionary for each Cell
        # Keys     - Cell keys
        # Values   - Cell objects
        self.cells = {}

        # Initialize dictionary for each Cell's fill
        # (e.g., Material, Universe or Lattice ID)
        # This dictionary is used later to link the fills with
        # the corresponding objects
        # Keys     - Cell keys
        # Values   - Filling Material, Universe or Lattice ID
        self._cell_fills = {}

        for key in self._f['geometry/cells'].keys():
            if key == 'n_cells':
                continue

            cell_id = int(key.lstrip('cell '))
            index = self._f['geometry/cells'][key]['index'][0]
            name = self._f['geometry/cells'][key]['name'][0]
            fill_type = self._f['geometry/cells'][key]['fill_type'][...][0]

            if fill_type == 'normal':
                fill = self._f['geometry/cells'][key]['material'][0]
            elif fill_type == 'universe':
                fill = self._f['geometry/cells'][key]['fill'][0]
            else:
                fill = self._f['geometry/cells'][key]['lattice'][0]

            if 'surfaces' in self._f['geometry/cells'][key].keys():
                surfaces = self._f['geometry/cells'][key]['surfaces'][...]
            else:
                surfaces = []

            # Create this Cell
            cell = openmc.Cell(cell_id=cell_id, name=name)

            if fill_type == 'universe':
                maps = self._f['geometry/cells'][key]['maps'][0]

                if maps > 0:
                    offset = self._f['geometry/cells'][key]['offset'][...]
                    cell.set_offset(offset)

                translated = self._f['geometry/cells'][key]['translated'][0]
                if translated:
                    translation = \
                      self._f['geometry/cells'][key]['translation'][...]
                    translation = np.asarray(translation, dtype=np.float64)
                    cell.translation = translation

                rotated = self._f['geometry/cells'][key]['rotated'][0]
                if rotated:
                    rotation = \
                      self._f['geometry/cells'][key]['rotation'][...]
                    rotation = np.asarray(rotation, dtype=np.int)
                    cell.rotation = rotation

            # Store Cell fill information for after Universe/Lattice creation
            self._cell_fills[index] = (fill_type, fill)

            # Iterate over all Surfaces and add them to the Cell
            for surface_halfspace in surfaces:

                halfspace = np.sign(surface_halfspace)
                surface_id = np.abs(surface_halfspace)
                surface = self.surfaces[surface_id]
                cell.add_surface(surface, halfspace)

            # Add the Cell to the global dictionary of all Cells
            self.cells[index] = cell

    def _read_universes(self):
        self.n_universes = self._f['geometry/n_universes'][0]

        # Initialize dictionary for each Universe
        # Keys     - Universe keys
        # Values   - Universe objects
        self.universes = {}

        for key in self._f['geometry/universes'].keys():
            if key == 'n_universes':
                continue

            universe_id = int(key.lstrip('universe '))
            index = self._f['geometry/universes'][key]['index'][0]
            cells = self._f['geometry/universes'][key]['cells'][...]

            # Create this Universe
            universe = openmc.Universe(universe_id=universe_id)

            # Add each Cell to the Universe
            for cell_id in cells:
                cell = self.cells[cell_id]
                universe.add_cell(cell)

            # Add the Universe to the global list of Universes
            self.universes[index] = universe

    def _read_lattices(self):
        self.n_lattices = self._f['geometry/n_lattices'][0]

        # Initialize lattices for each Lattice
        # Keys     - Lattice keys
        # Values   - Lattice objects
        self.lattices = {}

        for key in self._f['geometry/lattices'].keys():
            if key == 'n_lattices':
                continue

            lattice_id = int(key.lstrip('lattice '))
            index = self._f['geometry/lattices'][key]['index'][0]
            name = self._f['geometry/lattices'][key]['name'][0]
            lattice_type = self._f['geometry/lattices'][key]['type'][...][0]
            maps = self._f['geometry/lattices'][key]['maps'][0]
            offset_size = self._f['geometry/lattices'][key]['offset_size'][0]

            if offset_size > 0:
                offsets = self._f['geometry/lattices'][key]['offsets'][...]

            if lattice_type == 'rectangular':
                dimension = self._f['geometry/lattices'][key]['dimension'][...]
                lower_left = \
                     self._f['geometry/lattices'][key]['lower_left'][...]
                pitch = self._f['geometry/lattices'][key]['pitch'][...]
                outer = self._f['geometry/lattices'][key]['outer'][0]

                universe_ids = \
                     self._f['geometry/lattices'][key]['universes'][...]
                universe_ids = np.swapaxes(universe_ids, 0, 1)
                universe_ids = np.swapaxes(universe_ids, 1, 2)

                # Create the Lattice
                lattice = openmc.RectLattice(lattice_id=lattice_id, name=name)
                lattice.dimension = tuple(dimension)
                lattice.lower_left = lower_left
                lattice.pitch = pitch

                # If the Universe specified outer the Lattice is not void (-22)
                if outer != -22:
                    lattice.outer = self.universes[outer]

                # Build array of Universe pointers for the Lattice
                universes = \
                    np.ndarray(tuple(universe_ids.shape), dtype=openmc.Universe)

                for x in range(universe_ids.shape[0]):
                    for y in range(universe_ids.shape[1]):
                        for z in range(universe_ids.shape[2]):
                            universes[x, y, z] = \
                                 self.get_universe_by_id(universe_ids[x, y, z])

                # Transpose, reverse y-dimension for appropriate ordering
                shape = universes.shape
                universes = np.transpose(universes, (1, 0, 2))
                universes.shape = shape
                universes = universes[:, ::-1, :]
                lattice.universes = universes

                if offset_size > 0:
                    offsets = np.swapaxes(offsets, 0, 1)
                    offsets = np.swapaxes(offsets, 1, 2)
                    lattice.offsets = offsets

                # Add the Lattice to the global dictionary of all Lattices
                self.lattices[index] = lattice

            if lattice_type == 'hexagonal':
                n_rings = self._f['geometry/lattices'][key]['n_rings'][0]
                n_axial = self._f['geometry/lattices'][key]['n_axial'][0]
                center = self._f['geometry/lattices'][key]['center'][...]
                pitch = self._f['geometry/lattices'][key]['pitch'][...]
                outer = self._f['geometry/lattices'][key]['outer'][0]

                universe_ids = self._f[
                    'geometry/lattices'][key]['universes'][...]

                # Create the Lattice
                lattice = openmc.HexLattice(lattice_id=lattice_id, name=name)
                lattice.num_rings = n_rings
                lattice.num_axial = n_axial
                lattice.center = center
                lattice.pitch = pitch

                # If the Universe specified outer the Lattice is not void (-22)
                if outer != -22:
                    lattice.outer = self.universes[outer]

                # Build array of Universe pointers for the Lattice
                universes = \
                    np.ndarray(tuple(universe_ids.shape), dtype=openmc.Universe)

                for i in range(universe_ids.shape[0]):
                    for j in range(universe_ids.shape[1]):
                        for k in range(universe_ids.shape[2]):
                            if universe_ids[i, j, k] != -1:
                                universes[i, j, k] = self.get_universe_by_id(
                                    universe_ids[i, j, k])

                if offset_size > 0:
                    lattice.offsets = offsets

                # Add the Lattice to the global dictionary of all Lattices
                self.lattices[index] = lattice

    def _finalize_geometry(self):
        # Initialize Geometry object
        self.openmc_geometry = openmc.Geometry()

        # Iterate over all Cells and add fill Materials, Universes and Lattices
        for cell_key in self._cell_fills.keys():
            # Determine fill type ('normal', 'universe', or 'lattice') and ID
            fill_type = self._cell_fills[cell_key][0]
            fill_id = self._cell_fills[cell_key][1]

            # Retrieve the object corresponding to the fill type and ID
            if fill_type == 'normal':
                if fill_id > 0:
                    fill = self.get_material_by_id(fill_id)
                else:
                    fill = 'void'
            elif fill_type == 'universe':
                fill = self.get_universe_by_id(fill_id)
            else:
                fill = self.get_lattice_by_id(fill_id)

            # Set the fill for the Cell
            self.cells[cell_key].fill = fill

        # Set the root universe for the Geometry
        root_universe = self.get_universe_by_id(0)
        self.openmc_geometry.root_universe = root_universe

    def _read_tallies(self):
        # Initialize dictionaries for the Tallies
        # Keys     - Tally IDs
        # Values   - Tally objects
        self.tallies = {}

        # Read the number of tallies
        if 'tallies' not in self._f:
            self.n_tallies = 0
            return

        self.n_tallies = self._f['tallies/n_tallies'][0]

        # OpenMC Tally keys
        all_keys = self._f['tallies/'].keys()
        tally_keys = [key for key in all_keys if 'tally' in key]

        base = 'tallies/tally '

        # Iterate over all Tallies
        for tally_key in tally_keys:

            tally_id = int(tally_key.strip('tally '))
            subbase = '{0}{1}'.format(base, tally_id)

            # Read Tally name metadata
            name_size = self._f['{0}/name_size'.format(subbase)][0]
            if (name_size > 0):
                tally_name = self._f['{0}/name'.format(subbase)][0]
                tally_name = tally_name.lstrip('[\'')
                tally_name = tally_name.rstrip('\']')
            else:
                tally_name = ''

            # Create Tally object and assign basic properties
            tally = openmc.Tally(tally_id, tally_name)

            # Read score metadata
            score_bins = self._f['{0}/score_bins'.format(subbase)][...]
            for score_bin in score_bins:
                tally.add_score(openmc.SCORE_TYPES[score_bin])
            num_score_bins = self._f['{0}/n_score_bins'.format(subbase)][...]
            tally.num_score_bins = num_score_bins

            # Read filter metadata
            num_filters = self._f['{0}/n_filters'.format(subbase)][0]

            # Initialize all Filters
            for j in range(1, num_filters+1):

                subsubbase = '{0}/filter {1}'.format(subbase, j)

                # Read filter type (e.g., "cell", "energy", etc.)
                filter_type_code = self._f['{0}/type'.format(subsubbase)][0]
                filter_type = openmc.FILTER_TYPES[filter_type_code]

                # Read the filter bins
                num_bins = self._f['{0}/n_bins'.format(subsubbase)][0]
                bins = self._f['{0}/bins'.format(subsubbase)][...]

                # Create Filter object
                filter = openmc.Filter(filter_type, bins)
                filter.num_bins = num_bins

                # Add Filter to the Tally
                tally.add_filter(filter)

            # Add Tally to the global dictionary of all Tallies
            self.tallies[tally_id] = tally

    def make_opencg_geometry(self):
        """Create OpenCG geometry based on the information contained in the summary
        file. The geometry is stored as the 'opencg_geometry' attribute.

        """

        try:
            from openmc.opencg_compatible import get_opencg_geometry
        except ImportError:
            msg = 'Unable to import opencg which is needed ' \
                  'by Summary.make_opencg_geometry()'
            raise ImportError(msg)

        if self.opencg_geometry is None:
            self.opencg_geometry = get_opencg_geometry(self.openmc_geometry)

    def get_nuclide_by_zaid(self, zaid):
        """Return a Nuclide object given the 'zaid' identifier for the nuclide.

        Parameters
        ----------
        zaid : int
            1000*Z + A, where Z is the atomic number of the nuclide and A is the
            mass number. For example, the zaid for U-235 is 92235.

        Returns
        -------
        nuclide : openmc.nuclide.Nuclide or None
            Nuclide matching the specified zaid, or None if no matching object
            is found.

        """

        for index, nuclide in self.nuclides.items():
            if nuclide._zaid == zaid:
                return nuclide

        return None

    def get_material_by_id(self, material_id):
        """Return a Material object given the material id

        Parameters
        ----------
        id : int
            Unique identifier for the material

        Returns
        -------
        material : openmc.material.Material
            Material with given id

        """

        for index, material in self.materials.items():
            if material._id == material_id:
                return material

        return None

    def get_surface_by_id(self, surface_id):
        """Return a Surface object given the surface id

        Parameters
        ----------
        id : int
            Unique identifier for the surface

        Returns
        -------
        surface : openmc.surface.Surface
            Surface with given id

        """

        for index, surface in self.surfaces.items():
            if surface._id == surface_id:
                return surface

        return None

    def get_cell_by_id(self, cell_id):
        """Return a Cell object given the cell id

        Parameters
        ----------
        id : int
            Unique identifier for the cell

        Returns
        -------
        cell : openmc.universe.Cell
            Cell with given id

        """

        for index, cell in self.cells.items():
            if cell._id == cell_id:
                return cell

        return None

    def get_universe_by_id(self, universe_id):
        """Return a Universe object given the universe id

        Parameters
        ----------
        id : int
            Unique identifier for the universe

        Returns
        -------
        universe : openmc.universe.Universe
            Universe with given id

        """

        for index, universe in self.universes.items():
            if universe._id == universe_id:
                return universe

        return None

    def get_lattice_by_id(self, lattice_id):
        """Return a Lattice object given the lattice id

        Parameters
        ----------
        id : int
            Unique identifier for the lattice

        Returns
        -------
        lattice : openmc.universe.Lattice
            Lattice with given id

        """

        for index, lattice in self.lattices.items():
            if lattice._id == lattice_id:
                return lattice

        return None
