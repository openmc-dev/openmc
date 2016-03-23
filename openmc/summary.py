from collections import Iterable
import numpy as np
import re

import openmc
from openmc.region import Region


class Summary(object):
    """Information summarizing the geometry, materials, and tallies used in a
    simulation.

    Attributes
    ----------
    openmc_geometry : openmc.Geometry
        An OpenMC geometry object reconstructed from the summary file
    opencg_geometry : opencg.Geometry
        An OpenCG geometry object equivalent to the OpenMC geometry
        encapsulated by the summary file. Use of this attribute requires
        installation of the OpenCG Python module.

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
        self._openmc_geometry = None
        self._opencg_geometry = None

        self._read_metadata()
        self._read_geometry()
        self._read_tallies()

    @property
    def openmc_geometry(self):
        return self._openmc_geometry

    @property
    def opencg_geometry(self):
        if self._opencg_geometry is None:
            from openmc.opencg_compatible import get_opencg_geometry
            self._opencg_geometry = get_opencg_geometry(self.openmc_geometry)
        return self._opencg_geometry

    def _read_metadata(self):
        # Read OpenMC version
        self.version = [self._f['version_major'].value,
                         self._f['version_minor'].value,
                         self._f['version_release'].value]
        # Read date and time
        self.date_and_time = self._f['date_and_time'][...]

        # Read if continuous-energy or multi-group
        self.run_CE = (self._f['run_CE'].value == 1)

        self.n_batches = self._f['n_batches'].value
        self.n_particles = self._f['n_particles'].value
        self.n_active = self._f['n_active'].value
        self.n_inactive = self._f['n_inactive'].value
        self.gen_per_batch = self._f['gen_per_batch'].value
        self.n_procs = self._f['n_procs'].value

    def _read_geometry(self):
        # Read in and initialize the Materials and Geometry
        self._read_materials()
        self._read_surfaces()
        self._read_cells()
        self._read_universes()
        self._read_lattices()
        self._finalize_geometry()

    def _read_materials(self):
        self.n_materials = self._f['n_materials'].value

        # Initialize dictionary for each Material
        # Keys     - Material keys
        # Values   - Material objects
        self.materials = {}

        for key in self._f['materials'].keys():
            if key == 'n_materials':
                continue

            material_id = int(key.lstrip('material '))
            index = self._f['materials'][key]['index'].value
            name = self._f['materials'][key]['name'].value.decode()
            density = self._f['materials'][key]['atom_density'].value
            nuc_densities = self._f['materials'][key]['nuclide_densities'][...]
            nuclides = self._f['materials'][key]['nuclides'].value

            # Create the Material
            material = openmc.Material(material_id=material_id, name=name)

            # Read the names of the S(a,b) tables for this Material and add them
            if 'sab_names' in self._f['materials'][key]:
                sab_tables = self._f['materials'][key]['sab_names'].value
                for sab_table in sab_tables:
                    name, xs = sab_table.decode().split('.')
                    material.add_s_alpha_beta(name, xs)

            # Set the Material's density to atom/b-cm as used by OpenMC
            material.set_density(density=density, units='atom/b-cm')

            # Add all nuclides to the Material
            for fullname, density in zip(nuclides, nuc_densities):
                fullname = fullname.decode().strip()
                name, xs = fullname.split('.')

                if 'nat' in name:
                    material.add_element(openmc.Element(name=name, xs=xs),
                                         percent=density, percent_type='ao')
                else:
                    material.add_nuclide(openmc.Nuclide(name=name, xs=xs),
                                         percent=density, percent_type='ao')

            # Add the Material to the global dictionary of all Materials
            self.materials[index] = material

    def _read_surfaces(self):
        self.n_surfaces = self._f['geometry/n_surfaces'].value

        # Initialize dictionary for each Surface
        # Keys     - Surface keys
        # Values   - Surfacee objects
        self.surfaces = {}

        for key in self._f['geometry/surfaces'].keys():
            if key == 'n_surfaces':
                continue

            surface_id = int(key.lstrip('surface '))
            index = self._f['geometry/surfaces'][key]['index'].value
            name = self._f['geometry/surfaces'][key]['name'].value.decode()
            surf_type = self._f['geometry/surfaces'][key]['type'].value.decode()
            bc = self._f['geometry/surfaces'][key]['boundary_condition'].value.decode()
            coeffs = self._f['geometry/surfaces'][key]['coefficients'][...]

            # Create the Surface based on its type
            if surf_type == 'x-plane':
                x0 = coeffs[0]
                surface = openmc.XPlane(surface_id, bc, x0, name)

            elif surf_type == 'y-plane':
                y0 = coeffs[0]
                surface = openmc.YPlane(surface_id, bc, y0, name)

            elif surf_type == 'z-plane':
                z0 = coeffs[0]
                surface = openmc.ZPlane(surface_id, bc, z0, name)

            elif surf_type == 'plane':
                A = coeffs[0]
                B = coeffs[1]
                C = coeffs[2]
                D = coeffs[3]
                surface = openmc.Plane(surface_id, bc, A, B, C, D, name)

            elif surf_type == 'x-cylinder':
                y0 = coeffs[0]
                z0 = coeffs[1]
                R = coeffs[2]
                surface = openmc.XCylinder(surface_id, bc, y0, z0, R, name)

            elif surf_type == 'y-cylinder':
                x0 = coeffs[0]
                z0 = coeffs[1]
                R = coeffs[2]
                surface = openmc.YCylinder(surface_id, bc, x0, z0, R, name)

            elif surf_type == 'z-cylinder':
                x0 = coeffs[0]
                y0 = coeffs[1]
                R = coeffs[2]
                surface = openmc.ZCylinder(surface_id, bc, x0, y0, R, name)

            elif surf_type == 'sphere':
                x0 = coeffs[0]
                y0 = coeffs[1]
                z0 = coeffs[2]
                R = coeffs[3]
                surface = openmc.Sphere(surface_id, bc, x0, y0, z0, R, name)

            elif surf_type in ['x-cone', 'y-cone', 'z-cone']:
                x0 = coeffs[0]
                y0 = coeffs[1]
                z0 = coeffs[2]
                R2 = coeffs[3]

                if surf_type == 'x-cone':
                    surface = openmc.XCone(surface_id, bc, x0, y0, z0, R2, name)
                if surf_type == 'y-cone':
                    surface = openmc.YCone(surface_id, bc, x0, y0, z0, R2, name)
                if surf_type == 'z-cone':
                    surface = openmc.ZCone(surface_id, bc, x0, y0, z0, R2, name)

            elif surf_type == 'quadric':
                a, b, c, d, e, f, g, h, j, k = coeffs
                surface = openmc.Quadric(surface_id, bc, a, b, c, d, e, f,
                                         g, h, j, k, name)

            # Add Surface to global dictionary of all Surfaces
            self.surfaces[index] = surface

    def _read_cells(self):
        self.n_cells = self._f['geometry/n_cells'].value

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
            index = self._f['geometry/cells'][key]['index'].value
            name = self._f['geometry/cells'][key]['name'].value.decode()
            fill_type = self._f['geometry/cells'][key]['fill_type'].value.decode()

            if fill_type == 'normal':
                fill = self._f['geometry/cells'][key]['material'].value
            elif fill_type == 'universe':
                fill = self._f['geometry/cells'][key]['fill'].value
            else:
                fill = self._f['geometry/cells'][key]['lattice'].value

            if 'region' in self._f['geometry/cells'][key].keys():
                region = self._f['geometry/cells'][key]['region'].value.decode()
            else:
                region = []

            # Create this Cell
            cell = openmc.Cell(cell_id=cell_id, name=name)

            if fill_type == 'universe':
                if 'offset' in self._f['geometry/cells'][key]:
                    offset = self._f['geometry/cells'][key]['offset'][...]
                    cell.offsets = offset

                if 'translation' in self._f['geometry/cells'][key]:
                    translation = \
                      self._f['geometry/cells'][key]['translation'][...]
                    translation = np.asarray(translation, dtype=np.float64)
                    cell.translation = translation

                if 'rotation' in self._f['geometry/cells'][key]:
                    rotation = \
                      self._f['geometry/cells'][key]['rotation'][...]
                    rotation = np.asarray(rotation, dtype=np.int)
                    cell.rotation = rotation

            # Store Cell fill information for after Universe/Lattice creation
            self._cell_fills[index] = (fill_type, fill)

            # Generate Region object given infix expression
            if region:
                cell.region = Region.from_expression(
                    region, {s.id: s for s in self.surfaces.values()})

            # Get the distribcell index
            ind = self._f['geometry/cells'][key]['distribcell_index'].value
            if ind != 0:
               cell.distribcell_index = ind

            # Add the Cell to the global dictionary of all Cells
            self.cells[index] = cell

    def _read_universes(self):
        self.n_universes = self._f['geometry/n_universes'].value

        # Initialize dictionary for each Universe
        # Keys     - Universe keys
        # Values   - Universe objects
        self.universes = {}

        for key in self._f['geometry/universes'].keys():
            if key == 'n_universes':
                continue

            universe_id = int(key.lstrip('universe '))
            index = self._f['geometry/universes'][key]['index'].value
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
        self.n_lattices = self._f['geometry/n_lattices'].value

        # Initialize lattices for each Lattice
        # Keys     - Lattice keys
        # Values   - Lattice objects
        self.lattices = {}

        for key in self._f['geometry/lattices'].keys():
            if key == 'n_lattices':
                continue

            lattice_id = int(key.lstrip('lattice '))
            index = self._f['geometry/lattices'][key]['index'].value
            name = self._f['geometry/lattices'][key]['name'].value.decode()
            lattice_type = self._f['geometry/lattices'][key]['type'].value.decode()

            if 'offsets' in self._f['geometry/lattices'][key]:
                offsets = self._f['geometry/lattices'][key]['offsets'][...]
            else:
                offsets = None

            if lattice_type == 'rectangular':
                dimension = self._f['geometry/lattices'][key]['dimension'][...]
                lower_left = \
                     self._f['geometry/lattices'][key]['lower_left'][...]
                pitch = self._f['geometry/lattices'][key]['pitch'][...]
                outer = self._f['geometry/lattices'][key]['outer'].value
                universe_ids = \
                    self._f['geometry/lattices'][key]['universes'][...]

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

                for z in range(universe_ids.shape[0]):
                    for y in range(universe_ids.shape[1]):
                        for x in range(universe_ids.shape[2]):
                            universes[z, y, x] = \
                                 self.get_universe_by_id(universe_ids[z, y, x])

                # Use 2D NumPy array to store lattice universes for 2D lattices
                if len(dimension) == 2:
                    universes = np.squeeze(universes)
                    universes = np.atleast_2d(universes)

                # Set the universes for the lattice
                lattice.universes = universes

                # Set the distribcell offsets for the lattice
                if offsets is not None:
                    lattice.offsets = offsets[:, ::-1, :]

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

                # Build array of Universe pointers for the Lattice.  Note that
                # we need to convert between the HDF5's square array of
                # (x, alpha, z) to the Python API's format of a ragged nested
                # list of (z, ring, theta).
                universes = []
                for z in range(lattice.num_axial):
                    # Add a list for this axial level.
                    universes.append([])
                    x = lattice.num_rings - 1
                    a = 2*lattice.num_rings - 2
                    for r in range(lattice.num_rings - 1, 0, -1):
                        # Add a list for this ring.
                        universes[-1].append([])

                        # Climb down the top-right.
                        for i in range(r):
                            universes[-1][-1].append(universe_ids[z, a, x])
                            x += 1
                            a -= 1

                        # Climb down the right.
                        for i in range(r):
                            universes[-1][-1].append(universe_ids[z, a, x])
                            a -= 1

                        # Climb down the bottom-right.
                        for i in range(r):
                            universes[-1][-1].append(universe_ids[z, a, x])
                            x -= 1

                        # Climb up the bottom-left.
                        for i in range(r):
                            universes[-1][-1].append(universe_ids[z, a, x])
                            x -= 1
                            a += 1

                        # Climb up the left.
                        for i in range(r):
                            universes[-1][-1].append(universe_ids[z, a, x])
                            a += 1

                        # Climb up the top-left.
                        for i in range(r):
                            universes[-1][-1].append(universe_ids[z, a, x])
                            x += 1

                        # Move down to the next ring.
                        a -= 1

                        # Convert the ids into Universe objects.
                        universes[-1][-1] = [self.get_universe_by_id(u_id)
                                             for u_id in universes[-1][-1]]

                    # Handle the degenerate center ring separately.
                    u_id = universe_ids[z, a, x]
                    universes[-1].append([self.get_universe_by_id(u_id)])

                # Add the universes to the lattice.
                if len(pitch) == 2:
                    # Lattice is 3D
                    lattice.universes = universes
                else:
                    # Lattice is 2D; extract the only axial level
                    lattice.universes = universes[0]

                if offsets is not None:
                    lattice.offsets = offsets

                # Add the Lattice to the global dictionary of all Lattices
                self.lattices[index] = lattice

    def _finalize_geometry(self):
        # Initialize Geometry object
        self._openmc_geometry = openmc.Geometry()

        # Iterate over all Cells and add fill Materials, Universes and Lattices
        for cell_key in self._cell_fills.keys():
            # Determine fill type ('normal', 'universe', or 'lattice') and ID
            fill_type = self._cell_fills[cell_key][0]
            fill_id = self._cell_fills[cell_key][1]

            # Retrieve the object corresponding to the fill type and ID
            if fill_type == 'normal':
                if isinstance(fill_id, Iterable):
                    fill = [self.get_material_by_id(mat) if mat > 0 else 'void'
                            for mat in fill_id]
                else:
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

        self.n_tallies = self._f['tallies/n_tallies'].value

        # OpenMC Tally keys
        all_keys = self._f['tallies/'].keys()
        tally_keys = [key for key in all_keys if 'tally' in key]

        base = 'tallies/tally '

        # Iterate over all Tallies
        for tally_key in tally_keys:
            tally_id = int(tally_key.strip('tally '))
            subbase = '{0}{1}'.format(base, tally_id)

            # Read Tally name metadata
            tally_name = self._f['{0}/name'.format(subbase)].value.decode()

            # Create Tally object and assign basic properties
            tally = openmc.Tally(tally_id, tally_name)

            # Read scattering moment order strings (e.g., P3, Y1,2, etc.)
            moments = self._f['{0}/moment_orders'.format(subbase)].value

            # Read score metadata
            scores = self._f['{0}/score_bins'.format(subbase)].value
            for j, score in enumerate(scores):
                score = score.decode()

                # If this is a moment, use generic moment order
                pattern = r'-n$|-pn$|-yn$'
                score = re.sub(pattern, '-' + moments[j].decode(), score)
                tally.scores.append(score)

            # Read filter metadata
            num_filters = self._f['{0}/n_filters'.format(subbase)].value

            # Initialize all Filters
            for j in range(1, num_filters+1):
                subsubbase = '{0}/filter {1}'.format(subbase, j)

                # Read filter type (e.g., "cell", "energy", etc.)
                filter_type = self._f['{0}/type'.format(subsubbase)].value.decode()

                # Read the filter bins
                num_bins = self._f['{0}/n_bins'.format(subsubbase)].value
                bins = self._f['{0}/bins'.format(subsubbase)][...]

                # Create Filter object
                new_filter = openmc.Filter(filter_type, bins)
                new_filter.num_bins = num_bins

                # Read in distribcell paths
                if filter_type == 'distribcell':
                    paths = self._f['{0}/paths'.format(subsubbase)][...]
                    paths = [str(path.decode()) for path in paths]
                    new_filter.distribcell_paths = paths

                # Add Filter to the Tally
                tally.filters.append(new_filter)

            # Add Tally to the global dictionary of all Tallies
            self.tallies[tally_id] = tally

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
            if material.id == material_id:
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
            if surface.id == surface_id:
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
            if cell.id == cell_id:
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
            if universe.id == universe_id:
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
            if lattice.id == lattice_id:
                return lattice

        return None
