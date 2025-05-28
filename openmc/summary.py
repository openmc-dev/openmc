from collections.abc import Iterable
import warnings

import h5py
import numpy as np

import openmc
import openmc.checkvalue as cv
from .region import Region

_VERSION_SUMMARY = 6


class Summary:
    """Summary of model used in a simulation.

    Parameters
    ----------
    filename : str or path-like
        Path to file to load

    Attributes
    ----------
    date_and_time : str
        Date and time when simulation began
    geometry : openmc.Geometry
        The geometry reconstructed from the summary file
    materials : openmc.Materials
        The materials reconstructed from the summary file
    nuclides : dict
        Dictionary whose keys are nuclide names and values are atomic weight
        ratios.
    macroscopics : list
        Names of macroscopic data sets
    version: tuple of int
        Version of OpenMC

    """

    def __init__(self, filename):
        filename = str(filename)
        if not filename.endswith(('.h5', '.hdf5')):
            msg = f'Unable to open "{filename}" which is not an HDF5 summary file'
            raise ValueError(msg)

        self._f = h5py.File(filename, 'r')
        cv.check_filetype_version(self._f, 'summary', _VERSION_SUMMARY)

        self._geometry = openmc.Geometry()

        self._fast_materials = {}
        self._fast_surfaces = {}
        self._fast_cells = {}
        self._fast_universes = {}
        self._fast_lattices = {}

        self._materials = openmc.Materials()
        self._nuclides = {}
        self._macroscopics = []

        self._read_nuclides()
        self._read_macroscopics()
        self._read_geometry()

    @property
    def date_and_time(self):
        return self._f.attrs['date_and_time'].decode()

    @property
    def geometry(self):
        return self._geometry

    @property
    def materials(self):
        return self._materials

    @property
    def nuclides(self):
        return self._nuclides

    @property
    def macroscopics(self):
        return self._macroscopics

    @property
    def version(self):
        return tuple(self._f.attrs['openmc_version'])

    def _read_nuclides(self):
        if 'nuclides/names' in self._f:
            names = self._f['nuclides/names'][()]
            awrs = self._f['nuclides/awrs'][()]
            for name, awr in zip(names, awrs):
                self._nuclides[name.decode()] = awr

    def _read_macroscopics(self):
        if 'macroscopics/names' in self._f:
            names = self._f['macroscopics/names'][()]
            self._macroscopics = [name.decode() for name in names]

    def _read_geometry(self):
        with warnings.catch_warnings():
            # We expect that new objects will be created with the same IDs as
            # objects that might already exist in the Python process (if it was
            # also used to create the model), so silence ID warnings
            warnings.simplefilter("ignore", openmc.IDWarning)

            # Read in and initialize the Materials
            self._read_materials()

            # Read native geometry only
            if "dagmc" not in self._f['geometry'].attrs.keys():
                self._read_surfaces()
                cell_fills = self._read_cells()
                self._read_universes()
                self._read_lattices()
                self._finalize_geometry(cell_fills)

    def _read_materials(self):
        for group in self._f['materials'].values():
            material = openmc.Material.from_hdf5(group)

            # Add the material to the Materials collection
            self.materials.append(material)

            # Store in the dictionary of materials for fast queries
            self._fast_materials[material.id] = material

    def _read_surfaces(self):
        periodic_surface_ids = set()
        for group in self._f['geometry/surfaces'].values():
            surface = openmc.Surface.from_hdf5(group)
            # surface may be None for DAGMC surfaces
            if surface:
                self._fast_surfaces[surface.id] = surface
                if surface.boundary_type == "periodic":
                    periodic_surface_ids.add(surface.id)
                
        # assign periodic surfaces
        for surface_id in periodic_surface_ids:
            group = self._f[f'geometry/surfaces/surface {surface_id}']         
            surface = self._fast_surfaces[surface_id]
            periodic_surface_id = int(group['periodic_surface_id'][()].decode())
            surface.periodic_surface = self._fast_surfaces[periodic_surface_id]
                

    def _read_cells(self):

        # Initialize dictionary for each Cell's fill
        cell_fills = {}

        for key, group in self._f['geometry/cells'].items():
            cell_id = int(key.lstrip('cell '))
            name = group['name'][()].decode() if 'name' in group else ''
            fill_type = group['fill_type'][()].decode()

            if fill_type == 'material':
                fill_id = group['material'][()]
            elif fill_type == 'universe':
                fill_id = group['fill'][()]
            else:
                fill_id = group['lattice'][()]

            region = group['region'][()].decode() if 'region' in group else ''

            # Create this Cell
            cell = openmc.Cell(cell_id=cell_id, name=name)

            if fill_type == 'universe':
                if 'translation' in group:
                    translation = group['translation'][()]
                    translation = np.asarray(translation, dtype=np.float64)
                    cell.translation = translation

                if 'rotation' in group:
                    rotation = group['rotation'][()]
                    if rotation.size == 9:
                        rotation.shape = (3, 3)
                    cell.rotation = rotation

            elif fill_type == 'material':
                cell.temperature = group['temperature'][()]

            # Store Cell fill information for after Universe/Lattice creation
            cell_fills[cell.id] = (fill_type, fill_id)

            # Generate Region object given infix expression
            if region:
                cell.region = Region.from_expression(region, self._fast_surfaces)

            # Add the Cell to the global dictionary of all Cells
            self._fast_cells[cell.id] = cell

        return cell_fills

    def _read_universes(self):
        for group in self._f['geometry/universes'].values():
            geom_type = group.get('geom_type')
            if geom_type and geom_type[()].decode() == 'dagmc':
                universe = openmc.DAGMCUniverse.from_hdf5(group)
            else:
                universe = openmc.Universe.from_hdf5(group, self._fast_cells)
            self._fast_universes[universe.id] = universe

    def _read_lattices(self):
        for group in self._f['geometry/lattices'].values():
            lattice = openmc.Lattice.from_hdf5(group, self._fast_universes)
            self._fast_lattices[lattice.id] = lattice

    def _finalize_geometry(self, cell_fills):

        # Keep track of universes that are used as fills. That way, we can
        # determine which universe is NOT used as a fill (and hence is the root
        # universe)
        fill_univ_ids = set()

        # Iterate over all Cells and add fill Materials, Universes and Lattices
        for cell_id, (fill_type, fill_id) in cell_fills.items():
            # Retrieve the object corresponding to the fill type and ID
            if fill_type == 'material':
                if isinstance(fill_id, Iterable):
                    fill = [self._fast_materials[mat] if mat > 0 else None
                            for mat in fill_id]
                else:
                    fill = self._fast_materials[fill_id] if fill_id > 0 else None
            elif fill_type == 'universe':
                fill = self._fast_universes[fill_id]
                fill_univ_ids.add(fill_id)
            else:
                fill = self._fast_lattices[fill_id]
                for idx in fill._natural_indices:
                    univ = fill.get_universe(idx)
                    fill_univ_ids.add(univ.id)
                if fill.outer is not None:
                    fill_univ_ids.add(fill.outer.id)

            # Set the fill for the Cell
            self._fast_cells[cell_id].fill = fill

        # Determine root universe for geometry
        non_fill = set(self._fast_universes.keys()) - fill_univ_ids

        self.geometry.root_universe = self._fast_universes[non_fill.pop()]

    def add_volume_information(self, volume_calc):
        """Add volume information to the geometry within the summary file

        Parameters
        ----------
        volume_calc : openmc.VolumeCalculation
            Results from a stochastic volume calculation

        """
        if volume_calc.domain_type == "material" and self.materials:
            for material in self.materials:
                if material.id in volume_calc.volumes:
                    material.add_volume_information(volume_calc)
        else:
            self.geometry.add_volume_information(volume_calc)
