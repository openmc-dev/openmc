from collections import Iterable
import re

import numpy as np
import h5py

import openmc
import openmc.checkvalue as cv
from openmc.region import Region

_VERSION_SUMMARY = 5


class Summary(object):
    """Summary of geometry, materials, and tallies used in a simulation.

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
    version: tuple of int
        Version of OpenMC

    """

    def __init__(self, filename):
        openmc.reset_auto_ids()

        if not filename.endswith(('.h5', '.hdf5')):
            msg = 'Unable to open "{0}" which is not an HDF5 summary file'
            raise ValueError(msg)

        self._f = h5py.File(filename, 'r')
        cv.check_filetype_version(self._f, 'summary', _VERSION_SUMMARY)

        self._geometry = openmc.Geometry()
        self._materials = openmc.Materials()
        self._nuclides = {}

        self._read_nuclides()
        self._read_geometry()

    @property
    def date_and_time(self):
        return self._f.attrs['date_and_time']

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
    def version(self):
        return tuple(self._f.attrs['openmc_version'])

    def _read_nuclides(self):
        names = self._f['nuclides/names'].value
        awrs = self._f['nuclides/awrs'].value
        for name, awr in zip(names, awrs):
            self._nuclides[name.decode()] = awr

    def _read_geometry(self):
        # Read in and initialize the Materials and Geometry
        self._read_materials()
        surfaces = self._read_surfaces()
        cells, cell_fills = self._read_cells(surfaces)
        universes = self._read_universes(cells)
        lattices = self._read_lattices(universes)
        self._finalize_geometry(cells, cell_fills, universes, lattices)

    def _read_materials(self):
        for group in self._f['materials'].values():
            material = openmc.Material.from_hdf5(group)

            # Add the material to the Materials collection
            self.materials.append(material)

    def _read_surfaces(self):
        surfaces = {}
        for group in self._f['geometry/surfaces'].values():
            surface = openmc.Surface.from_hdf5(group)
            surfaces[surface.id] = surface

        return surfaces

    def _read_cells(self, surfaces):
        cells = {}

        # Initialize dictionary for each Cell's fill
        cell_fills = {}

        for key, group in self._f['geometry/cells'].items():
            cell_id = int(key.lstrip('cell '))
            name = group['name'].value.decode() if 'name' in group else ''
            fill_type = group['fill_type'].value.decode()

            if fill_type == 'material':
                fill = group['material'].value
            elif fill_type == 'universe':
                fill = group['fill'].value
            else:
                fill = group['lattice'].value

            region = group['region'].value.decode() if 'region' in group else ''

            # Create this Cell
            cell = openmc.Cell(cell_id=cell_id, name=name)

            if fill_type == 'universe':
                if 'offset' in group:
                    offset = group['offset'][...]
                    cell.offsets = offset

                if 'translation' in group:
                    translation = group['translation'][...]
                    translation = np.asarray(translation, dtype=np.float64)
                    cell.translation = translation

                if 'rotation' in group:
                    rotation = group['rotation'][...]
                    rotation = np.asarray(rotation, dtype=np.int)
                    cell._rotation = rotation

            elif fill_type == 'material':
                cell.temperature = group['temperature'][...]

            # Store Cell fill information for after Universe/Lattice creation
            cell_fills[cell.id] = (fill_type, fill)

            # Generate Region object given infix expression
            if region:
                cell.region = Region.from_expression(region, surfaces)

            # Get the distribcell data
            if 'distribcell_index' in group:
                ind = group['distribcell_index'].value
                cell.distribcell_index = ind
                paths = group['paths'][...]
                paths = [str(path.decode()) for path in paths]
                cell.distribcell_paths = paths

            # Add the Cell to the global dictionary of all Cells
            cells[cell.id] = cell

        return cells, cell_fills

    def _read_universes(self, cells):
        universes = {}
        for group in self._f['geometry/universes'].values():
            universe = openmc.Universe.from_hdf5(group, cells)
            universes[universe.id] = universe
        return universes

    def _read_lattices(self, universes):
        lattices = {}
        for group in self._f['geometry/lattices'].values():
            lattice = openmc.Lattice.from_hdf5(group, universes)
            lattices[lattice.id] = lattice
        return lattices

    def _finalize_geometry(self, cells, cell_fills, universes, lattices):
        materials = {m.id: m for m in self.materials}

        # Iterate over all Cells and add fill Materials, Universes and Lattices
        for cell_id, (fill_type, fill_id) in cell_fills.items():
            # Retrieve the object corresponding to the fill type and ID
            if fill_type == 'material':
                if isinstance(fill_id, Iterable):
                    fill = [materials[mat] if mat > 0 else None
                            for mat in fill_id]
                else:
                    fill = materials[fill_id] if fill_id > 0 else None
            elif fill_type == 'universe':
                fill = universes[fill_id]
            else:
                fill = lattices[fill_id]

            # Set the fill for the Cell
            cells[cell_id].fill = fill

        # Set the root universe for the Geometry
        self.geometry.root_universe = universes[0]

    def add_volume_information(self, volume_calc):
        """Add volume information to the geometry within the summary file

        Parameters
        ----------
        volume_calc : openmc.VolumeCalculation
            Results from a stochastic volume calculation

        """
        self.geometry.add_volume_information(volume_calc)
