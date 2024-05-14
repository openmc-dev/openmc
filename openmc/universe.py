import math
from abc import ABC, abstractmethod
from collections.abc import Iterable
from numbers import Integral, Real
from pathlib import Path
from tempfile import TemporaryDirectory
import warnings

import h5py
import lxml.etree as ET
import numpy as np

import openmc
import openmc.checkvalue as cv
from ._xml import get_text
from .checkvalue import check_type, check_value
from .mixin import IDManagerMixin
from .surface import _BOUNDARY_TYPES


class UniverseBase(ABC, IDManagerMixin):
    """A collection of cells that can be repeated.

    Attributes
    ----------
    id : int
        Unique identifier of the universe
    name : str
        Name of the universe
    """

    next_id = 1
    used_ids = set()

    def __init__(self, universe_id=None, name=''):
        # Initialize Universe class attributes
        self.id = universe_id
        self.name = name
        self._volume = None
        self._atoms = {}

        # Keys   - Cell IDs
        # Values - Cells
        self._cells = {}

    def __repr__(self):
        string = 'Universe\n'
        string += '{: <16}=\t{}\n'.format('\tID', self._id)
        string += '{: <16}=\t{}\n'.format('\tName', self._name)
        return string

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        if name is not None:
            cv.check_type('universe name', name, str)
            self._name = name
        else:
            self._name = ''

    @property
    def volume(self):
        return self._volume

    @volume.setter
    def volume(self, volume):
        if volume is not None:
            cv.check_type('universe volume', volume, Real)
        self._volume = volume

    def add_volume_information(self, volume_calc):
        """Add volume information to a universe.

        Parameters
        ----------
        volume_calc : openmc.VolumeCalculation
            Results from a stochastic volume calculation

        """
        if volume_calc.domain_type == 'universe':
            if self.id in volume_calc.volumes:
                self._volume = volume_calc.volumes[self.id].n
                self._atoms = volume_calc.atoms[self.id]
            else:
                raise ValueError(
                    'No volume information found for this universe.')
        else:
            raise ValueError('No volume information found for this universe.')

    def get_all_universes(self):
        """Return all universes that are contained within this one.

        Returns
        -------
        universes : dict
            Dictionary whose keys are universe IDs and values are
            :class:`Universe` instances

        """
        # Append all Universes within each Cell to the dictionary
        universes = {}
        for cell in self.get_all_cells().values():
            universes.update(cell.get_all_universes())

        return universes

    @abstractmethod
    def create_xml_subelement(self, xml_element, memo=None):
        """Add the universe xml representation to an incoming xml element

        Parameters
        ----------
        xml_element : lxml.etree._Element
            XML element to be added to

        memo : set or None
            A set of object id's representing geometry entities already
            written to the xml_element. This parameter is used internally
            and should not be specified by users.

        Returns
        -------
        None

        """

    @abstractmethod
    def _partial_deepcopy(self):
        """Deepcopy all parameters of an openmc.UniverseBase object except its cells.
        This should only be used from the openmc.UniverseBase.clone() context.

        """

    def clone(self, clone_materials=True, clone_regions=True, memo=None):
        """Create a copy of this universe with a new unique ID, and clones
        all cells within this universe.

        Parameters
        ----------
        clone_materials : bool
            Whether to create separates copies of the materials filling cells
            contained in this universe.
        clone_regions : bool
            Whether to create separates copies of the regions bounding cells
            contained in this universe.
        memo : dict or None
            A nested dictionary of previously cloned objects. This parameter
            is used internally and should not be specified by the user.

        Returns
        -------
        clone : openmc.Universe
            The clone of this universe

        """
        if memo is None:
            memo = {}

        # If no memoize'd clone exists, instantiate one
        if self not in memo:
            clone = self._partial_deepcopy()

            # Clone all cells for the universe clone
            clone._cells = {}
            for cell in self._cells.values():
                clone.add_cell(cell.clone(clone_materials, clone_regions,
                                          memo))

            # Memoize the clone
            memo[self] = clone

        return memo[self]


class Universe(UniverseBase):
    """A collection of cells that can be repeated.

    Parameters
    ----------
    universe_id : int, optional
        Unique identifier of the universe. If not specified, an identifier will
        automatically be assigned
    name : str, optional
        Name of the universe. If not specified, the name is the empty string.
    cells : Iterable of openmc.Cell, optional
        Cells to add to the universe. By default no cells are added.

    Attributes
    ----------
    id : int
        Unique identifier of the universe
    name : str
        Name of the universe
    cells : dict
        Dictionary whose keys are cell IDs and values are :class:`Cell`
        instances
    volume : float
        Volume of the universe in cm^3. This can either be set manually or
        calculated in a stochastic volume calculation and added via the
        :meth:`Universe.add_volume_information` method.
    bounding_box : openmc.BoundingBox
        Lower-left and upper-right coordinates of an axis-aligned bounding box
        of the universe.

    """

    def __init__(self, universe_id=None, name='', cells=None):
        super().__init__(universe_id, name)

        if cells is not None:
            self.add_cells(cells)

    def __repr__(self):
        string = super().__repr__()
        string += '{: <16}=\t{}\n'.format('\tGeom', 'CSG')
        string += '{: <16}=\t{}\n'.format('\tCells', list(self._cells.keys()))
        return string

    @property
    def cells(self):
        return self._cells

    @property
    def bounding_box(self):
        regions = [c.region for c in self.cells.values()
                   if c.region is not None]
        if regions:
            return openmc.Union(regions).bounding_box
        else:
            return openmc.BoundingBox.infinite()

    @classmethod
    def from_hdf5(cls, group, cells):
        """Create universe from HDF5 group

        Parameters
        ----------
        group : h5py.Group
            Group in HDF5 file
        cells : dict
            Dictionary mapping cell IDs to instances of :class:`openmc.Cell`.

        Returns
        -------
        openmc.Universe
            Universe instance

        """
        universe_id = int(group.name.split('/')[-1].lstrip('universe '))
        cell_ids = group['cells'][()]

        # Create this Universe
        universe = cls(universe_id)

        # Add each Cell to the Universe
        for cell_id in cell_ids:
            universe.add_cell(cells[cell_id])

        return universe

    def find(self, point):
        """Find cells/universes/lattices which contain a given point

        Parameters
        ----------
        point : 3-tuple of float
            Cartesian coordinates of the point

        Returns
        -------
        list
            Sequence of universes, cells, and lattices which are traversed to
            find the given point

        """
        p = np.asarray(point)
        for cell in self._cells.values():
            if p in cell:
                if cell.fill_type in ('material', 'distribmat', 'void'):
                    return [self, cell]
                elif cell.fill_type == 'universe':
                    if cell.translation is not None:
                        p -= cell.translation
                    if cell.rotation is not None:
                        p[:] = cell.rotation_matrix.dot(p)
                    return [self, cell] + cell.fill.find(p)
                else:
                    return [self, cell] + cell.fill.find(p)
        return []

    # default kwargs that are passed to plt.legend in the plot method below.
    _default_legend_kwargs = {'bbox_to_anchor': (
        1.05, 1), 'loc': 2, 'borderaxespad': 0.0}

    def plot(self, origin=None, width=None, pixels=40000,
             basis='xy', color_by='cell', colors=None, seed=None,
             openmc_exec='openmc', axes=None, legend=False, axis_units='cm',
             legend_kwargs=_default_legend_kwargs, outline=False,
             **kwargs):
        """Display a slice plot of the universe.

        Parameters
        ----------
        origin : iterable of float
            Coordinates at the origin of the plot. If left as None,
            universe.bounding_box.center will be used to attempt to ascertain
            the origin with infinite values being replaced by 0.
        width : iterable of float
            Width of the plot in each basis direction. If left as none then the
            universe.bounding_box.width() will be used to attempt to
            ascertain the plot width.  Defaults to (10, 10) if the bounding_box
            contains inf values
        pixels : Iterable of int or int
            If iterable of ints provided then this directly sets the number of
            pixels to use in each basis direction. If int provided then this
            sets the total number of pixels in the plot and the number of
            pixels in each basis direction is calculated from this total and
            the image aspect ratio.
        basis : {'xy', 'xz', 'yz'}
            The basis directions for the plot
        color_by : {'cell', 'material'}
            Indicate whether the plot should be colored by cell or by material
        colors : dict
            Assigns colors to specific materials or cells. Keys are instances of
            :class:`Cell` or :class:`Material` and values are RGB 3-tuples, RGBA
            4-tuples, or strings indicating SVG color names. Red, green, blue,
            and alpha should all be floats in the range [0.0, 1.0], for example:

            .. code-block:: python

               # Make water blue
               water = openmc.Cell(fill=h2o)
               universe.plot(..., colors={water: (0., 0., 1.))
        seed : int
            Seed for the random number generator
        openmc_exec : str
            Path to OpenMC executable.
        axes : matplotlib.Axes
            Axes to draw to

            .. versionadded:: 0.13.1
        legend : bool
            Whether a legend showing material or cell names should be drawn

            .. versionadded:: 0.14.0
        legend_kwargs : dict
            Keyword arguments passed to :func:`matplotlib.pyplot.legend`.

            .. versionadded:: 0.14.0
        outline : bool
            Whether outlines between color boundaries should be drawn

            .. versionadded:: 0.14.0
        axis_units : {'km', 'm', 'cm', 'mm'}
            Units used on the plot axis

            .. versionadded:: 0.14.0
        **kwargs
            Keyword arguments passed to :func:`matplotlib.pyplot.imshow`

        Returns
        -------
        matplotlib.axes.Axes
            Axes containing resulting image

        """
        import matplotlib.image as mpimg
        import matplotlib.patches as mpatches
        import matplotlib.pyplot as plt

        # Determine extents of plot
        if basis == 'xy':
            x, y = 0, 1
            xlabel, ylabel = f'x [{axis_units}]', f'y [{axis_units}]'
        elif basis == 'yz':
            x, y = 1, 2
            xlabel, ylabel = f'y [{axis_units}]', f'z [{axis_units}]'
        elif basis == 'xz':
            x, y = 0, 2
            xlabel, ylabel = f'x [{axis_units}]', f'z [{axis_units}]'

        bb = self.bounding_box
        # checks to see if bounding box contains -inf or inf values
        if np.isinf(bb.extent[basis]).any():
            if origin is None:
                origin = (0, 0, 0)
            if width is None:
                width = (10, 10)
        else:
            if origin is None:
                # if nan values in the bb.center they get replaced with 0.0
                # this happens when the bounding_box contains inf values
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    origin = np.nan_to_num(bb.center)
            if width is None:
                bb_width = bb.width
                x_width = bb_width['xyz'.index(basis[0])]
                y_width = bb_width['xyz'.index(basis[1])]
                width = (x_width, y_width)

        if isinstance(pixels, int):
            aspect_ratio = width[0] / width[1]
            pixels_y = math.sqrt(pixels / aspect_ratio)
            pixels = (int(pixels / pixels_y), int(pixels_y))

        axis_scaling_factor = {'km': 0.00001, 'm': 0.01, 'cm': 1, 'mm': 10}

        x_min = (origin[x] - 0.5*width[0]) * axis_scaling_factor[axis_units]
        x_max = (origin[x] + 0.5*width[0]) * axis_scaling_factor[axis_units]
        y_min = (origin[y] - 0.5*width[1]) * axis_scaling_factor[axis_units]
        y_max = (origin[y] + 0.5*width[1]) * axis_scaling_factor[axis_units]

        with TemporaryDirectory() as tmpdir:
            model = openmc.Model()
            model.geometry = openmc.Geometry(self)
            if seed is not None:
                model.settings.plot_seed = seed

            # Determine whether any materials contains macroscopic data and if
            # so, set energy mode accordingly
            for mat in self.get_all_materials().values():
                if mat._macroscopic is not None:
                    model.settings.energy_mode = 'multi-group'
                    break

            # Create plot object matching passed arguments
            plot = openmc.Plot()
            plot.origin = origin
            plot.width = width
            plot.pixels = pixels
            plot.basis = basis
            plot.color_by = color_by
            if colors is not None:
                plot.colors = colors
            model.plots.append(plot)

            # Run OpenMC in geometry plotting mode
            model.plot_geometry(False, cwd=tmpdir, openmc_exec=openmc_exec)

            # Read image from file
            img_path = Path(tmpdir) / f'plot_{plot.id}.png'
            if not img_path.is_file():
                img_path = img_path.with_suffix('.ppm')
            img = mpimg.imread(str(img_path))

            # Create a figure sized such that the size of the axes within
            # exactly matches the number of pixels specified
            if axes is None:
                px = 1/plt.rcParams['figure.dpi']
                fig, axes = plt.subplots()
                axes.set_xlabel(xlabel)
                axes.set_ylabel(ylabel)
                params = fig.subplotpars
                width = pixels[0]*px/(params.right - params.left)
                height = pixels[1]*px/(params.top - params.bottom)
                fig.set_size_inches(width, height)

            if outline:
                # Combine R, G, B values into a single int
                rgb = (img * 256).astype(int)
                image_value = (rgb[..., 0] << 16) + \
                    (rgb[..., 1] << 8) + (rgb[..., 2])

                axes.contour(
                    image_value,
                    origin="upper",
                    colors="k",
                    linestyles="solid",
                    linewidths=1,
                    levels=np.unique(image_value),
                    extent=(x_min, x_max, y_min, y_max),
                )

            # add legend showing which colors represent which material
            # or cell if that was requested
            if legend:
                if plot.colors == {}:
                    raise ValueError("Must pass 'colors' dictionary if you "
                                     "are adding a legend via legend=True.")

                if color_by == "cell":
                    expected_key_type = openmc.Cell
                else:
                    expected_key_type = openmc.Material

                patches = []
                for key, color in plot.colors.items():

                    if isinstance(key, int):
                        raise TypeError(
                            "Cannot use IDs in colors dict for auto legend.")
                    elif not isinstance(key, expected_key_type):
                        raise TypeError(
                            "Color dict key type does not match color_by")

                    # this works whether we're doing cells or materials
                    label = key.name if key.name != '' else key.id

                    # matplotlib takes RGB on 0-1 scale rather than 0-255. at
                    # this point PlotBase has already checked that 3-tuple
                    # based colors are already valid, so if the length is three
                    # then we know it just needs to be converted to the 0-1
                    # format.
                    if len(color) == 3 and not isinstance(color, str):
                        scaled_color = (
                            color[0]/255, color[1]/255, color[2]/255)
                    else:
                        scaled_color = color

                    key_patch = mpatches.Patch(color=scaled_color, label=label)
                    patches.append(key_patch)

                axes.legend(handles=patches, **legend_kwargs)

            # Plot image and return the axes
            axes.imshow(img, extent=(x_min, x_max, y_min, y_max), **kwargs)
            return axes

    def add_cell(self, cell):
        """Add a cell to the universe.

        Parameters
        ----------
        cell : openmc.Cell
            Cell to add

        """

        if not isinstance(cell, openmc.Cell):
            msg = f'Unable to add a Cell to Universe ID="{self._id}" since ' \
                  f'"{cell}" is not a Cell'
            raise TypeError(msg)

        cell_id = cell.id

        if cell_id not in self._cells:
            self._cells[cell_id] = cell

    def add_cells(self, cells):
        """Add multiple cells to the universe.

        Parameters
        ----------
        cells : Iterable of openmc.Cell
            Cells to add

        """

        if not isinstance(cells, Iterable):
            msg = f'Unable to add Cells to Universe ID="{self._id}" since ' \
                  f'"{cells}" is not iterable'
            raise TypeError(msg)

        for cell in cells:
            self.add_cell(cell)

    def remove_cell(self, cell):
        """Remove a cell from the universe.

        Parameters
        ----------
        cell : openmc.Cell
            Cell to remove

        """

        if not isinstance(cell, openmc.Cell):
            msg = f'Unable to remove a Cell from Universe ID="{self._id}" ' \
                  f'since "{cell}" is not a Cell'
            raise TypeError(msg)

        # If the Cell is in the Universe's list of Cells, delete it
        self._cells.pop(cell.id, None)

    def clear_cells(self):
        """Remove all cells from the universe."""

        self._cells.clear()

    def get_nuclides(self):
        """Returns all nuclides in the universe

        Returns
        -------
        nuclides : list of str
            List of nuclide names

        """

        nuclides = []

        # Append all Nuclides in each Cell in the Universe to the dictionary
        for cell in self.cells.values():
            for nuclide in cell.get_nuclides():
                if nuclide not in nuclides:
                    nuclides.append(nuclide)

        return nuclides

    def get_nuclide_densities(self):
        """Return all nuclides contained in the universe

        Returns
        -------
        nuclides : dict
            Dictionary whose keys are nuclide names and values are 2-tuples of
            (nuclide, density)

        """
        nuclides = {}

        if self._atoms:
            volume = self.volume
            for name, atoms in self._atoms.items():
                nuclide = openmc.Nuclide(name)
                density = 1.0e-24 * atoms.n/volume  # density in atoms/b-cm
                nuclides[name] = (nuclide, density)
        else:
            raise RuntimeError(
                'Volume information is needed to calculate microscopic cross '
                f'sections for universe {self.id}. This can be done by running '
                'a stochastic volume calculation via the '
                'openmc.VolumeCalculation object')

        return nuclides

    def get_all_cells(self, memo=None):
        """Return all cells that are contained within the universe

        Returns
        -------
        cells : dict
            Dictionary whose keys are cell IDs and values are :class:`Cell`
            instances

        """

        cells = {}

        if memo and self in memo:
            return cells

        if memo is not None:
            memo.add(self)

        # Add this Universe's cells to the dictionary
        cells.update(self._cells)

        # Append all Cells in each Cell in the Universe to the dictionary
        for cell in self._cells.values():
            cells.update(cell.get_all_cells(memo))

        return cells

    def get_all_materials(self, memo=None):
        """Return all materials that are contained within the universe

        Returns
        -------
        materials : dict
            Dictionary whose keys are material IDs and values are
            :class:`Material` instances

        """

        materials = {}

        # Append all Cells in each Cell in the Universe to the dictionary
        cells = self.get_all_cells(memo)
        for cell in cells.values():
            materials.update(cell.get_all_materials(memo))

        return materials

    def create_xml_subelement(self, xml_element, memo=None):
        # Iterate over all Cells
        for cell in self._cells.values():

            # If the cell was already written, move on
            if memo and cell in memo:
                continue

            if memo is not None:
                memo.add(cell)

            # Create XML subelement for this Cell
            cell_element = cell.create_xml_subelement(xml_element, memo)

            # Append the Universe ID to the subelement and add to Element
            cell_element.set("universe", str(self._id))
            xml_element.append(cell_element)

    def _determine_paths(self, path='', instances_only=False):
        """Count the number of instances for each cell in the universe, and
        record the count in the :attr:`Cell.num_instances` properties."""

        univ_path = path + f'u{self.id}'

        for cell in self.cells.values():
            cell_path = f'{univ_path}->c{cell.id}'
            fill = cell._fill
            fill_type = cell.fill_type

            # If universe-filled, recursively count cells in filling universe
            if fill_type == 'universe':
                fill._determine_paths(cell_path + '->', instances_only)

            # If lattice-filled, recursively call for all universes in lattice
            elif fill_type == 'lattice':
                latt = fill

                # Count instances in each universe in the lattice
                for index in latt._natural_indices:
                    latt_path = '{}->l{}({})->'.format(
                        cell_path, latt.id, ",".join(str(x) for x in index))
                    univ = latt.get_universe(index)
                    univ._determine_paths(latt_path, instances_only)

            else:
                if fill_type == 'material':
                    mat = fill
                elif fill_type == 'distribmat':
                    mat = fill[cell._num_instances]
                else:
                    mat = None

                if mat is not None:
                    mat._num_instances += 1
                    if not instances_only:
                        mat._paths.append(f'{cell_path}->m{mat.id}')

            # Append current path
            cell._num_instances += 1
            if not instances_only:
                cell._paths.append(cell_path)

    def _partial_deepcopy(self):
        """Clone all of the openmc.Universe object's attributes except for its cells,
        as they are copied within the clone function. This should only to be
        used within the openmc.UniverseBase.clone() context.
        """
        clone = openmc.Universe(name=self.name)
        clone.volume = self.volume
        return clone


class DAGMCUniverse(UniverseBase):
    """A reference to a DAGMC file to be used in the model.

    .. versionadded:: 0.13.0

    Parameters
    ----------
    filename : str
        Path to the DAGMC file used to represent this universe.
    universe_id : int, optional
        Unique identifier of the universe. If not specified, an identifier will
        automatically be assigned.
    name : str, optional
        Name of the universe. If not specified, the name is the empty string.
    auto_geom_ids : bool
        Set IDs automatically on initialization (True) or report overlaps in ID
        space between CSG and DAGMC (False)
    auto_mat_ids : bool
        Set IDs automatically on initialization (True)  or report overlaps in ID
        space between OpenMC and UWUW materials (False)

    Attributes
    ----------
    id : int
        Unique identifier of the universe
    name : str
        Name of the universe
    filename : str
        Path to the DAGMC file used to represent this universe.
    auto_geom_ids : bool
        Set IDs automatically on initialization (True) or report overlaps in ID
        space between CSG and DAGMC (False)
    auto_mat_ids : bool
        Set IDs automatically on initialization (True)  or report overlaps in ID
        space between OpenMC and UWUW materials (False)
    bounding_box : openmc.BoundingBox
        Lower-left and upper-right coordinates of an axis-aligned bounding box
        of the universe.

        .. versionadded:: 0.13.1
    material_names : list of str
        Return a sorted list of materials names that are contained within the
        DAGMC h5m file. This is useful when naming openmc.Material() objects
        as each material name present in the DAGMC h5m file must have a
        matching openmc.Material() with the same name.

        .. versionadded:: 0.13.2
    n_cells : int
        The number of cells in the DAGMC model. This is the number of cells at
        runtime and accounts for the implicit complement whether or not is it
        present in the DAGMC file.

        .. versionadded:: 0.13.2
    n_surfaces : int
        The number of surfaces in the model.

        .. versionadded:: 0.13.2

    """

    def __init__(self,
                 filename,
                 universe_id=None,
                 name='',
                 auto_geom_ids=False,
                 auto_mat_ids=False):
        super().__init__(universe_id, name)
        # Initialize class attributes
        self.filename = filename
        self.auto_geom_ids = auto_geom_ids
        self.auto_mat_ids = auto_mat_ids

    def __repr__(self):
        string = super().__repr__()
        string += '{: <16}=\t{}\n'.format('\tGeom', 'DAGMC')
        string += '{: <16}=\t{}\n'.format('\tFile', self.filename)
        return string

    @property
    def bounding_box(self):
        with h5py.File(self.filename) as dagmc_file:
            coords = dagmc_file['tstt']['nodes']['coordinates'][()]
            lower_left_corner = coords.min(axis=0)
            upper_right_corner = coords.max(axis=0)
            return openmc.BoundingBox(lower_left_corner, upper_right_corner)

    @property
    def filename(self):
        return self._filename

    @filename.setter
    def filename(self, val):
        cv.check_type('DAGMC filename', val, (Path, str))
        self._filename = val

    @property
    def auto_geom_ids(self):
        return self._auto_geom_ids

    @auto_geom_ids.setter
    def auto_geom_ids(self, val):
        cv.check_type('DAGMC automatic geometry ids', val, bool)
        self._auto_geom_ids = val

    @property
    def auto_mat_ids(self):
        return self._auto_mat_ids

    @auto_mat_ids.setter
    def auto_mat_ids(self, val):
        cv.check_type('DAGMC automatic material ids', val, bool)
        self._auto_mat_ids = val

    @property
    def material_names(self):
        dagmc_file_contents = h5py.File(self.filename)
        material_tags_hex = dagmc_file_contents['/tstt/tags/NAME'].get(
            'values')
        material_tags_ascii = []
        for tag in material_tags_hex:
            candidate_tag = tag.tobytes().decode().replace('\x00', '')
            # tags might be for temperature or reflective surfaces
            if candidate_tag.startswith('mat:'):
                # removes first 4 characters as openmc.Material name should be
                # set without the 'mat:' part of the tag
                material_tags_ascii.append(candidate_tag[4:])

        return sorted(set(material_tags_ascii))

    def get_all_cells(self, memo=None):
        return {}

    def get_all_materials(self, memo=None):
        return {}

    def _n_geom_elements(self, geom_type):
        """
        Helper function for retrieving the number geometric entities in a DAGMC
        file

        Parameters
        ----------
        geom_type : str
            The type of geometric entity to count. One of {'Volume', 'Surface'}. Returns
            the runtime number of voumes in the DAGMC model (includes implicit complement).

        Returns
        -------
        int
            Number of geometry elements of the specified type
        """
        cv.check_value('geometry type', geom_type, ('volume', 'surface'))

        def decode_str_tag(tag_val):
            return tag_val.tobytes().decode().replace('\x00', '')

        dagmc_filepath = Path(self.filename).resolve()
        with h5py.File(dagmc_filepath) as dagmc_file:
            category_data = dagmc_file['tstt/tags/CATEGORY/values']
            category_strs = map(decode_str_tag, category_data)
            n = sum([v == geom_type.capitalize() for v in category_strs])

            # check for presence of an implicit complement in the file and
            # increment the number of cells if it doesn't exist
            if geom_type == 'volume':
                name_data = dagmc_file['tstt/tags/NAME/values']
                name_strs = map(decode_str_tag, name_data)
                if not sum(['impl_complement' in n for n in name_strs]):
                    n += 1
        return n

    @property
    def n_cells(self):
        return self._n_geom_elements('volume')

    @property
    def n_surfaces(self):
        return self._n_geom_elements('surface')

    def create_xml_subelement(self, xml_element, memo=None):
        if memo and self in memo:
            return

        if memo is not None:
            memo.add(self)

        # Set xml element values
        dagmc_element = ET.Element('dagmc_universe')
        dagmc_element.set('id', str(self.id))

        if self.auto_geom_ids:
            dagmc_element.set('auto_geom_ids', 'true')
        if self.auto_mat_ids:
            dagmc_element.set('auto_mat_ids', 'true')
        dagmc_element.set('filename', str(self.filename))
        xml_element.append(dagmc_element)

    def bounding_region(
            self,
            bounded_type: str = 'box',
            boundary_type: str = 'vacuum',
            starting_id: int = 10000,
            padding_distance: float = 0.
        ):
        """Creates a either a spherical or box shaped bounding region around
        the DAGMC geometry.

        .. versionadded:: 0.13.1

        Parameters
        ----------
        bounded_type : str
            The type of bounding surface(s) to use when constructing the region.
            Options include a single spherical surface (sphere) or a rectangle
            made from six planes (box).
        boundary_type : str
            Boundary condition that defines the behavior for particles hitting
            the surface. Defaults to vacuum boundary condition. Passed into the
            surface construction.
        starting_id : int
            Starting ID of the surface(s) used in the region. For bounded_type
            'box', the next 5 IDs will also be used. Defaults to 10000 to reduce
            the chance of an overlap of surface IDs with the DAGMC geometry.
        padding_distance : float
            Distance between the bounding region surfaces and the minimal
            bounding box. Allows for the region to be larger than the DAGMC
            geometry.

        Returns
        -------
        openmc.Region
            Region instance
        """

        check_type('boundary type', boundary_type, str)
        check_value('boundary type', boundary_type, _BOUNDARY_TYPES)
        check_type('starting surface id', starting_id, Integral)
        check_type('bounded type', bounded_type, str)
        check_value('bounded type', bounded_type, ('box', 'sphere'))

        bbox = self.bounding_box.expand(padding_distance, True)

        if bounded_type == 'sphere':
            radius = np.linalg.norm(bbox.upper_right - bbox.center)
            bounding_surface = openmc.Sphere(
                surface_id=starting_id,
                x0=bbox.center[0],
                y0=bbox.center[1],
                z0=bbox.center[2],
                boundary_type=boundary_type,
                r=radius,
            )

            return -bounding_surface

        if bounded_type == 'box':
            # defines plane surfaces for all six faces of the bounding box
            lower_x = openmc.XPlane(bbox[0][0], surface_id=starting_id)
            upper_x = openmc.XPlane(bbox[1][0], surface_id=starting_id+1)
            lower_y = openmc.YPlane(bbox[0][1], surface_id=starting_id+2)
            upper_y = openmc.YPlane(bbox[1][1], surface_id=starting_id+3)
            lower_z = openmc.ZPlane(bbox[0][2], surface_id=starting_id+4)
            upper_z = openmc.ZPlane(bbox[1][2], surface_id=starting_id+5)

            region = +lower_x & -upper_x & +lower_y & -upper_y & +lower_z & -upper_z

            for surface in region.get_surfaces().values():
                surface.boundary_type = boundary_type

            return region

    def bounded_universe(self, bounding_cell_id=10000, **kwargs):
        """Returns an openmc.Universe filled with this DAGMCUniverse and bounded
        with a cell. Defaults to a box cell with a vacuum surface however this
        can be changed using the kwargs which are passed directly to
        DAGMCUniverse.bounding_region().

        Parameters
        ----------
        bounding_cell_id : int
            The cell ID number to use for the bounding cell, defaults to 10000 to reduce
            the chance of overlapping ID numbers with the DAGMC geometry.

        Returns
        -------
        openmc.Universe
            Universe instance
        """
        bounding_cell = openmc.Cell(
            fill=self, cell_id=bounding_cell_id, region=self.bounding_region(**kwargs))
        return openmc.Universe(cells=[bounding_cell])

    @classmethod
    def from_hdf5(cls, group):
        """Create DAGMC universe from HDF5 group

        Parameters
        ----------
        group : h5py.Group
            Group in HDF5 file

        Returns
        -------
        openmc.DAGMCUniverse
            DAGMCUniverse instance

        """
        id = int(group.name.split('/')[-1].lstrip('universe '))
        fname = group['filename'][()].decode()
        name = group['name'][()].decode() if 'name' in group else None

        out = cls(fname, universe_id=id, name=name)

        out.auto_geom_ids = bool(group.attrs['auto_geom_ids'])
        out.auto_mat_ids = bool(group.attrs['auto_mat_ids'])

        return out

    @classmethod
    def from_xml_element(cls, elem):
        """Generate DAGMC universe from XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            `<dagmc_universe>` element

        Returns
        -------
        openmc.DAGMCUniverse
            DAGMCUniverse instance

        """
        id = int(get_text(elem, 'id'))
        fname = get_text(elem, 'filename')

        out = cls(fname, universe_id=id)

        name = get_text(elem, 'name')
        if name is not None:
            out.name = name

        out.auto_geom_ids = bool(elem.get('auto_geom_ids'))
        out.auto_mat_ids = bool(elem.get('auto_mat_ids'))

        return out

    def _partial_deepcopy(self):
        """Clone all of the openmc.DAGMCUniverse object's attributes except for
        its cells, as they are copied within the clone function. This should
        only to be used within the openmc.UniverseBase.clone() context.
        """
        clone = openmc.DAGMCUniverse(name=self.name, filename=self.filename)
        clone.volume = self.volume
        clone.auto_geom_ids = self.auto_geom_ids
        clone.auto_mat_ids = self.auto_mat_ids
        return clone
