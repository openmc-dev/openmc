from collections.abc import Iterable, Mapping
from numbers import Integral, Real
from pathlib import Path

import h5py
import lxml.etree as ET
import numpy as np

import openmc
import openmc.checkvalue as cv
from openmc.checkvalue import PathLike

from ._xml import clean_indentation, get_elem_tuple, reorder_attributes, get_text
from .mixin import IDManagerMixin

_BASES = {'xy', 'xz', 'yz'}

_SVG_COLORS = {
    'aliceblue': (240, 248, 255),
    'antiquewhite': (250, 235, 215),
    'aqua': (0, 255, 255),
    'aquamarine': (127, 255, 212),
    'azure': (240, 255, 255),
    'beige': (245, 245, 220),
    'bisque': (255, 228, 196),
    'black': (0, 0, 0),
    'blanchedalmond': (255, 235, 205),
    'blue': (0, 0, 255),
    'blueviolet': (138, 43, 226),
    'brown': (165, 42, 42),
    'burlywood': (222, 184, 135),
    'cadetblue': (95, 158, 160),
    'chartreuse': (127, 255, 0),
    'chocolate': (210, 105, 30),
    'coral': (255, 127, 80),
    'cornflowerblue': (100, 149, 237),
    'cornsilk': (255, 248, 220),
    'crimson': (220, 20, 60),
    'cyan': (0, 255, 255),
    'darkblue': (0, 0, 139),
    'darkcyan': (0, 139, 139),
    'darkgoldenrod': (184, 134, 11),
    'darkgray': (169, 169, 169),
    'darkgreen': (0, 100, 0),
    'darkgrey': (169, 169, 169),
    'darkkhaki': (189, 183, 107),
    'darkmagenta': (139, 0, 139),
    'darkolivegreen': (85, 107, 47),
    'darkorange': (255, 140, 0),
    'darkorchid': (153, 50, 204),
    'darkred': (139, 0, 0),
    'darksalmon': (233, 150, 122),
    'darkseagreen': (143, 188, 143),
    'darkslateblue': (72, 61, 139),
    'darkslategray': (47, 79, 79),
    'darkslategrey': (47, 79, 79),
    'darkturquoise': (0, 206, 209),
    'darkviolet': (148, 0, 211),
    'deeppink': (255, 20, 147),
    'deepskyblue': (0, 191, 255),
    'dimgray': (105, 105, 105),
    'dimgrey': (105, 105, 105),
    'dodgerblue': (30, 144, 255),
    'firebrick': (178, 34, 34),
    'floralwhite': (255, 250, 240),
    'forestgreen': (34, 139, 34),
    'fuchsia': (255, 0, 255),
    'gainsboro': (220, 220, 220),
    'ghostwhite': (248, 248, 255),
    'gold': (255, 215, 0),
    'goldenrod': (218, 165, 32),
    'gray': (128, 128, 128),
    'green': (0, 128, 0),
    'greenyellow': (173, 255, 47),
    'grey': (128, 128, 128),
    'honeydew': (240, 255, 240),
    'hotpink': (255, 105, 180),
    'indianred': (205, 92, 92),
    'indigo': (75, 0, 130),
    'ivory': (255, 255, 240),
    'khaki': (240, 230, 140),
    'lavender': (230, 230, 250),
    'lavenderblush': (255, 240, 245),
    'lawngreen': (124, 252, 0),
    'lemonchiffon': (255, 250, 205),
    'lightblue': (173, 216, 230),
    'lightcoral': (240, 128, 128),
    'lightcyan': (224, 255, 255),
    'lightgoldenrodyellow': (250, 250, 210),
    'lightgray': (211, 211, 211),
    'lightgreen': (144, 238, 144),
    'lightgrey': (211, 211, 211),
    'lightpink': (255, 182, 193),
    'lightsalmon': (255, 160, 122),
    'lightseagreen': (32, 178, 170),
    'lightskyblue': (135, 206, 250),
    'lightslategray': (119, 136, 153),
    'lightslategrey': (119, 136, 153),
    'lightsteelblue': (176, 196, 222),
    'lightyellow': (255, 255, 224),
    'lime': (0, 255, 0),
    'limegreen': (50, 205, 50),
    'linen': (250, 240, 230),
    'magenta': (255, 0, 255),
    'maroon': (128, 0, 0),
    'mediumaquamarine': (102, 205, 170),
    'mediumblue': (0, 0, 205),
    'mediumorchid': (186, 85, 211),
    'mediumpurple': (147, 112, 219),
    'mediumseagreen': (60, 179, 113),
    'mediumslateblue': (123, 104, 238),
    'mediumspringgreen': (0, 250, 154),
    'mediumturquoise': (72, 209, 204),
    'mediumvioletred': (199, 21, 133),
    'midnightblue': (25, 25, 112),
    'mintcream': (245, 255, 250),
    'mistyrose': (255, 228, 225),
    'moccasin': (255, 228, 181),
    'navajowhite': (255, 222, 173),
    'navy': (0, 0, 128),
    'oldlace': (253, 245, 230),
    'olive': (128, 128, 0),
    'olivedrab': (107, 142, 35),
    'orange': (255, 165, 0),
    'orangered': (255, 69, 0),
    'orchid': (218, 112, 214),
    'palegoldenrod': (238, 232, 170),
    'palegreen': (152, 251, 152),
    'paleturquoise': (175, 238, 238),
    'palevioletred': (219, 112, 147),
    'papayawhip': (255, 239, 213),
    'peachpuff': (255, 218, 185),
    'peru': (205, 133, 63),
    'pink': (255, 192, 203),
    'plum': (221, 160, 221),
    'powderblue': (176, 224, 230),
    'purple': (128, 0, 128),
    'red': (255, 0, 0),
    'rosybrown': (188, 143, 143),
    'royalblue': (65, 105, 225),
    'saddlebrown': (139, 69, 19),
    'salmon': (250, 128, 114),
    'sandybrown': (244, 164, 96),
    'seagreen': (46, 139, 87),
    'seashell': (255, 245, 238),
    'sienna': (160, 82, 45),
    'silver': (192, 192, 192),
    'skyblue': (135, 206, 235),
    'slateblue': (106, 90, 205),
    'slategray': (112, 128, 144),
    'slategrey': (112, 128, 144),
    'snow': (255, 250, 250),
    'springgreen': (0, 255, 127),
    'steelblue': (70, 130, 180),
    'tan': (210, 180, 140),
    'teal': (0, 128, 128),
    'thistle': (216, 191, 216),
    'tomato': (255, 99, 71),
    'turquoise': (64, 224, 208),
    'violet': (238, 130, 238),
    'wheat': (245, 222, 179),
    'white': (255, 255, 255),
    'whitesmoke': (245, 245, 245),
    'yellow': (255, 255, 0),
    'yellowgreen': (154, 205, 50)
}

_PLOT_PARAMS = """
        Parameters
        ----------
        origin : iterable of float
            Coordinates at the origin of the plot. If left as None,
            the center of the bounding box will be used to attempt to ascertain
            the origin with infinite values being replaced by 0.
        width : iterable of float
            Width of the plot in each basis direction. If left as none then the
            width of the bounding box will be used to attempt to
            ascertain the plot width. Defaults to (10, 10) if the bounding box
            contains inf values.
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
        axis_units : {'km', 'm', 'cm', 'mm'}
            Units used on the plot axis

            .. versionadded:: 0.14.0
        outline : bool or str
            Whether outlines between color boundaries should be drawn. If set to
            'only', only outlines will be drawn.

            .. versionadded:: 0.14.0
        show_overlaps: bool
            Indicate whether or not overlapping regions are shown.
            Default is False.
        overlap_color: Iterable of int or str
            Color to apply to overlapping regions. Default is red.
        n_samples : int, optional
            The number of source particles to sample and add to plot. Defaults
            to None which doesn't plot any particles on the plot.
        plane_tolerance: float
            When plotting a plane the source locations within the plane +/-
            the plane_tolerance will be included and those outside of the
            plane_tolerance will not be shown
        legend_kwargs : dict
            Keyword arguments passed to :func:`matplotlib.pyplot.legend`.

            .. versionadded:: 0.14.0
        source_kwargs : dict, optional
            Keyword arguments passed to :func:`matplotlib.pyplot.scatter`.
        contour_kwargs : dict, optional
            Keyword arguments passed to :func:`matplotlib.pyplot.contour`.
        **kwargs
            Keyword arguments passed to :func:`matplotlib.pyplot.imshow`.

        Returns
        -------
        matplotlib.axes.Axes
            Axes containing resulting image
"""


# Decorator for consistently adding plot parameters to docstrings (Model.plot,
# Geometry.plot, Universe.plot, etc.)
def add_plot_params(func):
    func.__doc__ += _PLOT_PARAMS
    return func


def _get_plot_image(plot, cwd):
    from IPython.display import Image

    # Make sure .png file was created
    png_filename = plot.filename if plot.filename is not None else f'plot_{plot.id}'

    # Add file extension if not already present. The C++ code added it
    # automatically if it wasn't present.
    if Path(png_filename).suffix != ".png":
        png_filename += ".png"

    png_file = Path(cwd) / png_filename
    if not png_file.exists():
        raise FileNotFoundError(
            f"Could not find .png image for plot {plot.id}. Your version of "
            "OpenMC may not be built against libpng.")

    return Image(str(png_file))


def voxel_to_vtk(voxel_file: PathLike, output: PathLike = 'plot.vti'):
    """Converts a voxel HDF5 file to a VTK file

    .. versionadded:: 0.14.0

    Parameters
    ----------
    voxel_file : path-like
        Path of the input h5 to convert
    output : path-like
        Path of the output vti file produced

    Returns
    -------
    Path
        Path of the .vti file produced
    """

    # imported vtk only if used as vtk is an option dependency
    import vtk

    _min_version = (2, 0)

    # Read data from voxel file
    with h5py.File(voxel_file, "r") as fh:
        # check version
        version = tuple(fh.attrs["version"])
        if version < _min_version:
            old_version = ".".join(map(str, version))
            min_version = ".".join(map(str, _min_version))
            err_msg = (
                f"This voxel file's version is {old_version}. This function only "
                f" supports voxel files with version {min_version} or higher. "
                "Please generate a new voxel file using a newer version of OpenMC."
            )
            raise ValueError(err_msg)

        dimension = fh.attrs["num_voxels"]
        width = fh.attrs["voxel_width"]
        lower_left = fh.attrs["lower_left"]

        nx, ny, nz = dimension

        grid = vtk.vtkImageData()
        grid.SetDimensions(nx + 1, ny + 1, nz + 1)
        grid.SetOrigin(*lower_left)
        grid.SetSpacing(*width)

        # transpose data from OpenMC ordering (zyx) to VTK ordering (xyz)
        # and flatten to 1-D array
        h5data = fh["data"][...]

    data = vtk.vtkIntArray()
    data.SetName("id")
    # set the array using the h5data array
    data.SetArray(h5data, h5data.size, True)
    # add data to image grid
    grid.GetCellData().AddArray(data)

    writer = vtk.vtkXMLImageDataWriter()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        writer.SetInputData(grid)
    else:
        writer.SetInput(grid)
    output = str(output)
    if not output.endswith(".vti"):
        output += ".vti"
    writer.SetFileName(output)
    writer.Write()

    return output


class PlotBase(IDManagerMixin):
    """
    Parameters
    ----------
    plot_id : int
        Unique identifier for the plot
    name : str
        Name of the plot

    Attributes
    ----------
    id : int
        Unique identifier
    name : str
        Name of the plot
    pixels : Iterable of int
        Number of pixels to use in each direction
    filename : str
        Path to write the plot to
    color_by : {'cell', 'material'}
        Indicate whether the plot should be colored by cell or by material
    background : Iterable of int or str
        Color of the background
    mask_components : Iterable of openmc.Cell or openmc.Material or int
        The cells or materials (or corresponding IDs) to mask
    mask_background : Iterable of int or str
        Color to apply to all cells/materials listed in mask_components
    show_overlaps : bool
        Indicate whether or not overlapping regions are shown
    overlap_color : Iterable of int or str
        Color to apply to overlapping regions
    colors : dict
        Dictionary indicating that certain cells/materials should be
        displayed with a particular color. The keys can be of type
        :class:`~openmc.Cell`, :class:`~openmc.Material`, or int (ID for a
        cell/material).
    level : int
        Universe depth to plot at
    """

    next_id = 1
    used_ids = set()

    def __init__(self, plot_id=None, name=''):
        # Initialize Plot class attributes
        self.id = plot_id
        self.name = name
        self._pixels = [400, 400]
        self._filename = None
        self._color_by = 'cell'
        self._background = None
        self._mask_components = None
        self._mask_background = None
        self._show_overlaps = False
        self._overlap_color = None
        self._colors = {}
        self._level = None

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        cv.check_type('plot name', name, str)
        self._name = name

    @property
    def pixels(self):
        return self._pixels

    @pixels.setter
    def pixels(self, pixels):
        cv.check_type('plot pixels', pixels, Iterable, Integral)
        cv.check_length('plot pixels', pixels, 2, 3)
        for dim in pixels:
            cv.check_greater_than('plot pixels', dim, 0)
        self._pixels = pixels

    @property
    def filename(self):
        return self._filename

    @filename.setter
    def filename(self, filename):
        cv.check_type('filename', filename, str)
        self._filename = filename

    @property
    def color_by(self):
        return self._color_by

    @color_by.setter
    def color_by(self, color_by):
        cv.check_value('plot color_by', color_by, ['cell', 'material'])
        self._color_by = color_by

    @property
    def background(self):
        return self._background

    @background.setter
    def background(self, background):
        self._check_color('plot background', background)
        self._background = background

    @property
    def mask_components(self):
        return self._mask_components

    @mask_components.setter
    def mask_components(self, mask_components):
        cv.check_type('plot mask components', mask_components, Iterable,
                      (openmc.Cell, openmc.Material, Integral))
        self._mask_components = mask_components

    @property
    def mask_background(self):
        return self._mask_background

    @mask_background.setter
    def mask_background(self, mask_background):
        self._check_color('plot mask background', mask_background)
        self._mask_background = mask_background

    @property
    def show_overlaps(self):
        return self._show_overlaps

    @show_overlaps.setter
    def show_overlaps(self, show_overlaps):
        cv.check_type(f'Show overlaps flag for Plot ID="{self.id}"',
                      show_overlaps, bool)
        self._show_overlaps = show_overlaps

    @property
    def overlap_color(self):
        return self._overlap_color

    @overlap_color.setter
    def overlap_color(self, overlap_color):
        self._check_color('plot overlap color', overlap_color)
        self._overlap_color = overlap_color

    @property
    def colors(self):
        return self._colors

    @colors.setter
    def colors(self, colors):
        cv.check_type('plot colors', colors, Mapping)
        for key, value in colors.items():
            cv.check_type('plot color key', key,
                          (openmc.Cell, openmc.Material, Integral))
            self._check_color('plot color value', value)
        self._colors = colors

    @property
    def level(self):
        return self._level

    @level.setter
    def level(self, plot_level):
        cv.check_type('plot level', plot_level, Integral)
        cv.check_greater_than('plot level', plot_level, 0, equality=True)
        self._level = plot_level

    @staticmethod
    def _check_color(err_string, color):
        cv.check_type(err_string, color, Iterable)
        if isinstance(color, str):
            if color.lower() not in _SVG_COLORS:
                raise ValueError(f"'{color}' is not a valid color.")
        else:
            cv.check_length(err_string, color, 3)
            for rgb in color:
                cv.check_type(err_string, rgb, Real)
                cv.check_greater_than('RGB component', rgb, 0, True)
                cv.check_less_than('RGB component', rgb, 256)

    # Helper function that returns the domain ID given either a
    # Cell/Material object or the domain ID itself
    @staticmethod
    def _get_id(domain):
        return domain if isinstance(domain, Integral) else domain.id

    def colorize(self, geometry, seed=1):
        """Generate a color scheme for each domain in the plot.

        This routine may be used to generate random, reproducible color schemes.
        The colors generated are based upon cell/material IDs in the geometry.

        Parameters
        ----------
        geometry : openmc.Geometry
            The geometry for which the plot is defined
        seed : Integral
            The random number seed used to generate the color scheme

        """

        cv.check_type('geometry', geometry, openmc.Geometry)
        cv.check_type('seed', seed, Integral)
        cv.check_greater_than('seed', seed, 1, equality=True)

        # Get collections of the domains which will be plotted
        if self.color_by == 'material':
            domains = geometry.get_all_materials().values()
        else:
            domains = geometry.get_all_cells().values()

        rng = np.random.RandomState(seed)

        # Generate random colors for each feature
        for domain in domains:
            self.colors[domain] = rng.randint(0, 256, (3,))

    def _colors_to_xml(self, element):
        for domain, color in sorted(self._colors.items(),
                                    key=lambda x: self._get_id(x[0])):
            subelement = ET.SubElement(element, "color")
            subelement.set("id", str(self._get_id(domain)))
            if isinstance(color, str):
                color = _SVG_COLORS[color.lower()]
            subelement.set("rgb", ' '.join(str(x) for x in color))

    def to_xml_element(self):
        """Save common plot attributes to XML element

        Returns
        -------
        element : lxml.etree._Element
            XML element containing plot data

        """

        element = ET.Element("plot")
        element.set("id", str(self._id))
        if len(self._name) > 0:
            element.set("name", str(self.name))
        if self._filename is not None:
            element.set("filename", self._filename)
        element.set("color_by", self._color_by)

        subelement = ET.SubElement(element, "pixels")
        subelement.text = ' '.join(map(str, self._pixels))

        if self._background is not None:
            subelement = ET.SubElement(element, "background")
            color = self._background
            if isinstance(color, str):
                color = _SVG_COLORS[color.lower()]
            subelement.text = ' '.join(str(x) for x in color)

        if self._mask_components is not None:
            subelement = ET.SubElement(element, "mask")
            subelement.set("components", ' '.join(
                str(PlotBase._get_id(d)) for d in self._mask_components))
            color = self._mask_background
            if color is not None:
                if isinstance(color, str):
                    color = _SVG_COLORS[color.lower()]
                subelement.set("background", ' '.join(
                    str(x) for x in color))

        if self._level is not None:
            subelement = ET.SubElement(element, "level")
            subelement.text = str(self._level)

        return element


class Plot(PlotBase):
    """Definition of a finite region of space to be plotted.

    OpenMC is capable of generating two-dimensional slice plots, or
    three-dimensional voxel or projection plots. Colors that are used in plots can be given as
    RGB tuples, e.g. (255, 255, 255) would be white, or by a string indicating a
    valid `SVG color <https://www.w3.org/TR/SVG11/types.html#ColorKeywords>`_.

    Parameters
    ----------
    plot_id : int
        Unique identifier for the plot
    name : str
        Name of the plot

    Attributes
    ----------
    id : int
        Unique identifier
    name : str
        Name of the plot
    pixels : Iterable of int
        Number of pixels to use in each direction
    filename : str
        Path to write the plot to
    color_by : {'cell', 'material'}
        Indicate whether the plot should be colored by cell or by material
    background : Iterable of int or str
        Color of the background
    mask_components : Iterable of openmc.Cell or openmc.Material or int
        The cells or materials (or corresponding IDs) to mask
    mask_background : Iterable of int or str
        Color to apply to all cells/materials listed in mask_components
    show_overlaps : bool
        Indicate whether or not overlapping regions are shown
    overlap_color : Iterable of int or str
        Color to apply to overlapping regions
    colors : dict
        Dictionary indicating that certain cells/materials should be
        displayed with a particular color. The keys can be of type
        :class:`~openmc.Cell`, :class:`~openmc.Material`, or int (ID for a
        cell/material).
    level : int
        Universe depth to plot at
    width : Iterable of float
        Width of the plot in each basis direction
    origin : tuple or list of ndarray
        Origin (center) of the plot
    type : {'slice', 'voxel'}
        The type of the plot
    basis : {'xy', 'xz', 'yz'}
        The basis directions for the plot
    meshlines : dict
        Dictionary defining type, id, linewidth and color of a mesh to be
        plotted on top of a plot

    """

    def __init__(self, plot_id=None, name=''):
        super().__init__(plot_id, name)
        self._width = [4.0, 4.0]
        self._origin = [0., 0., 0.]
        self._type = 'slice'
        self._basis = 'xy'
        self._meshlines = None

    @property
    def width(self):
        return self._width

    @width.setter
    def width(self, width):
        cv.check_type('plot width', width, Iterable, Real)
        cv.check_length('plot width', width, 2, 3)
        self._width = width

    @property
    def origin(self):
        return self._origin

    @origin.setter
    def origin(self, origin):
        cv.check_type('plot origin', origin, Iterable, Real)
        cv.check_length('plot origin', origin, 3)
        self._origin = origin

    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, plottype):
        cv.check_value('plot type', plottype, ['slice', 'voxel'])
        self._type = plottype

    @property
    def basis(self):
        return self._basis

    @basis.setter
    def basis(self, basis):
        cv.check_value('plot basis', basis, _BASES)
        self._basis = basis

    @property
    def meshlines(self):
        return self._meshlines

    @meshlines.setter
    def meshlines(self, meshlines):
        cv.check_type('plot meshlines', meshlines, dict)
        if 'type' not in meshlines:
            msg = f'Unable to set the meshlines to "{meshlines}" which ' \
                'does not have a "type" key'
            raise ValueError(msg)

        elif meshlines['type'] not in ['tally', 'entropy', 'ufs', 'cmfd']:
            msg = f"Unable to set the meshlines with type \"{meshlines['type']}\""
            raise ValueError(msg)

        if 'id' in meshlines:
            cv.check_type('plot meshlines id', meshlines['id'], Integral)
            cv.check_greater_than('plot meshlines id', meshlines['id'], 0,
                                  equality=True)

        if 'linewidth' in meshlines:
            cv.check_type('plot mesh linewidth',
                          meshlines['linewidth'], Integral)
            cv.check_greater_than('plot mesh linewidth', meshlines['linewidth'],
                                  0, equality=True)

        if 'color' in meshlines:
            self._check_color('plot meshlines color', meshlines['color'])

        self._meshlines = meshlines

    def __repr__(self):
        string = 'Plot\n'
        string += '{: <16}=\t{}\n'.format('\tID', self._id)
        string += '{: <16}=\t{}\n'.format('\tName', self._name)
        string += '{: <16}=\t{}\n'.format('\tFilename', self._filename)
        string += '{: <16}=\t{}\n'.format('\tType', self._type)
        string += '{: <16}=\t{}\n'.format('\tBasis', self._basis)
        string += '{: <16}=\t{}\n'.format('\tWidth', self._width)
        string += '{: <16}=\t{}\n'.format('\tOrigin', self._origin)
        string += '{: <16}=\t{}\n'.format('\tPixels', self._pixels)
        string += '{: <16}=\t{}\n'.format('\tColor by', self._color_by)
        string += '{: <16}=\t{}\n'.format('\tBackground', self._background)
        string += '{: <16}=\t{}\n'.format('\tMask components',
                                          self._mask_components)
        string += '{: <16}=\t{}\n'.format('\tMask background',
                                          self._mask_background)
        string += '{: <16}=\t{}\n'.format('\tOverlap Color',
                                          self._overlap_color)
        string += '{: <16}=\t{}\n'.format('\tColors', self._colors)
        string += '{: <16}=\t{}\n'.format('\tLevel', self._level)
        string += '{: <16}=\t{}\n'.format('\tMeshlines', self._meshlines)
        return string

    @classmethod
    def from_geometry(cls, geometry, basis='xy', slice_coord=0.):
        """Return plot that encompasses a geometry.

        Parameters
        ----------
        geometry : openmc.Geometry
            The geometry to base the plot off of
        basis : {'xy', 'xz', 'yz'}
            The basis directions for the plot
        slice_coord : float
            The level at which the slice plot should be plotted. For example, if
            the basis is 'xy', this would indicate the z value used in the
            origin.

        """
        cv.check_type('geometry', geometry, openmc.Geometry)
        cv.check_value('basis', basis, _BASES)

        # Decide which axes to keep
        if basis == 'xy':
            pick_index = (0, 1)
            slice_index = 2
        elif basis == 'yz':
            pick_index = (1, 2)
            slice_index = 0
        elif basis == 'xz':
            pick_index = (0, 2)
            slice_index = 1

        # Get lower-left and upper-right coordinates for desired axes
        lower_left, upper_right = geometry.bounding_box
        lower_left = lower_left[np.array(pick_index)]
        upper_right = upper_right[np.array(pick_index)]

        if np.any(np.isinf((lower_left, upper_right))):
            raise ValueError('The geometry does not appear to be bounded '
                             f'in the {basis} plane.')

        plot = cls()
        plot.origin = np.insert((lower_left + upper_right)/2,
                                slice_index, slice_coord)
        plot.width = upper_right - lower_left
        plot.basis = basis
        return plot

    def highlight_domains(self, geometry, domains, seed=1,
                          alpha=0.5, background='gray'):
        """Use alpha compositing to highlight one or more domains in the plot.

        This routine generates a color scheme and applies alpha compositing to
        make all domains except the highlighted ones appear partially
        transparent.

        Parameters
        ----------
        geometry : openmc.Geometry
            The geometry for which the plot is defined
        domains : Iterable of openmc.Cell or openmc.Material
            A collection of the domain IDs to highlight in the plot
        seed : int
            The random number seed used to generate the color scheme
        alpha : float
            The value between 0 and 1 to apply in alpha compositing
        background : 3-tuple of int or str
            The background color to apply in alpha compositing

        """

        cv.check_type('domains', domains, Iterable,
                      (openmc.Cell, openmc.Material))
        cv.check_type('alpha', alpha, Real)
        cv.check_greater_than('alpha', alpha, 0., equality=True)
        cv.check_less_than('alpha', alpha, 1., equality=True)
        cv.check_type('background', background, Iterable)

        # Get a background (R,G,B) tuple to apply in alpha compositing
        if isinstance(background, str):
            if background.lower() not in _SVG_COLORS:
                raise ValueError(f"'{background}' is not a valid color.")
            background = _SVG_COLORS[background.lower()]

        # Generate a color scheme
        self.colorize(geometry, seed)

        # Apply alpha compositing to the colors for all domains
        # other than those the user wishes to highlight
        for domain, color in self.colors.items():
            if domain not in domains:
                if isinstance(color, str):
                    color = _SVG_COLORS[color.lower()]
                r, g, b = color
                r = int(((1-alpha) * background[0]) + (alpha * r))
                g = int(((1-alpha) * background[1]) + (alpha * g))
                b = int(((1-alpha) * background[2]) + (alpha * b))
                self._colors[domain] = (r, g, b)

    def to_xml_element(self):
        """Return XML representation of the slice/voxel plot

        Returns
        -------
        element : lxml.etree._Element
            XML element containing plot data

        """

        element = super().to_xml_element()
        element.set("type", self._type)

        if self._type == 'slice':
            element.set("basis", self._basis)

        subelement = ET.SubElement(element, "origin")
        subelement.text = ' '.join(map(str, self._origin))

        subelement = ET.SubElement(element, "width")
        subelement.text = ' '.join(map(str, self._width))

        if self._colors:
            self._colors_to_xml(element)

        if self._show_overlaps:
            subelement = ET.SubElement(element, "show_overlaps")
            subelement.text = "true"

            if self._overlap_color is not None:
                color = self._overlap_color
                if isinstance(color, str):
                    color = _SVG_COLORS[color.lower()]
                subelement = ET.SubElement(element, "overlap_color")
                subelement.text = ' '.join(str(x) for x in color)

        if self._meshlines is not None:
            subelement = ET.SubElement(element, "meshlines")
            subelement.set("meshtype", self._meshlines['type'])
            if 'id' in self._meshlines:
                subelement.set("id", str(self._meshlines['id']))
            if 'linewidth' in self._meshlines:
                subelement.set("linewidth", str(self._meshlines['linewidth']))
            if 'color' in self._meshlines:
                subelement.set("color", ' '.join(map(
                    str, self._meshlines['color'])))

        return element

    @classmethod
    def from_xml_element(cls, elem):
        """Generate plot object from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element

        Returns
        -------
        openmc.Plot
            Plot object

        """
        plot_id = int(elem.get("id"))
        name = get_text(elem, 'name', '')
        plot = cls(plot_id, name)
        if "filename" in elem.keys():
            plot.filename = elem.get("filename")
        plot.color_by = elem.get("color_by")
        plot.type = elem.get("type")
        if plot.type == 'slice':
            plot.basis = elem.get("basis")

        plot.origin = get_elem_tuple(elem, "origin", float)
        plot.width = get_elem_tuple(elem, "width", float)
        plot.pixels = get_elem_tuple(elem, "pixels")
        plot._background = get_elem_tuple(elem, "background")

        # Set plot colors
        colors = {}
        for color_elem in elem.findall("color"):
            uid = int(color_elem.get("id"))
            colors[uid] = tuple([int(x)
                                for x in color_elem.get("rgb").split()])
        plot.colors = colors

        # Set masking information
        mask_elem = elem.find("mask")
        if mask_elem is not None:
            plot.mask_components = [
                int(x) for x in mask_elem.get("components").split()]
            background = mask_elem.get("background")
            if background is not None:
                plot.mask_background = tuple(
                    [int(x) for x in background.split()])

        # show overlaps
        overlap_elem = elem.find("show_overlaps")
        if overlap_elem is not None:
            plot.show_overlaps = (overlap_elem.text in ('true', '1'))
        overlap_color = get_elem_tuple(elem, "overlap_color")
        if overlap_color is not None:
            plot.overlap_color = overlap_color

        # Set universe level
        level = elem.find("level")
        if level is not None:
            plot.level = int(level.text)

        # Set meshlines
        mesh_elem = elem.find("meshlines")
        if mesh_elem is not None:
            meshlines = {'type': mesh_elem.get('meshtype')}
            if 'id' in mesh_elem.keys():
                meshlines['id'] = int(mesh_elem.get('id'))
            if 'linewidth' in mesh_elem.keys():
                meshlines['linewidth'] = int(mesh_elem.get('linewidth'))
            if 'color' in mesh_elem.keys():
                meshlines['color'] = tuple(
                    [int(x) for x in mesh_elem.get('color').split()]
                )
            plot.meshlines = meshlines

        return plot

    def to_ipython_image(self, openmc_exec='openmc', cwd='.'):
        """Render plot as an image

        This method runs OpenMC in plotting mode to produce a .png file.

        .. versionchanged:: 0.13.0
            The *convert_exec* argument was removed since OpenMC now produces
            .png images directly.

        Parameters
        ----------
        openmc_exec : str
            Path to OpenMC executable
        cwd : str, optional
            Path to working directory to run in

        Returns
        -------
        IPython.display.Image
            Image generated

        """
        # Create plots.xml
        Plots([self]).export_to_xml(cwd)

        # Run OpenMC in geometry plotting mode
        openmc.plot_geometry(False, openmc_exec, cwd)

        # Return produced image
        return _get_plot_image(self, cwd)

    def to_vtk(self, output: PathLike | None = None,
               openmc_exec: str = 'openmc', cwd: str = '.'):
        """Render plot as an voxel image

        This method runs OpenMC in plotting mode to produce a .vti file.

        .. versionadded:: 0.14.0

        Parameters
        ----------
        output : path-like
            Path of the output .vti file produced
        openmc_exec : str
            Path to OpenMC executable
        cwd : str, optional
            Path to working directory to run in

        Returns
        -------
        Path
            Path of the .vti file produced

        """
        if self.type != 'voxel':
            raise ValueError(
                'Generating a VTK file only works for voxel plots')

        # Create plots.xml
        Plots([self]).export_to_xml(cwd)

        # Run OpenMC in geometry plotting mode and produces a h5 file
        openmc.plot_geometry(False, openmc_exec, cwd)

        h5_voxel_filename = self.filename if self.filename is not None else f'plot_{self.id}'

        # Add file extension if not already present
        if Path(h5_voxel_filename).suffix != ".h5":
            h5_voxel_filename += ".h5"

        h5_voxel_file = Path(cwd) / h5_voxel_filename
        if output is None:
            output = h5_voxel_file.with_suffix('.vti')

        return voxel_to_vtk(h5_voxel_file, output)


class RayTracePlot(PlotBase):
    """Definition of a camera's view of OpenMC geometry

    The camera projection may either by orthographic or perspective. Perspective
    projections are more similar to a pinhole camera, and orthographic
    projections preserve parallel lines and distances.

    This is an abstract base class that :class:`WireframeRayTracePlot` and
    :class:`SolidRayTracePlot` finish the implementation of.

    .. versionadded:: 0.15.1

    Parameters
    ----------
    plot_id : int
        Unique identifier for the plot
    name : str
        Name of the plot

    Attributes
    ----------
    horizontal_field_of_view : float
        Field of view horizontally, in units of degrees, defaults to 70.
    camera_position : tuple or list of ndarray
        Position of the camera in 3D space. Defaults to (1, 0, 0).
    look_at : tuple or list of ndarray
        The center of the camera's image points to this place in 3D space.
        Set to (0, 0, 0) by default.
    up : tuple or list of ndarray
        Which way is up for the camera. Must not be parallel to the
        line between look_at and camera_position. Set to (0, 0, 1) by default.
    orthographic_width : float
        If set to a nonzero value, an orthographic projection is used.
        All rays traced from the orthographic pixel array travel in the
        same direction. The width of the starting array must be specified,
        unlike with the default perspective projection. The height of the
        array is deduced from the ratio of pixel dimensions for the image.
        Defaults to zero, i.e. using perspective projection.
    """

    def __init__(self, plot_id=None, name=''):
        # Initialize Plot class attributes
        super().__init__(plot_id, name)
        self._horizontal_field_of_view = 70.0
        self._camera_position = (1.0, 0.0, 0.0)
        self._look_at = (0.0, 0.0, 0.0)
        self._up = (0.0, 0.0, 1.0)
        self._orthographic_width = 0.0

    @property
    def horizontal_field_of_view(self):
        return self._horizontal_field_of_view

    @horizontal_field_of_view.setter
    def horizontal_field_of_view(self, horizontal_field_of_view):
        cv.check_type('plot horizontal field of view', horizontal_field_of_view,
                      Real)
        assert horizontal_field_of_view > 0.0
        assert horizontal_field_of_view < 180.0
        self._horizontal_field_of_view = horizontal_field_of_view

    @property
    def camera_position(self):
        return self._camera_position

    @camera_position.setter
    def camera_position(self, camera_position):
        cv.check_type('plot camera position', camera_position, Iterable, Real)
        cv.check_length('plot camera position', camera_position, 3)
        self._camera_position = camera_position

    @property
    def look_at(self):
        return self._look_at

    @look_at.setter
    def look_at(self, look_at):
        cv.check_type('plot look at', look_at, Iterable, Real)
        cv.check_length('plot look at', look_at, 3)
        self._look_at = look_at

    @property
    def up(self):
        return self._up

    @up.setter
    def up(self, up):
        cv.check_type('plot up', up, Iterable, Real)
        cv.check_length('plot up', up, 3)
        self._up = up

    @property
    def orthographic_width(self):
        return self._orthographic_width

    @orthographic_width.setter
    def orthographic_width(self, orthographic_width):
        cv.check_type('plot orthographic width', orthographic_width, Real)
        assert orthographic_width >= 0.0
        self._orthographic_width = orthographic_width

    def _check_domains_consistent_with_color_by(self, domains):
        """Check domains are the same as the type we are coloring by"""
        for region in domains:
            # if an integer is passed, we have to assume it was a valid ID
            if isinstance(region, int):
                continue

            if self._color_by == 'material':
                if not isinstance(region, openmc.Material):
                    raise Exception('Domain list must be materials if '
                                    'color_by=material')
            else:
                if not isinstance(region, openmc.Cell):
                    raise Exception('Domain list must be cells if '
                                    'color_by=cell')

    def to_xml_element(self):
        """Return XML representation of the ray trace plot

        Returns
        -------
        element : lxml.etree._Element
            XML element containing plot data

        """

        element = super().to_xml_element()
        element.set("id", str(self._id))

        subelement = ET.SubElement(element, "camera_position")
        subelement.text = ' '.join(map(str, self._camera_position))

        subelement = ET.SubElement(element, "look_at")
        subelement.text = ' '.join(map(str, self._look_at))

        subelement = ET.SubElement(element, "horizontal_field_of_view")
        subelement.text = str(self._horizontal_field_of_view)

        # do not need to write if orthographic_width == 0.0
        if self._orthographic_width > 0.0:
            subelement = ET.SubElement(element, "orthographic_width")
            subelement.text = str(self._orthographic_width)

        return element

    def __repr__(self):
        string = ''
        string += '{: <16}=\t{}\n'.format('\tID', self._id)
        string += '{: <16}=\t{}\n'.format('\tName', self._name)
        string += '{: <16}=\t{}\n'.format('\tFilename', self._filename)
        string += '{: <16}=\t{}\n'.format('\tHorizontal FOV',
                                          self._horizontal_field_of_view)
        string += '{: <16}=\t{}\n'.format('\tOrthographic width',
                                          self._orthographic_width)
        string += '{: <16}=\t{}\n'.format('\tCamera position',
                                          self._camera_position)
        string += '{: <16}=\t{}\n'.format('\tLook at', self._look_at)
        string += '{: <16}=\t{}\n'.format('\tUp', self._up)
        string += '{: <16}=\t{}\n'.format('\tPixels', self._pixels)
        string += '{: <16}=\t{}\n'.format('\tColor by', self._color_by)
        string += '{: <16}=\t{}\n'.format('\tBackground', self._background)
        string += '{: <16}=\t{}\n'.format('\tColors', self._colors)
        string += '{: <16}=\t{}\n'.format('\tLevel', self._level)
        return string

    def _read_xml_attributes(self, elem):
        """Helper function called by from_xml_element
        of child classes. These are common vaues to be
        read by any ray traced plot.

        Returns
        -------
        None
        """

        if "filename" in elem.keys():
            self.filename = elem.get("filename")
        self.color_by = elem.get("color_by")

        horizontal_fov = elem.find("horizontal_field_of_view")
        if horizontal_fov is not None:
            self.horizontal_field_of_view = float(horizontal_fov.text)

        if (tmp := elem.find("orthographic_width")) is not None:
            self.orthographic_width = float(tmp)

        self.pixels = get_elem_tuple(elem, "pixels")
        self.camera_position = get_elem_tuple(elem, "camera_position", float)
        self.look_at = get_elem_tuple(elem, "look_at", float)

        if elem.find("background") is not None:
            self.background = get_elem_tuple(elem, "background")

        # Set masking information
        if (mask_elem := elem.find("mask")) is not None:
            mask_components = [int(x)
                               for x in mask_elem.get("components").split()]
            # TODO: set mask components(needs geometry information)
            background = mask_elem.get("background")
            if background is not None:
                self.mask_background = tuple(
                    [int(x) for x in background.split()])

        # Set universe level
        level = elem.find("level")
        if level is not None:
            self.level = int(level.text)


class WireframeRayTracePlot(RayTracePlot):
    """Plots wireframes of geometry with volume rendered colors

    Colors are defined in the same manner as the Plot class, but with the
    addition of a coloring parameter resembling a macroscopic cross section in
    units of inverse centimeters. The volume rendering technique is used to
    color regions of the model. An infinite cross section denotes a fully opaque
    region, and zero represents a transparent region which will expose the color
    of the regions behind it.

    .. versionchanged:: 0.15.1
        Renamed from ProjectionPlot to WireframeRayTracePlot

    Parameters
    ----------
    plot_id : int
        Unique identifier for the plot
    name : str
        Name of the plot

    Attributes
    ----------
    id : int
        Unique identifier
    name : str
        Name of the plot
    pixels : Iterable of int
        Number of pixels to use in each direction
    filename : str
        Path to write the plot to
    color_by : {'cell', 'material'}
        Indicate whether the plot should be colored by cell or by material
    background : Iterable of int or str
        Color of the background
    mask_components : Iterable of openmc.Cell or openmc.Material or int
        The cells or materials (or corresponding IDs) to mask
    mask_background : Iterable of int or str
        Color to apply to all cells/materials listed in mask_components
    show_overlaps : bool
        Indicate whether or not overlapping regions are shown
    overlap_color : Iterable of int or str
        Color to apply to overlapping regions
    colors : dict
        Dictionary indicating that certain cells/materials should be
        displayed with a particular color. The keys can be of type
        :class:`~openmc.Cell`, :class:`~openmc.Material`, or int (ID for a
        cell/material).
    level : int
        Universe depth to plot at
    horizontal_field_of_view : float
        Field of view horizontally, in units of degrees, defaults to 70.
    camera_position : tuple or list of ndarray
        Position of the camera in 3D space. Defaults to (1, 0, 0).
    look_at : tuple or list of ndarray
        The center of the camera's image points to this place in 3D space.
        Set to (0, 0, 0) by default.
    up : tuple or list of ndarray
        Which way is up for the camera. Must not be parallel to the
        line between look_at and camera_position. Set to (0, 0, 1) by default.
    orthographic_width : float
        If set to a nonzero value, an orthographic projection is used.
        All rays traced from the orthographic pixel array travel in the
        same direction. The width of the starting array must be specified,
        unlike with the default perspective projection. The height of the
        array is deduced from the ratio of pixel dimensions for the image.
        Defaults to zero, i.e. using perspective projection.
    wireframe_thickness : int
        Line thickness employed for drawing wireframes around cells or material
        regions. Can be set to zero for no wireframes at all. Defaults to one
        pixel.
    wireframe_color : tuple of ints
        RGB color of the wireframe lines. Defaults to black.
    wireframe_domains : iterable of either Material or Cells
        If provided, the wireframe is only drawn around these. If color_by is by
        material, it must be a list of materials, else cells.
    xs : dict
        A mapping from cell/material IDs to floats. The floating point values
        are macroscopic cross sections influencing the volume rendering opacity
        of each geometric region. Zero corresponds to perfect transparency, and
        infinity equivalent to opaque. These must be set by the user, but
        default values can be obtained using the :meth:`set_transparent` method.
    """

    def __init__(self, plot_id=None, name=''):
        super().__init__(plot_id, name)
        self._wireframe_thickness = 1
        self._wireframe_color = _SVG_COLORS['black']
        self._wireframe_domains = []
        self._xs = {}

    @property
    def wireframe_thickness(self):
        return self._wireframe_thickness

    @wireframe_thickness.setter
    def wireframe_thickness(self, wireframe_thickness):
        cv.check_type('plot wireframe thickness',
                      wireframe_thickness, Integral)
        assert wireframe_thickness >= 0
        self._wireframe_thickness = wireframe_thickness

    @property
    def wireframe_color(self):
        return self._wireframe_color

    @wireframe_color.setter
    def wireframe_color(self, wireframe_color):
        self._check_color('plot wireframe color', wireframe_color)
        self._wireframe_color = wireframe_color

    @property
    def wireframe_domains(self):
        return self._wireframe_domains

    @wireframe_domains.setter
    def wireframe_domains(self, wireframe_domains):
        self._wireframe_domains = wireframe_domains

    @property
    def xs(self):
        return self._xs

    @xs.setter
    def xs(self, xs):
        cv.check_type('plot xs', xs, Mapping)
        for key, value in xs.items():
            cv.check_type('plot xs key', key, (openmc.Cell, openmc.Material))
            cv.check_type('plot xs value', value, Real)
            assert value >= 0.0
        self._xs = xs

    def set_transparent(self, geometry):
        """Sets all volume rendering XS to zero for the model

        Parameters
        ----------
        geometry : openmc.Geometry
            The geometry for which the plot is defined
        """

        cv.check_type('geometry', geometry, openmc.Geometry)

        # Get collections of the domains which will be plotted
        if self.color_by == 'material':
            domains = geometry.get_all_materials().values()
        else:
            domains = geometry.get_all_cells().values()

        # Generate random colors for each feature
        for domain in domains:
            self.xs[domain] = 0.0

    def __repr__(self):
        string = 'Wireframe Ray-traced Plot\n'
        string += super().__repr__()
        string += '{: <16}=\t{}\n'.format('\tWireframe thickness',
                                          self._wireframe_thickness)
        string += '{: <16}=\t{}\n'.format('\tWireframe color',
                                          self._wireframe_color)
        string += '{: <16}=\t{}\n'.format('\tWireframe domains',
                                          self._wireframe_domains)
        string += '{: <16}=\t{}\n'.format('\tTransparencies', self._xs)
        return string

    def to_xml_element(self):
        """Return XML representation of the projection plot

        Returns
        -------
        element : lxml.etree._Element
            XML element containing plot data

        """
        element = super().to_xml_element()
        element.set("type", "wireframe_raytrace")

        subelement = ET.SubElement(element, "wireframe_thickness")
        subelement.text = str(self._wireframe_thickness)

        subelement = ET.SubElement(element, "wireframe_color")
        color = self._wireframe_color
        if isinstance(color, str):
            color = _SVG_COLORS[color.lower()]
        subelement.text = ' '.join(str(x) for x in color)

        self._check_domains_consistent_with_color_by(self.wireframe_domains)

        if self._wireframe_domains:
            id_list = [x.id for x in self._wireframe_domains]
            subelement = ET.SubElement(element, "wireframe_ids")
            subelement.text = ' '.join([str(x) for x in id_list])

        # note that this differs from the slice plot colors
        # in that "xs" must also be specified
        if self._colors:
            for domain, color in sorted(self._colors.items(),
                                        key=lambda x: x[0].id):
                subelement = ET.SubElement(element, "color")
                subelement.set("id", str(domain.id))
                if isinstance(color, str):
                    color = _SVG_COLORS[color.lower()]
                subelement.set("rgb", ' '.join(str(x) for x in color))
                subelement.set("xs", str(self._xs[domain]))

        return element

    @classmethod
    def from_xml_element(cls, elem):
        """Generate plot object from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element

        Returns
        -------
        openmc.WireframeRayTracePlot
            WireframeRayTracePlot object

        """

        plot_id = int(elem.get("id"))
        plot_name = get_text(elem, 'name', '')
        plot = cls(plot_id, plot_name)
        plot.type = "wireframe_raytrace"

        plot._read_xml_attributes(elem)

        # Attempt to get wireframe thickness.May not be present
        wireframe_thickness = elem.find("wireframe_thickness")
        if wireframe_thickness is not None:
            plot.wireframe_thickness = int(wireframe_thickness.text)
        wireframe_color = elem.get("wireframe_color")
        if wireframe_color:
            plot.wireframe_color = [int(item) for item in wireframe_color]

        # Set plot colors
        for color_elem in elem.findall("color"):
            uid = int(color_elem.get("id"))
            plot.colors[uid] = tuple(int(i)
                                     for i in get_text(color_elem, 'rgb').split())
            plot.xs[uid] = float(color_elem.get("xs"))

        return plot


class SolidRayTracePlot(RayTracePlot):
    """Phong shading-based rendering of an OpenMC geometry

    This class defines a plot that uses Phong shading to enhance the
    visualization of an OpenMC geometry by incorporating diffuse lighting and
    configurable opacity for certain regions. It extends :class:`RayTracePlot`
    by adding parameters related to lighting and transparency.

    .. versionadded:: 0.15.1

    Parameters
    ----------
    plot_id : int, optional
        Unique identifier for the plot
    name : str, optional
        Name of the plot

    Attributes
    ----------
    id : int
        Unique identifier
    name : str
        Name of the plot
    pixels : Iterable of int
        Number of pixels to use in each direction
    filename : str
        Path to write the plot to
    color_by : {'cell', 'material'}
        Indicate whether the plot should be colored by cell or by material
    overlap_color : Iterable of int or str
        Color to apply to overlapping regions
    colors : dict
        Dictionary indicating that certain cells/materials should be
        displayed with a particular color. The keys can be of type
        :class:`~openmc.Cell`, :class:`~openmc.Material`, or int (ID for a
        cell/material).
    horizontal_field_of_view : float
        Field of view horizontally, in units of degrees, defaults to 70.
    camera_position : tuple or list of ndarray
        Position of the camera in 3D space. Defaults to (1, 0, 0).
    look_at : tuple or list of ndarray
        The center of the camera's image points to this place in 3D space.
        Set to (0, 0, 0) by default.
    up : tuple or list of ndarray
        Which way is up for the camera. Must not be parallel to the
        line between look_at and camera_position. Set to (0, 0, 1) by default.
    orthographic_width : float
        If set to a nonzero value, an orthographic projection is used.
        All rays traced from the orthographic pixel array travel in the
        same direction. The width of the starting array must be specified,
        unlike with the default perspective projection. The height of the
        array is deduced from the ratio of pixel dimensions for the image.
        Defaults to zero, i.e. using perspective projection.
    light_position : tuple or list of float
        Position of the light source in 3D space. Defaults to None, which places
        the light at the camera position.
    diffuse_fraction : float
        Fraction of lighting that is diffuse (non-directional). Defaults to 0.1.
        Must be between 0 and 1.
    opaque_domains : list
        List of domains (e.g., cells or materials) that should be rendered as
        opaque rather than allowing transparency.
    """

    def __init__(self, plot_id=None, name=''):
        super().__init__(plot_id, name)
        self._light_position = None
        self._diffuse_fraction = 0.1
        self._opaque_domains = []

    @property
    def light_position(self):
        return self._light_position

    @light_position.setter
    def light_position(self, x):
        cv.check_type('plot light position', x, Iterable, Real)
        cv.check_length('plot light position', x, 3)
        self._light_position = x

    @property
    def diffuse_fraction(self):
        return self._diffuse_fraction

    @diffuse_fraction.setter
    def diffuse_fraction(self, x):
        cv.check_type('diffuse fraction', x, Real)
        cv.check_greater_than('diffuse fraction', x, 0.0, equality=True)
        cv.check_less_than('diffuse fraction', x, 1.0, equality=True)
        self._diffuse_fraction = x

    @property
    def opaque_domains(self):
        return self._opaque_domains

    @opaque_domains.setter
    def opaque_domains(self, x):
        # Note that _check_domains_consistent_with_color_by checks
        # the types within later. This is because we don't necessarily
        # know what types are acceptable until the user has set the
        # color_by attribute, too.
        cv.check_type('opaque domains', x, Iterable)
        self._opaque_domains = x

    def __repr__(self):
        string = 'Solid Ray-traced Plot\n'
        string += super().__repr__()
        string += '{: <16}=\t{}\n'.format('\tDiffuse Fraction',
                                          self._diffuse_fraction)
        string += '{: <16}=\t{}\n'.format('\tLight position',
                                          self._light_position)
        string += '{: <16}=\t{}\n'.format('\tOpaque domains',
                                          self._opaque_domains)
        return string

    def to_xml_element(self):
        """Return XML representation of the solid ray-traced plot

        Returns
        -------
        element : lxml.etree._Element
            XML element containing plot data

        """
        element = super().to_xml_element()
        element.set("type", "solid_raytrace")

        # no light position means put it at the camera
        if self._light_position:
            subelement = ET.SubElement(element, "light_position")
            subelement.text = ' '.join(map(str, self._light_position))

        # no diffuse fraction defaults to 0.1
        if self._diffuse_fraction:
            subelement = ET.SubElement(element, "diffuse_fraction")
            subelement.text = str(self._diffuse_fraction)

        self._check_domains_consistent_with_color_by(self.opaque_domains)
        subelement = ET.SubElement(element, "opaque_ids")

        # Extract all IDs, or use the integer value passed in
        # explicitly if that was given
        subelement.text = ' '.join(
            [str(domain) if isinstance(domain, int) else
             str(domain.id) for domain in self._opaque_domains])

        if self._colors:
            self._colors_to_xml(element)

        return element

    def _read_phong_attributes(self, elem):
        """Read attributes specific to the Phong plot from an XML element"""
        if elem.find('light_position') is not None:
            self.light_position = get_elem_tuple(elem, 'light_position', float)

        diffuse_fraction = elem.find('diffuse_fraction')
        if diffuse_fraction is not None:
            self.diffuse_fraction = float(diffuse_fraction.text)

        if elem.find('opaque_ids') is not None:
            self.opaque_domains = list(get_elem_tuple(elem, 'opaque_ids', int))

    @classmethod
    def from_xml_element(cls, elem):
        """Generate plot object from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element

        Returns
        -------
        openmc.WireframeRayTracePlot
            WireframeRayTracePlot object

        """

        plot_id = int(elem.get("id"))
        plot_name = get_text(elem, 'name', '')
        plot = cls(plot_id, plot_name)
        plot.type = "solid_raytrace"

        plot._read_xml_attributes(elem)
        plot._read_phong_attributes(elem)

        # Set plot colors
        for color_elem in elem.findall("color"):
            uid = color_elem.get("id")
            plot.colors[uid] = get_elem_tuple(color_elem, "rgb")

        return plot


class Plots(cv.CheckedList):
    """Collection of Plots used for an OpenMC simulation.

    This class corresponds directly to the plots.xml input file. It can be
    thought of as a normal Python list where each member is inherits from
    :class:`PlotBase`. It behaves like a list as the following example
    demonstrates:

    >>> xz_plot = openmc.Plot()
    >>> big_plot = openmc.Plot()
    >>> small_plot = openmc.Plot()
    >>> p = openmc.Plots((xz_plot, big_plot))
    >>> p.append(small_plot)
    >>> small_plot = p.pop()

    Parameters
    ----------
    plots : Iterable of openmc.PlotBase
        plots to add to the collection

    """

    def __init__(self, plots=None):
        super().__init__(PlotBase, 'plots collection')
        self._plots_file = ET.Element("plots")
        if plots is not None:
            self += plots

    def append(self, plot):
        """Append plot to collection

        Parameters
        ----------
        plot : openmc.PlotBase
            Plot to append

        """
        super().append(plot)

    def insert(self, index, plot):
        """Insert plot before index

        Parameters
        ----------
        index : int
            Index in list
        plot : openmc.Plot
            Plot to insert

        """
        super().insert(index, plot)

    def colorize(self, geometry, seed=1):
        """Generate a consistent color scheme for each domain in each plot.

        This routine may be used to generate random, reproducible color schemes.
        The colors generated are based upon cell/material IDs in the geometry.
        The color schemes will be consistent for all plots in "plots.xml".

        Parameters
        ----------
        geometry : openmc.Geometry
            The geometry for which the plots are defined
        seed : Integral
            The random number seed used to generate the color scheme

        """

        for plot in self:
            plot.colorize(geometry, seed)

    def highlight_domains(self, geometry, domains, seed=1,
                          alpha=0.5, background='gray'):
        """Use alpha compositing to highlight one or more domains in the plot.

        This routine generates a color scheme and applies alpha compositing to
        make all domains except the highlighted ones appear partially
        transparent.

        Parameters
        ----------
        geometry : openmc.Geometry
            The geometry for which the plot is defined
        domains : Iterable of openmc.Cell or openmc.Material
            A collection of the domain IDs to highlight in the plot
        seed : int
            The random number seed used to generate the color scheme
        alpha : float
            The value between 0 and 1 to apply in alpha compositing
        background : 3-tuple of int or str
            The background color to apply in alpha compositing

        """

        for plot in self:
            plot.highlight_domains(geometry, domains, seed, alpha, background)

    def _create_plot_subelements(self):
        for plot in self:
            xml_element = plot.to_xml_element()

            if len(plot.name) > 0:
                self._plots_file.append(ET.Comment(plot.name))

            self._plots_file.append(xml_element)

    def to_xml_element(self):
        """Create a 'plots' element to be written to an XML file.

        Returns
        -------
        element : lxml.etree._Element
            XML element containing all plot elements

        """
        # Reset xml element tree
        self._plots_file.clear()

        self._create_plot_subelements()

        # Clean the indentation in the file to be user-readable
        clean_indentation(self._plots_file)
        # TODO: Remove when support is Python 3.8+
        reorder_attributes(self._plots_file)

        return self._plots_file

    def export_to_xml(self, path='plots.xml'):
        """Export plot specifications to an XML file.

        Parameters
        ----------
        path : str
            Path to file to write. Defaults to 'plots.xml'.

        """
        # Check if path is a directory
        p = Path(path)
        if p.is_dir():
            p /= 'plots.xml'

        self.to_xml_element()
        # Write the XML Tree to the plots.xml file
        tree = ET.ElementTree(self._plots_file)
        tree.write(str(p), xml_declaration=True, encoding='utf-8')

    @classmethod
    def from_xml_element(cls, elem):
        """Generate plots collection from XML file

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element

        Returns
        -------
        openmc.Plots
            Plots collection

        """
        # Generate each plot
        plots = cls()
        for e in elem.findall('plot'):
            plot_type = e.get('type')
            if plot_type == 'wireframe_raytrace':
                plots.append(WireframeRayTracePlot.from_xml_element(e))
            elif plot_type == 'solid_raytrace':
                plots.append(SolidRayTracePlot.from_xml_element(e))
            elif plot_type in ('slice', 'voxel'):
                plots.append(Plot.from_xml_element(e))
            else:
                raise ValueError("Unknown plot type: {}".format(plot_type))
        return plots

    @classmethod
    def from_xml(cls, path='plots.xml'):
        """Generate plots collection from XML file

        Parameters
        ----------
        path : str, optional
            Path to plots XML file

        Returns
        -------
        openmc.Plots
            Plots collection

        """
        parser = ET.XMLParser(huge_tree=True)
        tree = ET.parse(path, parser=parser)
        root = tree.getroot()
        return cls.from_xml_element(root)
