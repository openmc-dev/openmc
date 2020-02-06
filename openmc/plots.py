from collections.abc import Iterable, Mapping
from numbers import Real, Integral
from pathlib import Path
import subprocess
import sys
import warnings
from xml.etree import ElementTree as ET

import numpy as np

import openmc
import openmc.checkvalue as cv
from openmc._xml import clean_indentation
from openmc.mixin import IDManagerMixin


_BASES = ['xy', 'xz', 'yz']

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


class Plot(IDManagerMixin):
    """Definition of a finite region of space to be plotted.

    OpenMC is capable of generating two-dimensional slice plots and
    three-dimensional voxel plots. Colors that are used in plots can be given as
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
    width : Iterable of float
        Width of the plot in each basis direction
    pixels : Iterable of int
        Number of pixels to use in each basis direction
    origin : tuple or list of ndarray
        Origin (center) of the plot
    filename :
        Path to write the plot to
    color_by : {'cell', 'material'}
        Indicate whether the plot should be colored by cell or by material
    type : {'slice', 'voxel'}
        The type of the plot
    basis : {'xy', 'xz', 'yz'}
        The basis directions for the plot
    background : Iterable of int or str
        Color of the background
    mask_components : Iterable of openmc.Cell or openmc.Material
        The cells or materials to plot
    mask_background : Iterable of int or str
        Color to apply to all cells/materials not listed in mask_components
    show_overlaps : bool
        Inidicate whether or not overlapping regions are shown
    overlap_color : Iterable of int or str
        Color to apply to overlapping regions
    colors : dict
        Dictionary indicating that certain cells/materials (keys) should be
        displayed with a particular color.
    level : int
        Universe depth to plot at
    meshlines : dict
        Dictionary defining type, id, linewidth and color of a mesh to be
        plotted on top of a plot

    """

    next_id = 1
    used_ids = set()

    def __init__(self, plot_id=None, name=''):
        # Initialize Plot class attributes
        self.id = plot_id
        self.name = name
        self._width = [4.0, 4.0]
        self._pixels = [400, 400]
        self._origin = [0., 0., 0.]
        self._filename = None
        self._color_by = 'cell'
        self._type = 'slice'
        self._basis = 'xy'
        self._background = None
        self._mask_components = None
        self._mask_background = None
        self._show_overlaps = False
        self._overlap_color = None
        self._colors = {}
        self._level = None
        self._meshlines = None

    @property
    def name(self):
        return self._name

    @property
    def width(self):
        return self._width

    @property
    def pixels(self):
        return self._pixels

    @property
    def origin(self):
        return self._origin

    @property
    def filename(self):
        return self._filename

    @property
    def color_by(self):
        return self._color_by

    @property
    def type(self):
        return self._type

    @property
    def basis(self):
        return self._basis

    @property
    def background(self):
        return self._background

    @property
    def mask_components(self):
        return self._mask_components

    @property
    def mask_background(self):
        return self._mask_background

    @property
    def show_overlaps(self):
        return self._show_overlaps

    @property
    def overlap_color(self):
        return self._overlap_color

    @property
    def colors(self):
        return self._colors

    @property
    def level(self):
        return self._level

    @property
    def meshlines(self):
        return self._meshlines

    @name.setter
    def name(self, name):
        cv.check_type('plot name', name, str)
        self._name = name

    @width.setter
    def width(self, width):
        cv.check_type('plot width', width, Iterable, Real)
        cv.check_length('plot width', width, 2, 3)
        self._width = width

    @origin.setter
    def origin(self, origin):
        cv.check_type('plot origin', origin, Iterable, Real)
        cv.check_length('plot origin', origin, 3)
        self._origin = origin

    @pixels.setter
    def pixels(self, pixels):
        cv.check_type('plot pixels', pixels, Iterable, Integral)
        cv.check_length('plot pixels', pixels, 2, 3)
        for dim in pixels:
            cv.check_greater_than('plot pixels', dim, 0)
        self._pixels = pixels

    @filename.setter
    def filename(self, filename):
        cv.check_type('filename', filename, str)
        self._filename = filename

    @color_by.setter
    def color_by(self, color_by):
        cv.check_value('plot color_by', color_by, ['cell', 'material'])
        self._color_by = color_by

    @type.setter
    def type(self, plottype):
        cv.check_value('plot type', plottype, ['slice', 'voxel'])
        self._type = plottype

    @basis.setter
    def basis(self, basis):
        cv.check_value('plot basis', basis, _BASES)
        self._basis = basis

    @background.setter
    def background(self, background):
        self._check_color('plot background', background)
        self._background = background

    @colors.setter
    def colors(self, colors):
        cv.check_type('plot colors', colors, Mapping)
        for key, value in colors.items():
            cv.check_type('plot color key', key, (openmc.Cell, openmc.Material))
            self._check_color('plot color value', value)
        self._colors = colors

    @mask_components.setter
    def mask_components(self, mask_components):
        cv.check_type('plot mask components', mask_components, Iterable,
                      (openmc.Cell, openmc.Material))
        self._mask_components = mask_components

    @mask_background.setter
    def mask_background(self, mask_background):
        self._check_color('plot mask background', mask_background)
        self._mask_background = mask_background

    @show_overlaps.setter
    def show_overlaps(self, show_overlaps):
        cv.check_type('Show overlaps flag for Plot ID="{}"'.format(self.id),
                      show_overlaps, bool)
        self._show_overlaps = show_overlaps

    @overlap_color.setter
    def overlap_color(self, overlap_color):
        self._check_color('plot overlap color', overlap_color)
        self._overlap_color = overlap_color

    @level.setter
    def level(self, plot_level):
        cv.check_type('plot level', plot_level, Integral)
        cv.check_greater_than('plot level', plot_level, 0, equality=True)
        self._level = plot_level

    @meshlines.setter
    def meshlines(self, meshlines):
        cv.check_type('plot meshlines', meshlines, dict)
        if 'type' not in meshlines:
            msg = 'Unable to set on plot the meshlines "{0}" which ' \
                  'does not have a "type" key'.format(meshlines)
            raise ValueError(msg)

        elif meshlines['type'] not in ['tally', 'entropy', 'ufs', 'cmfd']:
            msg = 'Unable to set the meshlines with ' \
                  'type "{0}"'.format(meshlines['type'])
            raise ValueError(msg)

        if 'id' in meshlines:
            cv.check_type('plot meshlines id', meshlines['id'], Integral)
            cv.check_greater_than('plot meshlines id', meshlines['id'], 0,
                                  equality=True)

        if 'linewidth' in meshlines:
            cv.check_type('plot mesh linewidth', meshlines['linewidth'], Integral)
            cv.check_greater_than('plot mesh linewidth', meshlines['linewidth'],
                                  0, equality=True)

        if 'color' in meshlines:
            self._check_color('plot meshlines color', meshlines['color'])

        self._meshlines = meshlines

    @staticmethod
    def _check_color(err_string, color):
        cv.check_type(err_string, color, Iterable)
        if isinstance(color, str):
            if color.lower() not in _SVG_COLORS:
                raise ValueError("'{}' is not a valid color.".format(color))
        else:
            cv.check_length(err_string, color, 3)
            for rgb in color:
                cv.check_type(err_string, rgb, Real)
                cv.check_greater_than('RGB component', rgb, 0, True)
                cv.check_less_than('RGB component', rgb, 256)

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
        string += '{: <16}=\t{}\n'.format('\Overlap Color',
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
                             'in the {} plane.'.format(basis))

        plot = cls()
        plot.origin = np.insert((lower_left + upper_right)/2,
                                slice_index, slice_coord)
        plot.width = upper_right - lower_left
        return plot

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

        # Set the seed for the random number generator
        np.random.seed(seed)

        # Generate random colors for each feature
        for domain in domains:
            self.colors[domain] = np.random.randint(0, 256, (3,))

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
            The value between 0 and 1 to apply in alpha compisiting
        background : 3-tuple of int or str
            The background color to apply in alpha compisiting

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
                raise ValueError("'{}' is not a valid color.".format(background))
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
        """Return XML representation of the plot

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing plot data

        """

        element = ET.Element("plot")
        element.set("id", str(self._id))
        if self._filename is not None:
            element.set("filename", self._filename)
        element.set("color_by", self._color_by)
        element.set("type", self._type)

        if self._type is 'slice':
            element.set("basis", self._basis)

        subelement = ET.SubElement(element, "origin")
        subelement.text = ' '.join(map(str, self._origin))

        subelement = ET.SubElement(element, "width")
        subelement.text = ' '.join(map(str, self._width))

        subelement = ET.SubElement(element, "pixels")
        subelement.text = ' '.join(map(str, self._pixels))

        if self._background is not None:
            subelement = ET.SubElement(element, "background")
            color = self._background
            if isinstance(color, str):
                color = _SVG_COLORS[color.lower()]
            subelement.text = ' '.join(str(x) for x in color)

        if self._colors:
            for domain, color in sorted(self._colors.items(),
                                        key=lambda x: x[0].id):
                subelement = ET.SubElement(element, "color")
                subelement.set("id", str(domain.id))
                if isinstance(color, str):
                    color = _SVG_COLORS[color.lower()]
                subelement.set("rgb", ' '.join(str(x) for x in color))

        if self._mask_components is not None:
            subelement = ET.SubElement(element, "mask")
            subelement.set("components", ' '.join(
                str(d.id) for d in self._mask_components))
            color = self._mask_background
            if color is not None:
                if isinstance(color, str):
                    color = _SVG_COLORS[color.lower()]
                subelement.set("background", ' '.join(
                    str(x) for x in color))

        if self._show_overlaps:
            subelement = ET.SubElement(element, "show_overlaps")
            subelement.text = "true"

            if self._overlap_color is not None:
                color = self._overlap_color
                if isinstance(color, str):
                    color = _SVG_COLORS[color.lower()]
                subelement = ET.SubElement(element, "overlap_color")
                subelement.text = ' '.join(str(x) for x in color)


        if self._level is not None:
            subelement = ET.SubElement(element, "level")
            subelement.text = str(self._level)

        if self._meshlines is not None:
            subelement = ET.SubElement(element, "meshlines")
            subelement.set("meshtype", self._meshlines['type'])
            if self._meshlines['id'] is not None:
                subelement.set("id", str(self._meshlines['id']))
            if self._meshlines['linewidth'] is not None:
                subelement.set("linewidth", str(self._meshlines['linewidth']))
            if self._meshlines['color'] is not None:
                subelement.set("color", ' '.join(map(
                    str, self._meshlines['color'])))

        return element

    def to_ipython_image(self, openmc_exec='openmc', cwd='.',
                         convert_exec='convert'):
        """Render plot as an image

        This method runs OpenMC in plotting mode to produce a bitmap image which
        is then converted to a .png file and loaded in as an
        :class:`IPython.display.Image` object. As such, it requires that your
        model geometry, materials, and settings have already been exported to
        XML.

        Parameters
        ----------
        openmc_exec : str
            Path to OpenMC executable
        cwd : str, optional
            Path to working directory to run in
        convert_exec : str, optional
            Command that can convert PPM files into PNG files

        Returns
        -------
        IPython.display.Image
            Image generated

        """
        from IPython.display import Image

        # Create plots.xml
        Plots([self]).export_to_xml()

        # Run OpenMC in geometry plotting mode
        openmc.plot_geometry(False, openmc_exec, cwd)

        # Convert to .png
        if self.filename is not None:
            ppm_file = '{}.ppm'.format(self.filename)
        else:
            ppm_file = 'plot_{}.ppm'.format(self.id)
        png_file = ppm_file.replace('.ppm', '.png')
        subprocess.check_call([convert_exec, ppm_file, png_file])

        return Image(png_file)


class Plots(cv.CheckedList):
    """Collection of Plots used for an OpenMC simulation.

    This class corresponds directly to the plots.xml input file. It can be
    thought of as a normal Python list where each member is a :class:`Plot`. It
    behaves like a list as the following example demonstrates:

    >>> xz_plot = openmc.Plot()
    >>> big_plot = openmc.Plot()
    >>> small_plot = openmc.Plot()
    >>> p = openmc.Plots((xz_plot, big_plot))
    >>> p.append(small_plot)
    >>> small_plot = p.pop()

    Parameters
    ----------
    plots : Iterable of openmc.Plot
        Plots to add to the collection

    """

    def __init__(self, plots=None):
        super().__init__(Plot, 'plots collection')
        self._plots_file = ET.Element("plots")
        if plots is not None:
            self += plots

    def append(self, plot):
        """Append plot to collection

        Parameters
        ----------
        plot : openmc.Plot
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
            The value between 0 and 1 to apply in alpha compisiting
        background : 3-tuple of int or str
            The background color to apply in alpha compisiting

        """

        for plot in self:
            plot.highlight_domains(geometry, domains, seed, alpha, background)

    def _create_plot_subelements(self):
        for plot in self:
            xml_element = plot.to_xml_element()

            if len(plot.name) > 0:
                self._plots_file.append(ET.Comment(plot.name))

            self._plots_file.append(xml_element)

    def export_to_xml(self, path='plots.xml'):
        """Export plot specifications to an XML file.

        Parameters
        ----------
        path : str
            Path to file to write. Defaults to 'plots.xml'.

        """
        # Reset xml element tree
        self._plots_file.clear()

        self._create_plot_subelements()

        # Clean the indentation in the file to be user-readable
        clean_indentation(self._plots_file)

        # Check if path is a directory
        p = Path(path)
        if p.is_dir():
            p /= 'plots.xml'

        # Write the XML Tree to the plots.xml file
        tree = ET.ElementTree(self._plots_file)
        tree.write(str(p), xml_declaration=True, encoding='utf-8')
