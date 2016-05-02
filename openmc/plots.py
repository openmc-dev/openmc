from collections import Iterable
from numbers import Real, Integral
from xml.etree import ElementTree as ET
import sys
import warnings

import numpy as np

import openmc
import openmc.checkvalue as cv
from openmc.clean_xml import *

if sys.version_info[0] >= 3:
    basestring = str

# A static variable for auto-generated Plot IDs
AUTO_PLOT_ID = 10000


def reset_auto_plot_id():
    global AUTO_PLOT_ID
    AUTO_PLOT_ID = 10000


BASES = ['xy', 'xz', 'yz']


class Plot(object):
    """Definition of a finite region of space to be plotted, either as a slice plot
    in two dimensions or as a voxel plot in three dimensions.

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
    color : {'cell', 'mat'}
        Indicate whether the plot should be colored by cell or by material
    type : {'slice', 'voxel'}
        The type of the plot
    basis : {'xy', 'xz', 'yz'}
        The basis directions for the plot
    background : tuple or list of ndarray
        Color of the background defined by RGB
    mask_components : Iterable of int
        Unique id numbers of the cells or materials to plot
    mask_background : Iterable of int
        Color to apply to all cells/materials not listed in mask_components
        defined by RGB
    col_spec : dict
        Dictionary indicating that certain cells/materials (keys) should be
        colored with a specific RGB (values)

    """

    def __init__(self, plot_id=None, name=''):
        # Initialize Plot class attributes
        self.id = plot_id
        self.name = name
        self._width = [4.0, 4.0]
        self._pixels = [1000, 1000]
        self._origin = [0., 0., 0.]
        self._filename = 'plot'
        self._color = 'cell'
        self._type = 'slice'
        self._basis = 'xy'
        self._background = [0, 0, 0]
        self._mask_components = None
        self._mask_background = None
        self._col_spec = None

    @property
    def id(self):
        return self._id

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
    def color(self):
        return self._color

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
    def col_spec(self):
        return self._col_spec

    @id.setter
    def id(self, plot_id):
        if plot_id is None:
            global AUTO_PLOT_ID
            self._id = AUTO_PLOT_ID
            AUTO_PLOT_ID += 1
        else:
            cv.check_type('plot ID', plot_id, Integral)
            cv.check_greater_than('plot ID', plot_id, 0, equality=True)
            self._id = plot_id

    @name.setter
    def name(self, name):
        cv.check_type('plot name', name, basestring)
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
        cv.check_type('filename', filename, basestring)
        self._filename = filename

    @color.setter
    def color(self, color):
        cv.check_type('plot color', color, basestring)
        cv.check_value('plot color', color, ['cell', 'mat'])
        self._color = color

    @type.setter
    def type(self, plottype):
        cv.check_type('plot type', plottype, basestring)
        cv.check_value('plot type', plottype, ['slice', 'voxel'])
        self._type = plottype

    @basis.setter
    def basis(self, basis):
        cv.check_type('plot basis', basis, basestring)
        cv.check_value('plot basis', basis, ['xy', 'xz', 'yz'])
        self._basis = basis

    @background.setter
    def background(self, background):
        cv.check_type('plot background', background, Iterable, Integral)
        cv.check_length('plot background', background, 3)
        for rgb in background:
            cv.check_greater_than('plot background',rgb, 0, True)
            cv.check_less_than('plot background', rgb, 256)
        self._background = background

    @col_spec.setter
    def col_spec(self, col_spec):
        cv.check_type('plot col_spec parameter', col_spec, dict, Integral)

        for key in col_spec:
            if key < 0:
                msg = 'Unable to create Plot ID="{0}" with col_spec ID "{1}" ' \
                      'which is less than 0'.format(self._id, key)
                raise ValueError(msg)

            elif not isinstance(col_spec[key], Iterable):
                msg = 'Unable to create Plot ID="{0}" with col_spec RGB values' \
                      ' "{1}" which is not iterable'.format(self._id, col_spec[key])
                raise ValueError(msg)

            elif len(col_spec[key]) != 3:
                msg = 'Unable to create Plot ID="{0}" with col_spec RGB ' \
                      'values of length "{1}" since 3 values must be ' \
                      'input'.format(self._id, len(col_spec[key]))
                raise ValueError(msg)

        self._col_spec = col_spec

    @mask_components.setter
    def mask_components(self, mask_components):
        cv.check_type('plot mask_components', mask_components, Iterable, Integral)
        for component in mask_components:
            cv.check_greater_than('plot mask_components', component, 0, True)
        self._mask_components = mask_components

    @mask_background.setter
    def mask_background(self, mask_background):
        cv.check_type('plot mask background', mask_background, Iterable, Integral)
        cv.check_length('plot mask background', mask_background, 3)
        for rgb in mask_background:
            cv.check_greater_than('plot mask background', rgb, 0, True)
            cv.check_less_than('plot mask background', rgb, 256)
        self._mask_background = mask_background

    def __repr__(self):
        string = 'Plot\n'
        string += '{0: <16}{1}{2}\n'.format('\tID', '=\t', self._id)
        string += '{0: <16}{1}{2}\n'.format('\tName', '=\t', self._name)
        string += '{0: <16}{1}{2}\n'.format('\tFilename', '=\t', self._filename)
        string += '{0: <16}{1}{2}\n'.format('\tType', '=\t', self._type)
        string += '{0: <16}{1}{2}\n'.format('\tBasis', '=\t', self._basis)
        string += '{0: <16}{1}{2}\n'.format('\tWidth', '=\t', self._width)
        string += '{0: <16}{1}{2}\n'.format('\tOrigin', '=\t', self._origin)
        string += '{0: <16}{1}{2}\n'.format('\tPixels', '=\t', self._origin)
        string += '{0: <16}{1}{2}\n'.format('\tColor', '=\t', self._color)
        string += '{0: <16}{1}{2}\n'.format('\tMask', '=\t',
                                            self._mask_components)
        string += '{0: <16}{1}{2}\n'.format('\tMask',    '=\t',
                                            self._mask_background)
        string += '{0: <16}{1}{2}\n'.format('\tCol Spec', '=\t', self._col_spec)
        return string

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
        if self.color is 'mat':
            domains = geometry.get_all_materials()
        else:
            domains = geometry.get_all_cells()

        # Set the seed for the random number generator
        np.random.seed(seed)

        # Generate random colors for each feature
        self.col_spec = {}
        for domain in domains:
            r = np.random.randint(0, 256)
            g = np.random.randint(0, 256)
            b = np.random.randint(0, 256)
            self.col_spec[domain] = (r, g, b)

    def highlight_domains(self, geometry, domains, seed=1,
                          alpha=0.5, background='gray'):
        """Use alpha compositing to highlight one or more domains in the plot.

        This routine generates a color scheme and applies alpha compositing
        to make all domains except the highlighted ones appear partially
        transparent.

        Parameters
        ----------
        geometry : openmc.Geometry
            The geometry for which the plot is defined
        domains : Iterable of Integral
            A collection of the domain IDs to highlight in the plot
        seed : Integral
            The random number seed used to generate the color scheme
        alpha : Real in [0,1]
            The value to apply in alpha compisiting
        background : 3-tuple of Integral or 'white' or 'black' or 'gray'
            The background color to apply in alpha compisiting

        """

        cv.check_iterable_type('domains', domains, Integral)
        cv.check_type('alpha', alpha, Real)
        cv.check_greater_than('alpha', alpha, 0., equality=True)
        cv.check_less_than('alpha', alpha, 1., equality=True)

        # Get a background (R,G,B) tuple to apply in alpha compositing
        if isinstance(background, basestring):
            if background == 'white':
                background = (255, 255, 255)
            elif background == 'black':
                background = (0, 0, 0)
            elif background == 'gray':
                background = (160, 160, 160)
            else:
                msg = 'The background "{}" is not defined'.format(background)
                raise ValueError(msg)

        cv.check_iterable_type('background', background, Integral)

        # Generate a color scheme
        self.colorize(geometry, seed)

        # Apply alpha compositing to the colors for all domains
        # other than those the user wishes to highlight
        for domain_id in self.col_spec:
            if domain_id not in domains:
                r, g, b = self.col_spec[domain_id]
                r = int(((1-alpha) * background[0]) + (alpha * r))
                g = int(((1-alpha) * background[1]) + (alpha * g))
                b = int(((1-alpha) * background[2]) + (alpha * b))
                self._col_spec[domain_id] = (r, g, b)

    def get_plot_xml(self):
        """Return XML representation of the plot

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing plot data

        """

        element = ET.Element("plot")
        element.set("id", str(self._id))
        element.set("filename", self._filename)
        element.set("color", self._color)
        element.set("type", self._type)

        if self._type is 'slice':
            element.set("basis", self._basis)

        subelement = ET.SubElement(element, "origin")
        subelement.text = ' '.join(map(str, self._origin))

        subelement = ET.SubElement(element, "width")
        subelement.text = ' '.join(map(str, self._width))

        subelement = ET.SubElement(element, "pixels")
        subelement.text = ' '.join(map(str, self._pixels))

        if self._mask_background is not None:
            subelement = ET.SubElement(element, "background")
            subelement.text = ' '.join(map(str, self._background))

        if self._col_spec is not None:
            for key in self._col_spec:
                subelement = ET.SubElement(element, "col_spec")
                subelement.set("id", str(key))
                subelement.set("rgb", ' '.join(map(
                    str, self._col_spec[key])))

        if self._mask_components is not None:
            subelement = ET.SubElement(element, "mask")
            subelement.set("components", ' '.join(map(
                str, self._mask_components)))
            subelement.set("background", ' '.join(map(
                str, self._mask_background)))

        return element


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
        super(Plots, self).__init__(Plot, 'plots collection')
        self._plots_file = ET.Element("plots")
        if plots is not None:
            self += plots

    def add_plot(self, plot):
        """Add a plot to the file.

        .. deprecated:: 0.8
            Use :meth:`Plots.append` instead.

        Parameters
        ----------
        plot : openmc.Plot
            Plot to add

        """
        warnings.warn("Plots.add_plot(...) has been deprecated and may be "
                      "removed in a future version. Use Plots.append(...) "
                      "instead.", DeprecationWarning)
        self.append(plot)

    def append(self, plot):
        """Append plot to collection

        Parameters
        ----------
        plot : openmc.Plot
            Plot to append

        """
        super(Plots, self).append(plot)

    def insert(self, index, plot):
        """Insert plot before index

        Parameters
        ----------
        index : int
            Index in list
        plot : openmc.Plot
            Plot to insert

        """
        super(Plots, self).insert(index, plot)

    def remove_plot(self, plot):
        """Remove a plot from the file.

        .. deprecated:: 0.8
            Use :meth:`Plots.remove` instead.

        Parameters
        ----------
        plot : openmc.Plot
            Plot to remove

        """
        warnings.warn("Plots.remove_plot(...) has been deprecated and may be "
                      "removed in a future version. Use Plots.remove(...) "
                      "instead.", DeprecationWarning)
        self.remove(plot)

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

        This routine generates a color scheme and applies alpha compositing
        to make all domains except the highlighted ones partially transparent.

        Parameters
        ----------
        geometry : openmc.Geometry
            The geometry for which the plot is defined
        domains : Iterable of Integral
            A collection of the domain IDs to highlight in the plot
        seed : Integral
            The random number seed used to generate the color scheme
        alpha : Real in [0,1]
            The value to apply in alpha compisiting
        background : 3-tuple of Integral or 'white' or 'black' or 'gray'
            The background color to apply in alpha compisiting

        """

        for plot in self:
            plot.highlight_domains(geometry, domains, seed, alpha, background)

    def _create_plot_subelements(self):
        for plot in self:
            xml_element = plot.get_plot_xml()

            if len(plot._name) > 0:
                self._plots_file.append(ET.Comment(plot._name))

            self._plots_file.append(xml_element)

    def export_to_xml(self):
        """Create a plots.xml file that can be used by OpenMC.

        """

        # Reset xml element tree
        self._plots_file.clear()

        self._create_plot_subelements()

        # Clean the indentation in the file to be user-readable
        clean_xml_indentation(self._plots_file)

        # Write the XML Tree to the plots.xml file
        tree = ET.ElementTree(self._plots_file)
        tree.write("plots.xml", xml_declaration=True,
                             encoding='utf-8', method="xml")
