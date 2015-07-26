from collections import Iterable
from numbers import Real, Integral
from xml.etree import ElementTree as ET
import sys

import numpy as np

from openmc.clean_xml import *
from openmc.checkvalue import (check_type, check_value, check_length,
                               check_greater_than, check_less_than)

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
    def mask_componenets(self):
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
            check_type('plot ID', plot_id, Integral)
            check_greater_than('plot ID', plot_id, 0)
            self._id = plot_id

    @name.setter
    def name(self, name):
        check_type('plot name', name, basestring)
        self._name = name

    @width.setter
    def width(self, width):
        check_type('plot width', width, Iterable, Real)
        check_length('plot width', width, 2, 3)
        self._width = width

    @origin.setter
    def origin(self, origin):
        check_type('plot origin', origin, Iterable, Real)
        check_length('plot origin', origin, 3)
        self._origin = origin

    @pixels.setter
    def pixels(self, pixels):
        check_type('plot pixels', pixels, Iterable, Integral)
        check_length('plot pixels', pixels, 2, 3)
        for dim in pixels:
            check_greater_than('plot pixels', dim, 0)
        self._pixels = pixels

    @filename.setter
    def filename(self, filename):
        check_type('filename', filename, basestring)
        self._filename = filename

    @color.setter
    def color(self, color):
        check_type('plot color', color, basestring)
        check_value('plot color', color, ['cell', 'mat'])
        self._color = color

    @type.setter
    def type(self, plottype):
        check_type('plot type', plottype, basestring)
        check_value('plot type', plottype, ['slice', 'voxel'])
        self._type = plottype

    @basis.setter
    def basis(self, basis):
        check_type('plot basis', basis, basestring)
        check_value('plot basis', basis, ['xy', 'xz', 'yz'])
        self._basis = basis

    @background.setter
    def background(self, background):
        check_type('plot background', background, Iterable, Integral)
        check_length('plot background', background, 3)
        for rgb in background:
            check_greater_than('plot background',rgb, 0, True)
            check_less_than('plot background', rgb, 256)
        self._background = background

    @col_spec.setter
    def col_spec(self, col_spec):
        check_type('plot col_spec parameter', col_spec, dict, Integral)

        for key in col_spec:
            if key < 0:
                msg = 'Unable to create Plot ID={0} with col_spec ID {1} ' \
                      'which is less than 0'.format(self._id, key)
                raise ValueError(msg)

            elif not isinstance(col_spec[key], Iterable):
                msg = 'Unable to create Plot ID={0} with col_spec RGB values' \
                      ' {1} which is not iterable'.format(self._id, col_spec[key])
                raise ValueError(msg)

            elif len(col_spec[key]) != 3:
                msg = 'Unable to create Plot ID={0} with col_spec RGB ' \
                      'values of length {1} since 3 values must be ' \
                      'input'.format(self._id, len(col_spec[key]))
                raise ValueError(msg)

        self._col_spec = col_spec

    @mask_componenets.setter
    def mask_components(self, mask_components):
        check_type('plot mask_components', mask_components, Iterable, Integral)
        for component in mask_components:
            check_greater_than('plot mask_components', component, 0, True)
        self._mask_components = mask_components

    @mask_background.setter
    def mask_background(self, mask_background):
        check_type('plot mask background', mask_background, Iterable, Integral)
        check_length('plot mask background', mask_background, 3)
        for rgb in mask_background:
            check_greater_than('plot mask background', rgb, 0, True)
            check_less_than('plot mask background', rgb, 256)
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


class PlotsFile(object):
    """Plots file used for an OpenMC simulation. Corresponds directly to the
    plots.xml input file.

    """

    def __init__(self):
        # Initialize PlotsFile class attributes
        self._plots = []
        self._plots_file = ET.Element("plots")

    def add_plot(self, plot):
        """Add a plot to the file.

        Parameters
        ----------
        plot : Plot
            Plot to add

        """

        if not isinstance(plot, Plot):
            msg = 'Unable to add a non-Plot {0} to the PlotsFile'.format(plot)
            raise ValueError(msg)

        self._plots.append(plot)

    def remove_plot(self, plot):
        """Remove a plot from the file.

        Parameters
        ----------
        plot : Plot
            Plot to remove

        """

        self._plots.remove(plot)

    def _create_plot_subelements(self):
        for plot in self._plots:
            xml_element = plot.get_plot_xml()

            if len(plot._name) > 0:
                self._plots_file.append(ET.Comment(plot._name))

            self._plots_file.append(xml_element)

    def export_to_xml(self):
        """Create a plots.xml file that can be used by OpenMC.

        """

        self._create_plot_subelements()

        # Clean the indentation in the file to be user-readable
        clean_xml_indentation(self._plots_file)

        # Write the XML Tree to the plots.xml file
        tree = ET.ElementTree(self._plots_file)
        tree.write("plots.xml", xml_declaration=True,
                             encoding='utf-8', method="xml")
