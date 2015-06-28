from xml.etree import ElementTree as ET

import numpy as np

from openmc.checkvalue import *
from openmc.clean_xml import *


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
    width : tuple or list or ndarray
        Width of the plot in each basis direction
    pixels : tuple or list or ndarray
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
    mask_components : tuple or list or ndarray
        Unique id numbers of the cells or materials to plot
    mask_background : tuple or list or ndarray
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

        # Check that the ID is an integer and wasn't already used
        elif not is_integer(plot_id):
            msg = 'Unable to set a non-integer Plot ID {0}'.format(plot_id)
            raise ValueError(msg)

        elif plot_id < 0:
            msg = 'Unable to set Plot ID to {0} since it must be a ' \
                  'non-negative integer'.format(plot_id)
            raise ValueError(msg)

        else:
            self._id = plot_id

    @name.setter
    def name(self, name):

        if not is_string(name):
            msg = 'Unable to set name for Plot ID={0} with a non-string ' \
                  'value {1}'.format(self._id, name)
            raise ValueError(msg)

        else:
            self._name = name

    @width.setter
    def width(self, width):
        if not isinstance(width, (tuple, list, np.ndarray)):
            msg = 'Unable to create Plot ID={0} with width {1} which is not ' \
                   'a Python tuple/list or NumPy array'.format(self._id, width)
            raise ValueError(msg)

        elif len(width) != 2 and len(width) != 3:
            msg = 'Unable to create Plot ID={0} with width {1} since only 2D ' \
                  'and 3D plots are supported'.format(self._id, width)
            raise ValueError(msg)

        for dim in width:
            if not is_integer(dim) and not is_float(dim):
                msg = 'Unable to create Plot ID={0} with width {1} since ' \
                      'each element must be a floating point value or ' \
                      'integer'.format(self._id, width)
                raise ValueError(msg)

        self._width = width

    @origin.setter
    def origin(self, origin):
        if not isinstance(origin, (tuple, list, np.ndarray)):
            msg = 'Unable to create Plot ID={0} with origin {1} which is not ' \
                  'a Python tuple/list or NumPy array'.format(self._id, origin)
            raise ValueError(msg)

        elif len(origin) != 3:
            msg = 'Unable to create Plot ID={0} with origin {1} since only ' \
                  'a 3D coordinate must be input'.format(self._id, origin)
            raise ValueError(msg)

        for dim in origin:
            if not is_integer(dim) and not is_float(dim):
                msg = 'Unable to create Plot ID={0} with origin {1} since ' \
                      'each element must be a floating point value or ' \
                      'integer'.format(self._id, origin)
                raise ValueError(msg)

        self._origin = origin

    @pixels.setter
    def pixels(self, pixels):
        if not isinstance(pixels, (tuple, list, np.ndarray)):
            msg = 'Unable to create Plot ID={0} with pixels {1} which is not ' \
                  'a Python tuple/list or NumPy array'.format(self._id, pixels)
            raise ValueError(msg)

        elif len(pixels) != 2 and len(pixels) != 3:
            msg = 'Unable to create Plot ID={0} with pixels {1} since ' \
                  'only 2D and 3D plots are supported'.format(self._id, pixels)
            raise ValueError(msg)

        for dim in pixels:
            if not is_integer(dim):
                msg = 'Unable to create Plot ID={0} with pixel value {1} ' \
                      'which is not an integer'.format(self._id, dim)
                raise ValueError(msg)

            elif dim < 0:
                msg = 'Unable to create Plot ID={0} with pixel value {1} ' \
                      'which is less than 0'.format(self._id, dim)
                raise ValueError(msg)

        self._pixels = pixels

    @filename.setter
    def filename(self, filename):
        if not is_string(filename):
            msg = 'Unable to create Plot ID={0} with filename {1} which is ' \
                  'not a string'.format(self._id, filename)
            raise ValueError(msg)

        self._filename = filename

    @color.setter
    def color(self, color):
        if not is_string(color):
            msg = 'Unable to create Plot ID={0} with color {1} which is not ' \
                  'a string'.format(self._id, color)
            raise ValueError(msg)

        elif color not in ['cell', 'mat']:
            msg = 'Unable to create Plot ID={0} with color {1} which is not ' \
                  'a cell or mat'.format(self._id, color)
            raise ValueError(msg)

        self._color = color

    @type.setter
    def type(self, type):
        if not is_string(type):
            msg = 'Unable to create Plot ID={0} with type {1} which is not ' \
                  'a string'.format(self._id, type)
            raise ValueError(msg)

        elif type not in ['slice', 'voxel']:
            msg = 'Unable to create Plot ID={0} with type {1} which is not ' \
                  'slice or voxel'.format(self._id, type)
            raise ValueError(msg)

        self._type = type

    @basis.setter
    def basis(self, basis):
        if not is_string(basis):
            msg = 'Unable to create Plot ID={0} with basis {1} which is not ' \
                  'a string'.format(self._id, basis)
            raise ValueError(msg)

        elif basis not in ['xy', 'xz', 'yz']:
            msg = 'Unable to create Plot ID={0} with basis {1} which is not ' \
                  'xy, xz, or yz'.format(self._id, basis)
            raise ValueError(msg)

        self._basis = basis

    @background.setter
    def background(self, background):
        if not isinstance(background, (tuple, list, np.ndarray)):
            msg = 'Unable to create Plot ID={0} with background {1} ' \
                  'which is not a Python tuple/list or NumPy ' \
                  'array'.format(self._id, background)
            raise ValueError(msg)

        elif len(background) != 3:
            msg = 'Unable to create Plot ID={0} with background {1} ' \
                  'which is not 3 integer RGB ' \
                  'values'.format(self._id, background)
            raise ValueError(msg)

        for rgb in background:
            if not is_integer(rgb):
                msg = 'Unable to create Plot ID={0} with background RGB ' \
                      'value {1} which is not an integer'.format(self._id, rgb)
                raise ValueError(msg)

            elif rgb < 0 or rgb > 255:
                msg = 'Unable to create Plot ID={0} with background RGB value ' \
                      '{1} which is not between 0 and 255'.format(self._id, rgb)
                raise ValueError(msg)

        self._background = background

    @col_spec.setter
    def col_spec(self, col_spec):
        if not isinstance(col_spec, dict):
            msg = 'Unable to create Plot ID={0} with col_spec parameter {1} ' \
                  'which is not a Python dictionary of IDs to ' \
                  'pixels'.format(self._id, col_spec)
            raise ValueError(msg)

        for key in col_spec:
            if not is_integer(key):
                msg = 'Unable to create Plot ID={0} with col_spec ID {1} ' \
                      'which is not an integer'.format(self._id, key)
                raise ValueError(msg)

            elif key < 0:
                msg = 'Unable to create Plot ID={0} with col_spec ID {1} ' \
                      'which is less than 0'.format(self._id, key)
                raise ValueError(msg)

            elif not isinstance(col_spec[key], (tuple, list, np.ndarray)):
                msg = 'Unable to create Plot ID={0} with col_spec RGB ' \
                      'values {1} which is not a Python tuple/list or NumPy ' \
                      'array'.format(self._id, col_spec[key])
                raise ValueError(msg)

            elif len(col_spec[key]) != 3:
                msg = 'Unable to create Plot ID={0} with col_spec RGB ' \
                      'values of length {1} since 3 values must be ' \
                      'input'.format(self._id, len(col_spec[key]))
                raise ValueError(msg)

        self._col_spec = col_spec

    @mask_componenets.setter
    def mask_components(self, mask_components):
        if not isinstance(mask_components, (list, tuple, np.ndarray)):
            msg = 'Unable to create Plot ID={0} with mask components {1} ' \
                  'which is not a Python tuple/list or NumPy ' \
                  'array'.format(self._id, mask_components)
            raise ValueError(msg)

        for component in mask_components:
            if not is_integer(component):
                msg = 'Unable to create Plot ID={0} with mask component {1} ' \
                      'which is not an integer'.format(self._id, component)
                raise ValueError(msg)

            elif component < 0:
                msg = 'Unable to create Plot ID={0} with mask component {1} ' \
                      'which is less than 0'.format(self._id, component)
                raise ValueError(msg)

        self._mask_components = mask_components

    @mask_background.setter
    def mask_background(self, mask_background):
        if not isinstance(mask_background, (list, tuple, np.ndarray)):
            msg = 'Unable to create Plot ID={0} with mask background {1} ' \
                  'which is not a Python tuple/list or NumPy ' \
                  'array'.format(self._id, mask_background)
            raise ValueError(msg)

        elif len(mask_background) != 3 and len(mask_background) != 0:
            msg = 'Unable to create Plot ID={0} with mask background ' \
                  '{1} since 3 RGB values must be ' \
                  'input'.format(self._id, mask_background)
            raise ValueError(msg)

        for rgb in mask_background:

            if not is_integer(rgb):
                msg = 'Unable to create Plot ID={0} with mask background RGB ' \
                      'value {1} which is not an integer'.format(self._id, rgb)
                raise ValueError(msg)

            elif rgb < 0 or rgb > 255:
                msg = 'Unable to create Plot ID={0} with mask bacground ' \
                      'RGB value {1} which is not between 0 and ' \
                      '255'.format(self._id, rgb)
                raise ValueError(msg)

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
