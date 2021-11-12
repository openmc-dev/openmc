from abc import ABC
from collections.abc import Iterable
from numbers import Real, Integral
import warnings
from xml.etree import ElementTree as ET

import numpy as np

import openmc.checkvalue as cv
import openmc
from ._xml import get_text
from .mixin import IDManagerMixin
from .surface import _BOUNDARY_TYPES


class WeightWindowMesh(Mesh):
    """A Class to handle the creation of a specific weight window mesh for
    a particular problem
    """

    def __init__(self,Mesh):
        self.mesh = Mesh
        return

    @property
    def mesh(self):
        return self.mesh

class WeightWindow():
    """ A class to handle the creation of a set of specific weight window
    paramaters - a variance reduction class may have several of these
    """

    def __init__(self, particle_types="",
                 upper_bound = 5,
                 survival = 3,
                 max_split = 10,
                 multiplier = 1.0
                 weight_cutoff = 1.e-38,
                 energy_bins = [1.e-50,1e50],
                 lower_ww_bounds = [],
                 upper_ww_bounds = [],
                 mesh = None):

        self.particle_types = particle_types
        self.upper_bound = upper_bound
        self.survival = survival
        self.max_split = max_split
        self.weight_cutoff = weight_cutoff
        self.energy_bins = energy_bins
        self.lower_ww_bounds = lower_ww_bounds
        self.upper_ww_bounds = upper_ww_bounds

        self.ww_mesh = WeightWindowMesh(mesh)

        return

    def

class VarianceReduction():
    """ A class to handle various Variance Reduction operations
    """
    def __init__(self, weight_windows = None):
        self.weight_windows = weight_windows
        return

    def export_to_xml(self, path='variance_reduction.xml'):
        """ Exports the variance reduction parameters to the
        the variance_reduction file

        Parameters
        ----------
        path : str
            Path to file to write. Defaults to 'variance_reduction.xml'.

            .. versionadded:: 0.13

        """

        # Create XML representation
        root_element = ET.Element("variance_reduction")
        self.root_element.create_xml_subelement(root_element, memo=set())

        # add the weight windows to the xml file
        if not self.weight_windows:
             return

        for window in self.weight_windows():
            self.root_element.add_child(window.to_xml_element())

        # Check if path is a directory
        p = Path(path)
        if p.is_dir():
            p /= 'variance_reduction.xml'

        # Write the XML Tree to the geometry.xml file
        xml.reorder_attributes(root_element)  # TODO: Remove when support is Python 3.8+
        tree = ET.ElementTree(root_element)
        tree.write(str(p), xml_declaration=True, encoding='utf-8')

        return



class WeightWindow():
    """A class to handle the creation of weight windows

    Parameters
    ----------
    mesh_id : int
        Unique identifier for the mesh
    name : str
        Name of the mesh

    Attributes
    ----------
    id : int
        Unique identifier for the mesh
    name : str
        Name of the mesh
    type : int
        weight window input file type

	--- mesh cell boundary ---
	origin : float
	    lower-left coordinates for mesh
	xmesh : float
	    the locations of the coarse meshes in the X direction
	xints : float
	    the number of fine meshes within corresponding coars meshes in the X direction
	ymesh : float
	    the locations of the coarse meshes in the Y direction
	yints : float
	    the number of fine meshes within corresponding coars meshes in the Y direction
	zmesh : float
	    the locations of the coarse meshes in the Z direction
	zints : float
	    the number of fine meshes within corresponding coars meshes in the Z direction
	--- mesh cell boundary ---

	--- energy group ---
	n_energy_group : Iterable of float
	    energy group for neutron
	p_energy_group : Iterable of float
	    energy group for photon
	--- energy group ---

	--- lower weight window ---
	lower_ww : Iterable of float
	    lower weight window for mesh for neutron and/or photon
	--- lower weight window ---

	--- neutron WWP ---
	n_upper_ratio : float
        upper weight window = upper_ratio * lower weight window
    n_survival_ratio : float
        survival weight = survival_ratio * lower weight window
    n_max_split : int
        max number of split particles
    n_multiplier : float
        multiplier for weight window lower bounds
	--- neutron WWP ---

	--- photon WWP ---
	p_upper_ratio : float
        upper weight window = upper_ratio * lower weight window
    p_survival_ratio : float
        survival weight = survival_ratio * lower weight window
    p_max_split : int
        max number of split particles
    p_multiplier : float
        multiplier for weight window lower bounds
	--- photon WWP ---

	--- source weight biasing ---
    biasing_energy : Iterable of float
	    energy group for weight biasing
	origin_probability : Iterable of float
	    original probability for each group
	biasing : Iterable of float
	    biasing for each energy group
    --- source weight biasing ---

    """

    def __init__(self, mesh_id=None, name=''):
        super().__init__(mesh_id, name)

		self._type = None
		self._origin = None
		self._xmesh = None
		self._xints = None
		self._ymesh = None
		self._yints = None
		self._zmesh = None
		self._zints = None

		self._n_energy_group = None
		self._p_energy_group = None
        self._lower_ww = None

        self._n_upper_ratio = None
		self._n_survival_ratio = None
		self._n_max_split = None
		self._n_multiplier = None

        self._p_upper_ratio = None
		self._p_survival_ratio = None
		self._p_max_split = None
		self._p_multiplier = None

		self._biasing_energy = None
		self._origin_probability = None
		self._biasing = None

    @property
	def mesh(self):
	    return self._mesh

    @property
    def type(self):
        return self._type

    @origin.setter
    def origin(self, origin):
        cv.check_type('origin', origin, Iterable, Real)
        cv.check_length('origin', origin, 1, 3)
        self._origin = origin

    @xmesh.setter
    def xmesh(self, xmesh):
        cv.check_type('xmesh', xmesh, Iterable, Real)
        self._xmesh = xmesh

    @xints.setter
    def xints(self, xints):
        cv.check_type('xints', xints, Iterable, Integral)
        self._xints = xints

    @ymesh.setter
    def ymesh(self, ymesh):
        cv.check_type('ymesh', ymesh, Iterable, Real)
        self._ymesh = ymesh

    @yints.setter
    def yints(self, yints):
        cv.check_type('yints', yints, Iterable, Integral)
        self._yints = yints

    @zmesh.setter
    def zmesh(self, zmesh):
        cv.check_type('zmesh', zmesh, Iterable, Real)
        self._zmesh = zmesh

    @zints.setter
    def zints(self, zints):
        cv.check_type('zints', zints, Iterable, Integral)
        self._zints = zints

    @property
    def n_energy_group(self):
        return self._n_energy_group

    @property
    def p_energy_group(self):
        return self._p_energy_group

    @property
    def lower_ww(self):
        return self._lower_ww

	@property
    def n_upper_ratio(self):
        return self._n_upper_ratio

	@property
    def n_survival_ratio(self):
        return self._n_survival_ratio

	@property
    def n_max_split(self):
        return self._n_max_split

	@property
    def n_multiplier(self):
        return self._n_multiplier

	@property
    def p_upper_ratio(self):
        return self._p_upper_ratio

	@property
    def p_survival_ratio(self):
        return self._p_survival_ratio

	@property
    def p_max_split(self):
        return self._p_max_split

	@property
    def p_multiplier(self):
        return self._p_multiplier


	@mesh.setter
    def mesh(self,grid,grid,grid):
        cv

    @type.setter
    def type(self, type):
        cv.check_type('type', type, Integral)
        self._type = type

	@n_energy_group.setter
	def n_energy_group(self, n_energy_group):
	    cv.check_type('n_energy_group', n_energy_group, Iterable, Real)
        self._n_energy_group = n_energy_group

	@p_energy_group.setter
	def p_energy_group(self, p_energy_group):
	    cv.check_type('p_energy_group', p_energy_group, Iterable, Real)
        self._p_energy_group = p_energy_group

	@lower_ww.setter
	def lower_ww(self, lower_ww):
	    cv.check_type('lower_ww', lower_ww, Iterable, Real)
        self._lower_ww = lower_ww

	@n_upper_ratio.setter
    def n_upper_ratio(self, n_upper_ratio):
        cv.check_type('n_upper_ratio', n_upper_ratio, Real)
		cv.check_greater_than('n_upper_ratio', n_upper_ratio, 1.0)
        self._n_upper_ratio = n_upper_ratio

	@n_survival_ratio.setter
    def n_survival_ratio(self, n_survival_ratio):
        cv.check_type('n_survival_ratio', n_survival_ratio, Real)
		cv.check_greater_than('n_survival_ratio', n_survival_ratio, 1.0)
        self._n_survival_ratio = n_survival_ratio

	@n_max_split.setter
    def n_max_split(self, n_max_split):
        cv.check_type('n_max_split', n_max_split, Integral)
		cv.check_greater_than('n_max_split', n_max_split, 1.0)
        self._n_max_split = n_max_split

	@n_multiplier.setter
    def n_multiplier(self, n_multiplier):
        cv.check_type('n_multiplier', n_multiplier, Real)
		cv.check_greater_than('n_multiplier', n_multiplier, 0.0)
        self._n_multiplier = n_multiplier

	@p_upper_ratio.setter
    def p_upper_ratio(self, p_upper_ratio):
        cv.check_type('p_upper_ratio', p_upper_ratio, Real)
		cv.check_greater_than('p_upper_ratio', p_upper_ratio, 1.0)
        self._p_upper_ratio = p_upper_ratio

	@p_survival_ratio.setter
    def p_survival_ratio(self, p_survival_ratio):
        cv.check_type('p_survival_ratio', p_survival_ratio, Real)
		cv.check_greater_than('p_survival_ratio', p_survival_ratio, 1.0)
        self._p_survival_ratio = p_survival_ratio

	@p_max_split.setter
    def p_max_split(self, p_max_split):
        cv.check_type('p_max_split', p_max_split, Integral)
		cv.check_greater_than('p_max_split', p_max_split, 1.0)
        self._p_max_split = p_max_split

	@p_multiplier.setter
    def p_multiplier(self, p_multiplier):
        cv.check_type('p_multiplier', p_multiplier, Real)
		cv.check_greater_than('p_multiplier', p_multiplier, 0.0)
        self._p_multiplier = p_multiplier

	@biasing_energy.setter
	def biasing_energy(self, biasing_energy):
	    cv.check_type('biasing_energy', biasing_energy, Iterable, Real)
        self._biasing_energy = biasing_energy

	@origin_probability.setter
	def origin_probability(self, origin_probability):
	    cv.check_type('origin_probability', origin_probability, Iterable, Real)
        self._origin_probability = origin_probability

	@biasing.setter
	def biasing(self, biasing):
	    cv.check_type('biasing', biasing, Iterable, Real)
        self._biasing = biasing

    @classmethod
    def to_xml_element(self):
        """Return XML representation of the WeightWindowMesh

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing WeightWindowMesh data

        """

        element = ET.Element("weightwindow")

        if self._type is not None:
            subelement = ET.SubElement(element, "type")
            subelement.text = str(self._type)

        if self._origin is not None:
            subelement = ET.SubElement(element, "origin")
            subelement.text = ' '.join(map(str, self._origin))

        if self._xmesh is not None:
            subelement = ET.SubElement(element, "xmesh")
            subelement.text = ' '.join(map(str, self._xmesh))

        if self._xints is not None:
            subelement = ET.SubElement(element, "xints")
            subelement.text = ' '.join(map(str, self._xints))

        if self._ymesh is not None:
            subelement = ET.SubElement(element, "ymesh")
            subelement.text = ' '.join(map(str, self._ymesh))

        if self._yints is not None:
            subelement = ET.SubElement(element, "yints")
            subelement.text = ' '.join(map(str, self._yints))

        if self._zmesh is not None:
            subelement = ET.SubElement(element, "zmesh")
            subelement.text = ' '.join(map(str, self._zmesh))

        if self._zints is not None:
            subelement = ET.SubElement(element, "zints")
            subelement.text = ' '.join(map(str, self._zints))

        if self._n_energy_group is not None:
            subelement = ET.SubElement(element, "n_energy_group")
            subelement.text = ' '.join(map(str, self._n_energy_group))

        if self._p_energy_group is not None:
            subelement = ET.SubElement(element, "p_energy_group")
            subelement.text = ' '.join(map(str, self._p_energy_group))

        if self._lower_ww is not None:
            subelement = ET.SubElement(element, "lower_ww")
            subelement.text = ' '.join(map(str, self._lower_ww))

        if self._n_upper_ratio is not None:
            subelement = ET.SubElement(element, "n_upper_ratio")
            subelement.text = self._n_upper_ratio.value

        if self._n_survival_ratio is not None:
            subelement = ET.SubElement(element, "n_survival_ratio")
            subelement.text = self._n_survival_ratio.value

        if self._n_max_split is not None:
            subelement = ET.SubElement(element, "n_max_split")
            subelement.text = self._n_max_split.value

        if self._n_multiplier is not None:
            subelement = ET.SubElement(element, "n_multiplier")
            subelement.text = self._n_multiplier.value

        if self._p_upper_ratio is not None:
            subelement = ET.SubElement(element, "p_upper_ratio")
            subelement.text = self._p_upper_ratio.value

        if self._p_survival_ratio is not None:
            subelement = ET.SubElement(element, "p_survival_ratio")
            subelement.text = self._p_survival_ratio.value

        if self._p_max_split is not None:
            subelement = ET.SubElement(element, "p_max_split")
            subelement.text = self._p_max_split.value

        if self._p_multiplier is not None:
            subelement = ET.SubElement(element, "p_multiplier")
            subelement.text = self._p_multiplier.value

        if self._biasing_energy is not None:
            subelement = ET.SubElement(element, "biasing_energy")
            subelement.text = ' '.join(map(str, self._biasing_energy))

        if self._origin_probability is not None:
            subelement = ET.SubElement(element, "origin_probability")
            subelement.text = ' '.join(map(str, self._origin_probability))

        if self._biasing is not None:
            subelement = ET.SubElement(element, "biasing")
            subelement.text = ' '.join(map(str, self._biasing))


        return element

    def from_xml_element(cls, elem):
        """Generate WeightWindowMesh from an XML element

        Parameters
        ----------
        elem : xml.etree.ElementTree.Element
            XML element

        Returns
        -------
        openmc.WeightWindowMesh
            WeightWindowMesh generated from XML element

        """

        weightwindowmesh = cls()

        type = get_text(elem, 'type')
        if type is not None:
            weightwindowmesh._type = int(type)

        origin = get_text(elem, 'origin')
        if origin is not None:
            weightwindowmesh._origin = [float(x) for x in origin.split()]

        xmesh = get_text(elem, 'xmesh')
        if xmesh is not None:
            weightwindowmesh._xmesh = [float(x) for x in xmesh.split()]

        xints = get_text(elem, 'xints')
        if xints is not None:
            weightwindowmesh._xints = [int(x) for x in xints.split()]

        ymesh = get_text(elem, 'ymesh')
        if ymesh is not None:
            weightwindowmesh._ymesh = [float(x) for x in ymesh.split()]

        yints = get_text(elem, 'yints')
        if yints is not None:
            weightwindowmesh._yints = [int(x) for x in yints.split()]

        zmesh = get_text(elem, 'zmesh')
        if zmesh is not None:
            weightwindowmesh._zmesh = [float(x) for x in zmesh.split()]

        zints = get_text(elem, 'zints')
        if zints is not None:
            weightwindowmesh._zints = [int(x) for x in zints.split()]

        n_energy_group = get_text(elem, 'n_energy_group')
        if n_energy_group is not None:
            weightwindowmesh._n_energy_group = [float(x) for x in n_energy_group.split()]

        p_energy_group = get_text(elem, 'p_energy_group')
        if p_energy_group is not None:
            weightwindowmesh._p_energy_group = [float(x) for x in p_energy_group.split()]

        lower_ww = get_text(elem, 'lower_ww')
        if lower_ww is not None:
            weightwindowmesh._lower_ww = [float(x) for x in lower_ww.split()]

        n_upper_ratio = get_text(elem, 'n_upper_ratio')
        if n_upper_ratio is not None:
            weightwindowmesh._n_upper_ratio = float(n_upper_ratio)

        n_survival_ratio = get_text(elem, 'n_survival_ratio')
        if n_survival_ratio is not None:
            weightwindowmesh._n_survival_ratio = float(n_survival_ratio)

        n_max_split = get_text(elem, 'n_max_split')
        if n_max_split is not None:
            weightwindowmesh._n_max_split = int(n_max_split)

        n_multiplier = get_text(elem, 'n_multiplier')
        if n_multiplier is not None:
            weightwindowmesh._n_multiplier = float(n_multiplier)

        p_upper_ratio = get_text(elem, 'p_upper_ratio')
        if p_upper_ratio is not None:
            weightwindowmesh._p_upper_ratio = float(p_upper_ratio)

        p_survival_ratio = get_text(elem, 'p_survival_ratio')
        if p_survival_ratio is not None:
            weightwindowmesh._p_survival_ratio = float(p_survival_ratio)

        p_max_split = get_text(elem, 'p_max_split')
        if p_max_split is not None:
            weightwindowmesh._p_max_split = int(p_max_split)

        p_multiplier = get_text(elem, 'p_multiplier')
        if p_multiplier is not None:
            weightwindowmesh._p_multiplier = float(p_multiplier)

        biasing_energy = get_text(elem, 'biasing_energy')
        if biasing_energy is not None:
            weightwindowmesh._biasing_energy = [float(x) for x in biasing_energy.split()]

        origin_probability = get_text(elem, 'origin_probability')
        if origin_probability is not None:
            weightwindowmesh._origin_probability = [float(x) for x in origin_probability.split()]

        biasing = get_text(elem, 'biasing')
        if biasing is not None:
            weightwindowmesh._biasing = [float(x) for x in biasing.split()]

        return weightwindowmesh
