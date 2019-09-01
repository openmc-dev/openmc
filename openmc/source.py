from collections.abc import Iterable
from numbers import Real
import sys
from xml.etree import ElementTree as ET

from openmc._xml import get_text
from openmc.stats.univariate import Univariate
from openmc.stats.multivariate import UnitSphere, Spatial
import openmc.checkvalue as cv


class Source(object):
    """Distribution of phase space coordinates for source sites.

    Parameters
    ----------
    space : openmc.stats.Spatial, optional
        Spatial distribution of source sites
    angle : openmc.stats.UnitSphere, optional
        Angular distribution of source sites
    energy : openmc.stats.Univariate, optional
        Energy distribution of source sites
    filename : str, optional
        Source file from which sites should be sampled
    strength : Real
        Strength of the source
    particle : {'neutron', 'photon'}
        Source particle type
    pyne_source_mode : int
        Determine the mode of the PyNE source sampler. Required if a PyNE
        source file is used.
    pyne_source_e_bounds : iterable of float
        The energy boundaries of the PyNE source in [eV].
 
    Attributes
    ----------
    space : openmc.stats.Spatial or None
        Spatial distribution of source sites
    angle : openmc.stats.UnitSphere or None
        Angular distribution of source sites
    energy : openmc.stats.Univariate or None
        Energy distribution of source sites
    file : str or None
        Source file from which sites should be sampled
    strength : Real
        Strength of the source
    particle : {'neutron', 'photon'}
        Source particle type
    pyne_source_mode : int
        Determine the mode of the PyNE source sampler. Required if a PyNE
        source file is used.
    pyne_source_e_bounds : iterable of float
        The energy boundaries of the PyNE source in [eV].
 
    """

    def __init__(self, space=None, angle=None, energy=None, filename=None,
                 strength=1.0, particle='neutron'):
        self._space = None
        self._angle = None
        self._energy = None
        self._file = None
        self._pyne_source_mode = None
        self._pyne_source_e_bounds = None

        if space is not None:
            self.space = space
        if angle is not None:
            self.angle = angle
        if energy is not None:
            self.energy = energy
        if filename is not None:
            self.file = filename
        self.strength = strength
        self.particle = particle

    @property
    def file(self):
        return self._file

    @property
    def space(self):
        return self._space

    @property
    def angle(self):
        return self._angle

    @property
    def energy(self):
        return self._energy

    @property
    def strength(self):
        return self._strength

    @property
    def particle(self):
        return self._particle

    @property
    def pyne_source_mode(self):
        return self._pyne_source_mode

    @property
    def pyne_source_e_bounds(self):
        return self._pyne_source_e_bounds
 
    @file.setter
    def file(self, filename):
        cv.check_type('source file', filename, str)
        self._file = filename

    @space.setter
    def space(self, space):
        cv.check_type('spatial distribution', space, Spatial)
        self._space = space

    @angle.setter
    def angle(self, angle):
        cv.check_type('angular distribution', angle, UnitSphere)
        self._angle = angle

    @energy.setter
    def energy(self, energy):
        cv.check_type('energy distribution', energy, Univariate)
        self._energy = energy

    @strength.setter
    def strength(self, strength):
        cv.check_type('source strength', strength, Real)
        cv.check_greater_than('source strength', strength, 0.0, True)
        self._strength = strength

    @particle.setter
    def particle(self, particle):
        cv.check_value('source particle', particle, ['neutron', 'photon'])
        self._particle = particle

    @pyne_source_mode.setter
    def pyne_source_mode(self, pyne_source_mode):
        cv.check_type('pyne source mode', pyne_source_mode, int)
        self._pyne_source_mode = pyne_source_mode

    @pyne_source_e_bounds.setter
    def pyne_source_e_bounds(self, pyne_source_e_bounds):
        cv.check_type('pyne source e_bounds', pyne_source_e_bounds, Iterable,
                Real)
        cv.check_length('pyne source e_bounds minimum length',
                pyne_source_e_bounds, 2, float('inf'))
        cv.check_value('pyne source lowest energy boundary',
                pyne_source_e_bounds[0], [0.0])
        # e_bounds should increse
        for i in range(1, len(pyne_source_e_bounds)):
            cv.check_greater_than('pyne source energy',
                    pyne_source_e_bounds[i], pyne_source_e_bounds[i-1])
        self._pyne_source_e_bounds = pyne_source_e_bounds


    def _create_pyne_source_mode_subelement(self, element):
        if self._pyne_source_mode is not None:
            sub_elem = ET.SubElement(element, 'pyne_source_mode')
            sub_elem.text = str(self._pyne_source_mode).lower()

    def _create_pyne_source_e_bounds_subelement(self, element):
        if self._pyne_source_e_bounds is not None:
            sub_elem = ET.SubElement(element, 'pyne_source_e_bounds')
            sub_elem.text = ' '.join(map(str, self._pyne_source_e_bounds))

    def to_xml_element(self):
        """Return XML representation of the source

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing source data

        """
        element = ET.Element("source")
        element.set("strength", str(self.strength))
        if self.particle != 'neutron':
            element.set("particle", self.particle)
        if self.file is not None:
            element.set("file", self.file)
        if self.space is not None:
            element.append(self.space.to_xml_element())
        if self.angle is not None:
            element.append(self.angle.to_xml_element())
        if self.energy is not None:
            element.append(self.energy.to_xml_element('energy'))
        if self.pyne_source_mode is not None:
            self._create_pyne_source_mode_subelement(element)
        if self.pyne_source_e_bounds is not None:
            self._create_pyne_source_e_bounds_subelement(element)
        return element

    @classmethod
    def from_xml_element(cls, elem):
        """Generate source from an XML element

        Parameters
        ----------
        elem : xml.etree.ElementTree.Element
            XML element

        Returns
        -------
        openmc.Source
            Source generated from XML element

        """
        source = cls()

        strength = get_text(elem, 'strength')
        if strength is not None:
            source.strength = float(strength)

        particle = get_text(elem, 'particle')
        if particle is not None:
            source.particle = particle

        filename = get_text(elem, 'file')
        if filename is not None:
            source.file = filename

        space = elem.find('space')
        if space is not None:
            source.space = Spatial.from_xml_element(space)

        angle = elem.find('angle')
        if angle is not None:
            source.angle = UnitSphere.from_xml_element(angle)

        energy = elem.find('energy')
        if energy is not None:
            source.energy = Univariate.from_xml_element(energy)

        pyne_source_mode = get_text(elem, 'pyne_source_mode')
        if pyne_source_mode is not None:
            self.pyne_source_mode = int(pyne_source_mode)

        pyne_source_e_bounds = get_text(elem, 'pyne_source_e_bounds')
        if pyne_source_e_bounds is not None:
            self.pyne_source_e_bounds = [float(x) for x in 
                    pyne_source_e_bounds.split()]

        return source
