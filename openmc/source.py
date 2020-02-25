from numbers import Real
import sys
from xml.etree import ElementTree as ET

from openmc._xml import get_text
from openmc.stats.univariate import Univariate
from openmc.stats.multivariate import UnitSphere, Spatial
import openmc.checkvalue as cv


class Source:
    """Distribution of phase space coordinates for source sites.

    Parameters
    ----------
    space : openmc.stats.Spatial
        Spatial distribution of source sites
    angle : openmc.stats.UnitSphere
        Angular distribution of source sites
    energy : openmc.stats.Univariate
        Energy distribution of source sites
    filename : str
        Source file from which sites should be sampled
    library : str
        Path to a custom source library
    strength : float
        Strength of the source
    particle : {'neutron', 'photon'}
        Source particle type

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
    library : str or None
        Path to a custom source library
    strength : float
        Strength of the source
    particle : {'neutron', 'photon'}
        Source particle type

    """

    def __init__(self, space=None, angle=None, energy=None, filename=None,
                 library=None, strength=1.0, particle='neutron'):
        self._space = None
        self._angle = None
        self._energy = None
        self._file = None
        self._library = None

        if space is not None:
            self.space = space
        if angle is not None:
            self.angle = angle
        if energy is not None:
            self.energy = energy
        if filename is not None:
            self.file = filename
        if library is not None:
            self.library = library
        self.strength = strength
        self.particle = particle

    @property
    def file(self):
        return self._file

    @property
    def library(self):
        return self._library

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

    @file.setter
    def file(self, filename):
        cv.check_type('source file', filename, str)
        self._file = filename

    @library.setter
    def library(self, library_name):
        cv.check_type('library', library_name, str)
        self._library = library_name

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
        if self.library is not None:
            element.set("library", self.library)
        if self.space is not None:
            element.append(self.space.to_xml_element())
        if self.angle is not None:
            element.append(self.angle.to_xml_element())
        if self.energy is not None:
            element.append(self.energy.to_xml_element('energy'))
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

        library = get_text(elem, 'library')
        if library is not None:
            source.library = library

        space = elem.find('space')
        if space is not None:
            source.space = Spatial.from_xml_element(space)

        angle = elem.find('angle')
        if angle is not None:
            source.angle = UnitSphere.from_xml_element(angle)

        energy = elem.find('energy')
        if energy is not None:
            source.energy = Univariate.from_xml_element(energy)

        return source
