from numbers import Real
import sys
from xml.etree import ElementTree as ET

from openmc.stats.univariate import (Univariate, Discrete, Uniform, Maxwell,
                                     Watt, Normal, Muir, Tabular)
from openmc.stats.multivariate import (UnitSphere, Spatial, PolarAzimuthal,
                                       Isotropic, Monodirectional, Box, Point,
                                       CartesianIndependent)
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

    """

    def __init__(self, space=None, angle=None, energy=None, filename=None,
                 strength=1.0, particle='neutron'):
        self._space = None
        self._angle = None
        self._energy = None
        self._file = None

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

        strength = elem.find('strength')
        if strength is not None:
            source.strength = float(strength.text)

        particle = elem.find('particle')
        if particle is not None:
            source.particle = particle.text

        filename = elem.find('file')
        if filename is not None:
            source.file = filename.text

        space = elem.find('space')
        if space is not None:
            space_type = space.get('type')
            if space_type == 'cartesian':
                source.space = CartesianIndependent.from_xml_element(space)
            elif space_type == 'box' or space_type == 'fission':
                source.space = Box.from_xml_element(space)
            elif space_type == 'point':
                source.space = Point.from_xml_element(space)

        angle = elem.find('angle')
        if angle is not None:
            angle_type = angle.get('type')
            if angle_type == 'mu-phi':
                source.angle = PolarAzimuthal.from_xml_element(angle)
            elif angle_type == 'isotropic':
                source.angle = Isotropic.from_xml_element(angle)
            elif angle_type == 'monodirectional':
                source.angle = Monodirectional.from_xml_element(angle)

        energy = elem.find('energy')
        if energy is not None:
            energy_type = energy.get('type')
            if energy_type == 'discrete':
                source.energy = Discrete.from_xml_element(energy)
            elif energy_type == 'uniform':
                source.energy = Uniform.from_xml_element(energy)
            elif energy_type == 'maxwell':
                source.energy = Maxwell.from_xml_element(energy)
            elif energy_type == 'watt':
                source.energy = Watt.from_xml_element(energy)
            elif energy_type == 'normal':
                source.energy = Normal.from_xml_element(energy)
            elif energy_type == 'muir':
                source.energy = Muir.from_xml_element(energy)
            elif energy_type == 'tabular':
                source.energy = Tabular.from_xml_element(energy)

        return source
