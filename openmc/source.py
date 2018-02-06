from numbers import Real
import sys
from xml.etree import ElementTree as ET

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

    """

    def __init__(self, space=None, angle=None, energy=None, filename=None, strength=1.0):
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

    def to_xml_element(self):
        """Return XML representation of the source

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing source data

        """
        element = ET.Element("source")
        element.set("strength", str(self.strength))
        if self.file is not None:
            element.set("file", self.file)
        if self.space is not None:
            element.append(self.space.to_xml_element())
        if self.angle is not None:
            element.append(self.angle.to_xml_element())
        if self.energy is not None:
            element.append(self.energy.to_xml_element('energy'))
        return element
