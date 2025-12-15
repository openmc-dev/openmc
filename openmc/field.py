from abc import ABC, abstractmethod

import openmc

from .mixin import IDManagerMixin


class ScalarField(IDManagerMixin, ABC):
    """
    """
    def to_xml_element(self):
        """Return XML representation of the field

        Returns
        -------
        element : lxml.etree._Element
            XML element containing mesh data

        """
        elem = ET.Element("field")

        elem.set("id", str(self._id))
        if self.name:
            elem.set("name", self.name)

        return elem


class TemperatureField(ScalarField):
    """
    """
    def __init__(self, mesh, values):
        """
        """
        self.mesh = mesh
        self.values = values

    @classmethod
    def from_mesh_and_values(cls, mesh, values):
        """
        """
        return cls(mesh, values)

    def to_xml_element(self):
        """Return XML representation of the mesh

        Returns
        -------
        element : lxml.etree._Element
            XML element containing mesh data

        """
        element = super().to_xml_element()

        subelement = ET.SubElement(element, "mesh")
        subelement.text = self.mesh.to_xml_element()

        return element
