#!/usr/bin/env python

import openmc
from openmc.clean_xml import *
from xml.etree import ElementTree as ET

def reset_auto_ids():
    openmc.reset_auto_material_id()
    openmc.reset_auto_surface_id()
    openmc.reset_auto_cell_id()
    openmc.reset_auto_universe_id()


class Geometry(object):

    def __init__(self):

        # Initialize Geometry class attributes
        self._root_universe = None
        self._offsets = dict()


    def getOffset(self, path, filter_offset):
        """
        Returns the corresponding location in the results array for a given
        path and filter number.

        Parameters
        ----------
        path : list
                A list of IDs that form the path to the target. It should begin
                with 0 for the base universe, and should cover every universe,
                cell, and lattice passed through. For the case of the lattice,
                a tuple should be provided to indicate which coordinates in the
                lattice should be entered. This should be in the
                form: (lat_id, i_x, i_y, i_z)

        filter_offset : int
                An integer that specifies which offset map the filter is using

        """

        # Return memoize'd offset if possible
        if (path, filter_offset) in self._offsets:
            offset = self._offsets[(path, filter_offset)]

        # Begin recursive call to compute offset starting with the base Universe
        else:
            offset = self._root_universe.getOffset(path, filter_offset)
            self._offsets[(path, filter_offset)] = offset

        # Return the final offset
        return offset


    def getAllCells(self):
        return self._root_universe.getAllCells()


    def getAllUniverses(self):
        return self._root_universe.getAllUniverses()


    def getAllNuclides(self):

        nuclides = dict()
        materials = self.getAllMaterials()

        for material in materials:
            nuclides.update(material.getAllNuclides())

        return nuclides


    def getAllMaterials(self):

        material_cells = self.getAllMaterialCells()
        materials = set()

        for cell in material_cells:
            materials.add(cell._fill)

        return list(materials)


    def getAllMaterialCells(self):

        all_cells = self.getAllCells()
        material_cells = set()

        for cell_id, cell in all_cells.items():
            if cell._type == 'normal':
                material_cells.add(cell)

        return list(material_cells)


    def getAllMaterialUniverses(self):

        all_universes = self.getAllUniverses()
        material_universes = set()

        for universe_id, universe in all_universes.items():

            cells = universe._cells

            for cell_id, cell in cells.items():
                if cell._type == 'normal':
                    material_universes.add(universe)

        return list(material_universes)


    def setRootUniverse(self, root_universe):

        if not isinstance(root_universe, openmc.Universe):
            msg = 'Unable to add root Universe {0} to Geometry since ' \
                  'it is not a Universe'.format(root_universe)
            raise ValueError(msg)

        elif root_universe._id != 0:
            msg = 'Unable to add root Universe {0} to Geometry since ' \
                  'it has ID={1} instead of ' \
                  'ID=0'.format(root_universe, root_universe._id)
            raise ValueError(msg)

        self._root_universe = root_universe



class GeometryFile(object):

    def __init__(self):

        # Initialize GeometryFile class attributes
        self._geometry = None
        self._geometry_file = ET.Element("geometry")


    def setGeometry(self, geometry):

        if not isinstance(geometry, Geometry):
            msg = 'Unable to set the Geometry to {0} for the GeometryFile ' \
                        'since it is not a Geometry object'.format(geometry)
            raise ValueError(msg)

        self._geometry = geometry


    def exportToXML(self):

        root_universe = self._geometry._root_universe
        root_universe.createXMLSubElement(self._geometry_file)

        # Clean the indentation in the file to be user-readable
        sort_xml_elements(self._geometry_file)
        clean_xml_indentation(self._geometry_file)

        # Write the XML Tree to the materials.xml file
        tree = ET.ElementTree(self._geometry_file)
        tree.write("geometry.xml", xml_declaration=True,
                             encoding='utf-8', method="xml")
