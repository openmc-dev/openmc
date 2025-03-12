import os
import pathlib

import h5py
import lxml.etree as ET

import openmc
from openmc._xml import clean_indentation, reorder_attributes


class DataLibrary(list):
    """Collection of cross section data libraries.

    This class behaves like a list where each item is a dictionary summarizing
    cross section data from a single file. The dictionary has keys 'path',
    'type', and 'materials'.

    .. versionchanged:: 0.14.0
        This class now behaves like a list rather than requiring you to access
        the list of libraries through a special attribute.

    """

    def __init__(self):
        super().__init__()

    @property
    def libraries(self):
        # For backwards compatibility
        return self

    def get_by_material(self, name, data_type='neutron'):
        """Return the library dictionary containing a given material.

        Parameters
        ----------
        name : str
            Name of material, e.g. 'Am241'
        data_type : str
            Name of data type, e.g. 'neutron', 'photon', 'wmp', or 'thermal'

            .. versionadded:: 0.12

        Returns
        -------
        library : dict or None
            Dictionary summarizing cross section data from a single file;
            the dictionary has keys 'path', 'type', and 'materials'.

        """
        for library in self:
            if name in library['materials'] and data_type in library['type']:
                return library
        return None

    def remove_by_material(self, name: str, data_type='neutron'):
        """Remove the library dictionary containing a specific material

        Parameters
        ----------
        name : str
            Name of material, e.g. 'Am241'
        data_type : str
            Name of data type, e.g. 'neutron', 'photon', 'wmp', or 'thermal'

        """
        library = self.get_by_material(name, data_type)
        if library is not None:
            self.remove(library)

    def register_file(self, filename):
        """Register a file with the data library.

        Parameters
        ----------
        filename : str or Path
            Path to the file to be registered.
            If an ``xml`` file, treat as the depletion chain file without
            materials.

        """
        if not isinstance(filename, pathlib.Path):
            path = pathlib.Path(filename)
        else:
            path = filename

        if path.suffix == '.xml':
            filetype = 'depletion_chain'
            materials = []
        elif path.suffix == '.h5':
            with h5py.File(path, 'r') as h5file:
                filetype = h5file.attrs['filetype'].decode()[5:]
                materials = list(h5file)
        else:
            raise ValueError(
                f"File type {path.name} not supported by {self.__class__.__name__}")

        library = {'path': str(path), 'type': filetype, 'materials': materials}
        self.append(library)

    def export_to_xml(self, path='cross_sections.xml'):
        """Export cross section data library to an XML file.

        Parameters
        ----------
        path : str
            Path to file to write. Defaults to 'cross_sections.xml'.

        """
        root = ET.Element('cross_sections')

        # Determine common directory for library paths
        common_dir = os.path.dirname(os.path.commonprefix(
            [lib['path'] for lib in self]))
        if common_dir == '':
            common_dir = '.'

        if os.path.relpath(common_dir, os.path.dirname(str(path))) != '.':
            dir_element = ET.SubElement(root, "directory")
            dir_element.text = os.path.realpath(common_dir)

        for library in self:
            if library['type'] == "depletion_chain":
                lib_element = ET.SubElement(root, "depletion_chain")
            else:
                lib_element = ET.SubElement(root, "library")
                lib_element.set('materials', ' '.join(library['materials']))
            lib_element.set('path', os.path.relpath(library['path'], common_dir))
            lib_element.set('type', library['type'])

        # Clean the indentation to be user-readable
        clean_indentation(root)

        # Write XML file
        reorder_attributes(root)  # TODO: Remove when support is Python 3.8+
        tree = ET.ElementTree(root)
        tree.write(str(path), xml_declaration=True, encoding='utf-8',
                   method='xml')

    @classmethod
    def from_xml(cls, path=None):
        """Read cross section data library from an XML file.

        Parameters
        ----------
        path : str, optional
            Path to XML file to read. If not provided,
            openmc.config['cross_sections'] will be used.

        Returns
        -------
        data : openmc.data.DataLibrary
            Data library object initialized from the provided XML

        """

        data = cls()

        # If path is None, get the cross sections from the global configuration
        if path is None:
            path = openmc.config.get('cross_sections')

        # Check to make sure we picked up cross sections
        if path is None:
            raise ValueError("Either path or openmc.config['cross_sections'] "
                             "must be set")

        tree = ET.parse(path)
        root = tree.getroot()
        if root.find('directory') is not None:
            directory = root.find('directory').text
        else:
            directory = os.path.dirname(path)

        for lib_element in root.findall('library'):
            filename = os.path.join(directory, lib_element.attrib['path'])
            filetype = lib_element.attrib['type']
            materials = lib_element.attrib['materials'].split()
            library = {'path': filename, 'type': filetype,
                       'materials': materials}
            data.libraries.append(library)

        # get depletion chain data
        dep_node = root.find("depletion_chain")
        if dep_node is not None:
            filename = os.path.join(directory, dep_node.attrib['path'])
            library = {'path': filename, 'type': 'depletion_chain',
                       'materials': []}
            data.libraries.append(library)

        return data
