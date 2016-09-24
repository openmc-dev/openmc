import os
import xml.etree.ElementTree as ET

import h5py

from openmc.mixin import EqualityMixin
from openmc.clean_xml import clean_xml_indentation


class DataLibrary(EqualityMixin):
    """Collection of cross section data libraries.

    Attributes
    ----------
    libraries : list of dict
        List in which each item is a dictionary summarizing cross section data
        from a single file. The dictionary has keys 'path', 'type', and
        'materials'.

    """

    def __init__(self):
        self.libraries = []

    def register_file(self, filename):
        """Register a file with the data library.

        Parameters
        ----------
        filename : str
            Path to the file to be registered.

        """
        h5file = h5py.File(filename, 'r')

        materials = []
        filetype = 'neutron'
        for name in h5file:
            if name.startswith('c_'):
                filetype = 'thermal'
            materials.append(name)

        library = {'path': filename, 'type': filetype, 'materials': materials}
        self.libraries.append(library)

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
            [lib['path'] for lib in self.libraries]))
        if common_dir == '':
            common_dir = '.'

        directory = os.path.relpath(common_dir, os.path.dirname(path))
        if directory != '.':
            dir_element = ET.SubElement(root, "directory")
            dir_element.text = directory

        for library in self.libraries:
            lib_element = ET.SubElement(root, "library")
            lib_element.set('materials', ' '.join(library['materials']))
            lib_element.set('path', os.path.relpath(library['path'], common_dir))
            lib_element.set('type', library['type'])

        # Clean the indentation to be user-readable
        clean_xml_indentation(root)

        # Write XML file
        tree = ET.ElementTree(root)
        tree.write(path, xml_declaration=True, encoding='utf-8',
                   method='xml')
