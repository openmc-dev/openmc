import os
import xml.etree.ElementTree as ET

import h5py

from openmc.mixin import EqualityMixin
from openmc.clean_xml import clean_xml_indentation


class DataLibrary(EqualityMixin):
    def __init__(self):
        self.libraries = []

    def register_file(self, filename, filetype='neutron'):
        h5file = h5py.File(filename, 'r')

        materials = []
        for name in h5file:
            materials.append(name)

        library = {'path': filename, 'type': filetype, 'materials': materials}
        self.libraries.append(library)

    def export_to_xml(self, path='cross_sections.xml'):
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
