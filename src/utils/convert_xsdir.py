#!/usr/bin/env python2

import os
import sys
from xml.dom.minidom import getDOMImplementation

elements = [None, "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na",
            "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V",
            "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se",
            "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",
            "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba",
            "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho",
            "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt",
            "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",
            "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
            "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg",
            "Cn"]


class Xsdir(object):

    def __init__(self, filename):
        self.f = open(filename, 'r')
        self.filename = os.path.abspath(filename)
        self.directory = os.path.dirname(filename)
        self.awr = {}
        self.tables = []

        self.filetype = set()
        self.recordlength = set()
        self.entries = set()

        # Read first section (DATAPATH)
        line = self.f.readline()
        words = line.split()
        if words:
            if words[0].lower().startswith('datapath'):
                if '=' in words[0]:
                    index = line.index('=')
                    self.datapath = line[index+1:].strip()
                else:
                    if len(line.strip()) > 8:
                        self.datapath = line[8:].strip()
            else:
                self.f.seek(0)

        # Read second section
        line = self.f.readline()
        words = line.split()
        assert len(words) == 3
        assert words[0].lower() == 'atomic'
        assert words[1].lower() == 'weight'
        assert words[2].lower() == 'ratios'

        while True:
            line = self.f.readline()
            words = line.split()

            # Check for end of second section
            if len(words) % 2 != 0 or words[0] == 'directory':
                break

            for zaid, awr in zip(words[::2], words[1::2]):
                self.awr[zaid] = awr

        # Read third section
        while words[0] != 'directory':
            words = self.f.readline().split()

        while True:
            words = self.f.readline().split()
            if not words:
                break

            # Handle continuation lines
            while words[-1] == '+':
                extraWords = self.f.readline().split()
                words = words[:-1] + extraWords
            assert len(words) >= 7

            # Create XsdirTable object and add to line
            table = XsdirTable(self.directory)
            self.tables.append(table)

            # All tables have at least 7 attributes
            table.name = words[0]
            table.awr = float(words[1])
            table.filename = words[2]
            table.access = words[3]
            table.filetype = int(words[4])
            table.location = int(words[5])
            table.length = int(words[6])

            self.filetype.add(table.filetype)

            if len(words) > 7:
                table.recordlength = int(words[7])
                self.recordlength.add(table.recordlength)
            if len(words) > 8:
                table.entries = int(words[8])
                self.entries.add(table.entries)
            if len(words) > 9:
                table.temperature = float(words[9])
            if len(words) > 10:
                table.ptable = (words[10] == 'ptable')

        if len(self.filetype) == 1:
            if 1 in self.filetype:
                self.filetype = 'ascii'
            elif 2 in self.filetype:
                self.filetype = 'binary'
        else:
            self.filetype = None

        if len(self.recordlength) == 1:
            self.recordlength = list(self.recordlength)[0]
        else:
            self.recordlength = None
        if len(self.entries) == 1:
            self.entries = list(self.entries)[0]
        else:
            self.recordlength = None

    def to_xml(self):
        # Create XML document
        impl = getDOMImplementation()
        doc = impl.createDocument(None, "cross_sections", None)

        # Get root element
        root = doc.documentElement

        # Add a directory node
        if self.directory:
            directoryNode = doc.createElement("directory")
            text = doc.createTextNode(self.directory)
            directoryNode.appendChild(text)
            root.appendChild(directoryNode)

            for table in self.tables:
                table.path = os.path.basename(table.path)

        # Add filetype, record_length and entries nodes
        if self.filetype:
            node = doc.createElement("filetype")
            text = doc.createTextNode(self.filetype)
            node.appendChild(text)
            root.appendChild(node)
        if self.recordlength:
            node = doc.createElement("record_length")
            text = doc.createTextNode(str(self.recordlength))
            node.appendChild(text)
            root.appendChild(node)
        if self.entries:
            node = doc.createElement("entries")
            text = doc.createTextNode(str(self.entries))
            node.appendChild(text)
            root.appendChild(node)

        # Add a node for each table
        for table in self.tables:
            if table.name[-1] in ['e', 'p', 'u', 'h', 'g', 'm', 'd']:
                continue
            node = table.to_xml_node(doc)
            root.appendChild(node)

        return doc


class XsdirTable(object):

    def __init__(self, directory=None):
        self.directory = None
        self.name = None
        self.awr = None
        self.filename = None
        self.access = None
        self.filetype = None
        self.location = None
        self.length = None
        self.recordlength = None
        self.entries = None
        self.temperature = None
        self.ptable = False

    @property
    def path(self):
        if self.directory:
            return os.path.join(self.directory, self.filename)
        else:
            return self.filename

    @path.setter
    def path(self, value):
        self.diretory = ''
        self.filename = value

    @property
    def metastable(self):
        # Only valid for neutron cross-sections
        if not self.name.endswith('c'):
            return

        # Handle special case of Am-242 and Am-242m
        if self.zaid == '95242':
            return 1
        elif self.zaid == '95642':
            return 0

        # All other cases
        A = int(self.zaid) % 1000
        if A > 300:
            return 1
        else:
            return 0

    @property
    def alias(self):
        zaid = self.zaid
        if zaid:
            Z = int(zaid[:-3])
            A = zaid[-3:]

            if A == '000':
                s = 'Nat'
            elif zaid == '95242':
                s = '242m'
            elif zaid == '95642':
                s = '242'
            elif int(A) > 300:
                s = str(int(A) - 400) + "m"
            else:
                s = str(int(A))

            return "{0}-{1}.{2}".format(elements[Z], s, self.xs)
        else:
            return None

    @property
    def zaid(self):
        if self.name.endswith('c'):
            return self.name[:self.name.find('.')]
        else:
            return 0

    @property
    def xs(self):
        return self.name[self.name.find('.')+1:]

    def to_xml_node(self, doc):
        node = doc.createElement("ace_table")
        node.setAttribute("name", self.name)
        for attribute in ["alias", "zaid", "type", "metastable", "awr",
                          "temperature", "path", "location"]:
            if hasattr(self, attribute):
                string = str(getattr(self, attribute))

                # Skip metastable and binary if 0
                if attribute == "metastable" and self.metastable == 0:
                    continue

                # Skip any attribute that is none
                if getattr(self, attribute) is None:
                    continue

                # Create attribute node
                node.setAttribute(attribute, string)

        return node


if __name__ == '__main__':
    # Read command line arguments
    if len(sys.argv) < 3:
        sys.exit("Usage: convert_xsdir.py  xsdirFile  xmlFile")
    xsdirFile = sys.argv[1]
    xmlFile = sys.argv[2]

    # Read xsdata and create XML document object
    xsdirObject = Xsdir(xsdirFile)
    doc = xsdirObject.to_xml()

    # Reduce number of lines
    lines = doc.toprettyxml(indent='  ')

    # Write document in pretty XML to specified file
    f = open(xmlFile, 'w')
    f.write(lines)
    f.close()
