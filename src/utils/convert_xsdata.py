#!/usr/bin/env python2

import os
import sys
from xml.dom.minidom import getDOMImplementation

types = {1: "neutron", 2: "dosimetry", 3: "thermal"}


class Xsdata(object):

    def __init__(self, filename):
        self._table_dict = {}
        self.tables = []

        for line in open(filename, 'r'):
            words = line.split()

            # If this listing is just an alias listing, only assign the alias
            # attribute
            name = words[1]
            alias = words[0]
            table = self.find_table(name)
            if table:
                if name not in table.alias:
                    table.alias.append(alias)
                continue

            table = XsdataTable()
            table.name = name
            table.type = types[int(words[2])]
            table.zaid = int(words[3])
            table.metastable = int(words[4])
            table.awr = float(words[5])
            table.temperature = 8.6173423e-11 * float(words[6])
            table.binary = int(words[7])
            table.path = words[8]

            self.tables.append(table)
            self._table_dict[name] = table

        # Check for common directory
        self.directory = os.path.dirname(self.tables[0].path)
        for table in self.tables:
            if not table.path.startswith(self.directory):
                self.directory = None
                break

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

        # Add a node for each table
        for table in self.tables:
            node = table.to_xml_node(doc)
            root.appendChild(node)

        return doc

    def find_table(self, name):
        if name in self._table_dict:
            return self._table_dict[name]
        else:
            return None


class XsdataTable(object):

    def __init__(self):
        self.alias = []

    def to_xml_node(self, doc):
        node = doc.createElement("ace_table")
        node.setAttribute("name", self.name)
        for attribute in ["alias", "zaid", "type", "metastable",
                          "awr", "temperature", "binary", "path"]:
            if hasattr(self, attribute):
                # Join string for alias attribute
                if attribute == "alias":
                    if not self.alias:
                        continue
                    string = " ".join(self.alias)
                else:
                    string = "{0}".format(getattr(self, attribute))

                # Skip metastable and binary if 0
                if attribute == "metastable" and self.metastable == 0:
                    continue
                if attribute == "binary" and self.binary == 0:
                    continue

                # Create attribute node
                # nodeAttr = doc.createElement(attribute)
                # text = doc.createTextNode(string)
                # nodeAttr.appendChild(text)
                # node.appendChild(nodeAttr)
                node.setAttribute(attribute, string)
        return node


if __name__ == '__main__':
    # Read command line arguments
    if len(sys.argv) < 3:
        sys.exit("Usage: convert_xsdata.py  xsdataFile  xmlFile")
    xsdataFile = sys.argv[1]
    xmlFile = sys.argv[2]

    # Read xsdata and create XML document object
    xsdataObject = Xsdata(xsdataFile)
    doc = xsdataObject.to_xml()

    # Reduce number of lines
    lines = doc.toprettyxml(indent='  ')
    lines = lines.replace('<alias>\n      ', '<alias>')
    lines = lines.replace('\n    </alias>', '</alias>')
    lines = lines.replace('<zaid>\n      ', '<zaid>')
    lines = lines.replace('\n    </zaid>', '</zaid>')
    lines = lines.replace('<type>\n      ', '<type>')
    lines = lines.replace('\n    </type>', '</type>')
    lines = lines.replace('<awr>\n      ', '<awr>')
    lines = lines.replace('\n    </awr>', '</awr>')
    lines = lines.replace('<temperature>\n      ', '<temperature>')
    lines = lines.replace('\n    </temperature>', '</temperature>')
    lines = lines.replace('<path>\n      ', '<path>')
    lines = lines.replace('\n    </path>', '</path>')
    lines = lines.replace('<metastable>\n      ', '<metastable>')
    lines = lines.replace('\n    </metastable>', '</metastable>')
    lines = lines.replace('<binary>\n      ', '<binary>')
    lines = lines.replace('\n    </binary>', '</binary>')

    # Write document in pretty XML to specified file
    f = open(xmlFile, 'w')
    f.write(lines)
    f.close()
