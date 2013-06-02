#!/usr/bin/env python2

from xml.dom.minidom import parse

class Geometry(object):

    def __init__(self, filename):
        dom = parse(filename)
        rootElement = dom.firstChild

        cells = rootElement.getElementsByTagName('cell')
        surfaces = rootElement.getElementsByTagName('surfaces')
        lattices = rootElement.getElementsByTagName('lattices')

        self.cells = [Cell(elem) for elem in cells]
        self.surfaces = [Surface(elem) for elem in surfaces]
        self.lattices = [Lattice(elem) for elem in lattices]


class Cell(object):

    def __init__(self, domElement):
        self.parse(domElement)

    def parse(self, element):
        for attribute in ['uid', 'universe', 'fill', 'material', 'surfaces']:
            if element.hasAttribute(attribute):
                setattr(self, attribute, element.getAttribute(attribute))

        # Split strings into lists where necessary
        if hasattr(self, 'surfaces'):
            self.surfaces = self.surfaces.split()


class Surface(object):

    def __init__(self, domElement):
        self.parse(domElement)

    def parse(self, element):
        for attribute in ['uid', 'type', 'coeffs', 'boundary']:
            if element.hasAttribute(attribute):
                setattr(self, attribute, element.getAttribute(attribute))

        # Split strings into lists where necessary
        if hasattr(self, 'coeffs'):
            self.surfaces = self.surfaces.split()
