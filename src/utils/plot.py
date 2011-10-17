#!/usr/bin/env python

import random
import struct
import sys

import matplotlib.pyplot as plt


class Plot(object):

    def __init__(self):
        self.segments = []
        self.cells = set()

    def plot(self):
        # Create axes
        ax = plt.axes()

        # Create polygon object for each particle track segment
        for s in self.segments:
            if s.cell == 0:
                continue
            ax.add_patch(s.polygon(self.pixel))

        # Set range of x- and y-axes
        plt.axis([self.origin[0] - self.width[0]/2, self.origin[0] + self.width[0]/2,
                  self.origin[1] - self.width[1]/2, self.origin[1] + self.width[1]/2])

        # Display plot
        plt.show()

    def load_file(self, filename):
        # Create binary reader for plot.out file
        plotFile = BinaryReader(filename)

        # Read plot origin, width, basis, and pixel width
        self.origin = plotFile.get_double(3)
        self.width = plotFile.get_double(2)
        self.basis = plotFile.get_double(6)
        self.pixel = plotFile.get_double()

        # Determine bounding x coordinate
        lastXCoord = self.origin[0] + self.width[0]/2

        # Read initial starting coordinate (top-left)
        startingCoord = plotFile.get_double(3)

        while True:
            # Read next coordinate and cell
            nextCoord = plotFile.get_double(3)
            cell = plotFile.get_int()

            # Add segment
            thisSegment = Segment(startingCoord, nextCoord, cell)
            self.segments.append(thisSegment)

            # We also need a set of all the cells to convert to colors later
            self.cells.add(cell)

            # Set the starting coordinate for the next segment as the ending
            # coordinate of the last segment
            startingCoord = nextCoord

            # If the ending coordinate is the bounding value, move to next
            # horizontal ray
            if nextCoord[0] == lastXCoord:
                try:
                    startingCoord = plotFile.get_double(3)
                except BinaryReaderError:
                    break

        # Assign random colors to each cell
        colorDict = {}
        for c in self.cells:
            # Cell 0 indicates a void, so assign white to that 
            if c == 0:
                color = (1,1,1,0)
            else:
                color = (random.random(),
                         random.random(),
                         random.random())
            colorDict[c] = color

        # Set the color of each segment
        for s in self.segments:
            s.color = colorDict[s.cell]


class Segment(object):
    def __init__(self, start, end, cell):
        self.start = start
        self.end = end
        self.cell = cell

    def polygon(self, pixel):
        return plt.Polygon([[self.start[0], self.start[1] + pixel/2],
                            [self.start[0], self.start[1] - pixel/2],
                            [self.end[0], self.end[1] - pixel/2],
                            [self.end[0], self.end[1] + pixel/2]],
                           closed=True, color=self.color)


class BinaryReader(object):

    def __init__(self, filename):
        """
        Initialize instance of Record object.
        """

        fh = open(filename, 'rb')
        self.data = fh.read()
        self.numBytes = len(self.data)

        self.reset()
        self.intSize    = struct.calcsize('i')
        self.longSize   = struct.calcsize('q')
        self.floatSize  = struct.calcsize('f')
        self.doubleSize = struct.calcsize('d')

    def get_data(self, n, typeCode, itemSize):
        """
        Returns one or more items of a specified type at the current
        position within the data list. If more than one item is read,
        the items are returned in a list.
        """
        if self.pos >= self.numBytes:
            raise BinaryReaderError("Already read all data from record")
        
        values = struct.unpack('{0}{1}'.format(n,typeCode), self.data[self.pos:self.pos+itemSize*n])
        self.pos += itemSize * n
        if n == 1:
            return values[0]
        else:
            return list(values)
        

    def get_int(self, n=1):
        """
        Returns one or more 4-byte integers.
        """
        return self.get_data(n,'i',self.intSize)
                             
    def get_long(self, n=1):
        """
        Returns one or more 8-byte integers.
        """
        return self.get_data(n,'q',self.longSize)
                             
    def get_float(self, n=1):
        """
        Returns one or more floats.
        """
        return self.get_data(n,'f',self.floatSize)
                             
    def get_double(self, n=1):
        """
        Returns one or more double
        """
        return self.get_data(n,'d',self.doubleSize)
                             
    def get_string(self, length, n=1):
        """
        Returns a string of a specified length starting at the current
        position in the data list.
        """

        if self.pos >= self.numBytes:
            raise BinaryReaderError("Already read all data from record")
        
        relevantData = self.data[self.pos:self.pos+length*n]
        (s,) = struct.unpack('{0}s'.format(length*n), relevantData)
        self.pos += length*n
        if n == 1:
            return s
        else:
            return [s[i*length:(i+1)*length] for i in range(n)]

    def reset(self):
        self.pos = 0


class BinaryReaderError(Exception):
    """Case class for all binary reader errors."""

    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return repr(self.msg)


if __name__ == '__main__':
    if len(sys.argv) >= 1:
        filename = sys.argv[1]

        p = Plot()
        p.load_file(filename)
        p.plot()
