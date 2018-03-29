from numbers import Integral, Real
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd

import openmc.checkvalue as cv
from . import Filter


class LegendreFilter(Filter):
    r"""Score Legendre expansion moments up to specified order.

    This filter allows scores to be multiplied by Legendre polynomials of the
    change in particle angle ($\mu$) up to a user-specified order.

    Parameters
    ----------
    order : int
        Maximum Legendre polynomial order
    filter_id : int or None
        Unique identifier for the filter

    Attributes
    ----------
    order : int
        Maximum Legendre polynomial order
    id : int
        Unique identifier for the filter
    num_bins : int
        The number of filter bins

    """

    def __init__(self, order, filter_id=None):
        self.order = order
        self.bins = ['P{}'.format(i) for i in range(order + 1)]
        self.id = filter_id

    def __hash__(self):
        string = type(self).__name__ + '\n'
        string += '{: <16}=\t{}\n'.format('\tOrder', self.order)
        return hash(string)

    def __repr__(self):
        string = type(self).__name__ + '\n'
        string += '{: <16}=\t{}\n'.format('\tOrder', self.order)
        string += '{: <16}=\t{}\n'.format('\tID', self.id)
        return string

    @property
    def order(self):
        return self._order

    @order.setter
    def order(self, order):
        cv.check_type('Legendre order', order, Integral)
        cv.check_greater_than('Legendre order', order, 0, equality=True)
        self._order = order

    @classmethod
    def from_hdf5(cls, group, **kwargs):
        if group['type'].value.decode() != cls.short_name.lower():
            raise ValueError("Expected HDF5 data for filter type '"
                             + cls.short_name.lower() + "' but got '"
                             + group['type'].value.decode() + " instead")

        filter_id = int(group.name.split('/')[-1].lstrip('filter '))

        out = cls(group['order'].value, filter_id)

        return out

    def to_xml_element(self):
        """Return XML Element representing the filter.

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing Legendre filter data

        """
        element = ET.Element('filter')
        element.set('id', str(self.id))
        element.set('type', self.short_name.lower())

        subelement = ET.SubElement(element, 'order')
        subelement.text = str(self.order)

        return element


class SpatialLegendreFilter(Filter):
    r"""Score Legendre expansion moments in space up to specified order.

    This filter allows scores to be multiplied by Legendre polynomials of the
    the particle's position along a particular axis, normalized to a given
    range, up to a user-specified order.

    Parameters
    ----------
    order : int
        Maximum Legendre polynomial order
    axis : {'x', 'y', 'z'}
        Axis along which to take the expansion
    minimum : float
        Minimum value along selected axis
    maximum : float
        Maximum value along selected axis
    filter_id : int or None
        Unique identifier for the filter

    Attributes
    ----------
    order : int
        Maximum Legendre polynomial order
    axis : {'x', 'y', 'z'}
        Axis along which to take the expansion
    minimum : float
        Minimum value along selected axis
    maximum : float
        Maximum value along selected axis
    id : int
        Unique identifier for the filter
    num_bins : int
        The number of filter bins

    """

    def __init__(self, order, axis, minimum, maximum, filter_id=None):
        self.order = order
        self.axis = axis
        self.minimum = minimum
        self.maximum = maximum
        self.bins = ['P{}'.format(i) for i in range(order + 1)]
        self.id = filter_id

    def __hash__(self):
        string = type(self).__name__ + '\n'
        string += '{: <16}=\t{}\n'.format('\tOrder', self.order)
        string += '{: <16}=\t{}\n'.format('\tAxis', self.axis)
        string += '{: <16}=\t{}\n'.format('\tMin', self.minimum)
        string += '{: <16}=\t{}\n'.format('\tMax', self.maximum)
        return hash(string)

    def __repr__(self):
        string = type(self).__name__ + '\n'
        string += '{: <16}=\t{}\n'.format('\tOrder', self.order)
        string += '{: <16}=\t{}\n'.format('\tAxis', self.axis)
        string += '{: <16}=\t{}\n'.format('\tMin', self.minimum)
        string += '{: <16}=\t{}\n'.format('\tMax', self.maximum)
        string += '{: <16}=\t{}\n'.format('\tID', self.id)
        return string

    @property
    def order(self):
        return self._order

    @order.setter
    def order(self, order):
        cv.check_type('Legendre order', order, Integral)
        cv.check_greater_than('Legendre order', order, 0, equality=True)
        self._order = order

    @property
    def axis(self):
        return self._axis

    @axis.setter
    def axis(self, axis):
        cv.check_value('axis', axis, ('x', 'y', 'z'))
        self._axis = axis

    @property
    def minimum(self):
        return self._minimum

    @minimum.setter
    def minimum(self, minimum):
        cv.check_type('minimum', minimum, Real)
        self._minimum = minimum

    @property
    def maximum(self):
        return self._maximum

    @maximum.setter
    def maximum(self, maximum):
        cv.check_type('maximum', maximum, Real)
        self._maximum = maximum

    @classmethod
    def from_hdf5(cls, group, **kwargs):
        if group['type'].value.decode() != cls.short_name.lower():
            raise ValueError("Expected HDF5 data for filter type '"
                             + cls.short_name.lower() + "' but got '"
                             + group['type'].value.decode() + " instead")

        filter_id = int(group.name.split('/')[-1].lstrip('filter '))
        order = group['order'].value
        axis = group['axis'].value.decode()
        min_, max_ = group['min'].value, group['max'].value

        return cls(order, axis, min_, max_, filter_id)

    def to_xml_element(self):
        """Return XML Element representing the filter.

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing Legendre filter data

        """
        element = ET.Element('filter')
        element.set('id', str(self.id))
        element.set('type', self.short_name.lower())

        subelement = ET.SubElement(element, 'order')
        subelement.text = str(self.order)
        subelement = ET.SubElement(element, 'axis')
        subelement.text = self.axis
        subelement = ET.SubElement(element, 'min')
        subelement.text = str(self.minimum)
        subelement = ET.SubElement(element, 'max')
        subelement.text = str(self.maximum)

        return element


class SphericalHarmonicsFilter(Filter):
    r"""Score spherical harmonic expansion moments up to specified order.

    Parameters
    ----------
    order : int
        Maximum spherical harmonics order
    filter_id : int or None
        Unique identifier for the filter

    Attributes
    ----------
    order : int
        Maximum spherical harmonics order
    id : int
        Unique identifier for the filter
    cosine : {'scatter', 'particle'}
        How to handle the cosine term.
    num_bins : int
        The number of filter bins

    """

    def __init__(self, order, filter_id=None):
        self.order = order
        self.id = filter_id
        self.bins = ['Y{},{}'.format(n, m)
                     for n in range(order + 1)
                     for m in range(-n, n + 1)]
        self._cosine = 'particle'

    def __hash__(self):
        string = type(self).__name__ + '\n'
        string += '{: <16}=\t{}\n'.format('\tOrder', self.order)
        string += '{: <16}=\t{}\n'.format('\tCosine', self.cosine)
        return hash(string)

    def __repr__(self):
        string = type(self).__name__ + '\n'
        string += '{: <16}=\t{}\n'.format('\tOrder', self.order)
        string += '{: <16}=\t{}\n'.format('\tCosine', self.cosine)
        string += '{: <16}=\t{}\n'.format('\tID', self.id)
        return string

    @property
    def order(self):
        return self._order

    @order.setter
    def order(self, order):
        cv.check_type('spherical harmonics order', order, Integral)
        cv.check_greater_than('spherical harmonics order', order, 0, equality=True)
        self._order = order

    @property
    def cosine(self):
        return self._cosine

    @cosine.setter
    def cosine(self, cosine):
        cv.check_value('Spherical harmonics cosine treatment', cosine,
                       ('scatter', 'particle'))
        self._cosine = cosine

    @classmethod
    def from_hdf5(cls, group, **kwargs):
        if group['type'].value.decode() != cls.short_name.lower():
            raise ValueError("Expected HDF5 data for filter type '"
                             + cls.short_name.lower() + "' but got '"
                             + group['type'].value.decode() + " instead")

        filter_id = int(group.name.split('/')[-1].lstrip('filter '))

        out = cls(group['order'].value, filter_id)
        out.cosine = group['cosine'].value.decode()

        return out

    def to_xml_element(self):
        """Return XML Element representing the filter.

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing spherical harmonics filter data

        """
        element = ET.Element('filter')
        element.set('id', str(self.id))
        element.set('type', self.short_name.lower())
        element.set('cosine', self.cosine)

        subelement = ET.SubElement(element, 'order')
        subelement.text = str(self.order)

        return element


class ZernikeFilter(Filter):
    r"""Score Zernike expansion moments in space up to specified order.

    This filter allows scores to be multiplied by Zernike polynomials of the the
    particle's position along a particular axis, normalized to a given unit
    circle, up to a user-specified order.

    Parameters
    ----------
    order : int
        Maximum Zernike polynomial order
    x : float
        x-coordinate of center of circle for normalization
    y : float
        y-coordinate of center of circle for normalization
    r : int or None
        Radius of circle for normalization

    Attributes
    ----------
    order : int
        Maximum Zernike polynomial order
    x : float
        x-coordinate of center of circle for normalization
    y : float
        y-coordinate of center of circle for normalization
    r : int or None
        Radius of circle for normalization
    id : int
        Unique identifier for the filter
    num_bins : int
        The number of filter bins

    """

    def __init__(self, order, x, y, r, filter_id=None):
        self.order = order
        self.x = x
        self.y = y
        self.r = r
        self.bins = ['Z{},{}'.format(n, m)
                     for n in range(order + 1)
                     for m in range(-n, n + 1, 2)]
        self.id = filter_id

    def __hash__(self):
        string = type(self).__name__ + '\n'
        string += '{: <16}=\t{}\n'.format('\tOrder', self.order)
        return hash(string)

    def __repr__(self):
        string = type(self).__name__ + '\n'
        string += '{: <16}=\t{}\n'.format('\tOrder', self.order)
        string += '{: <16}=\t{}\n'.format('\tID', self.id)
        return string

    @property
    def order(self):
        return self._order

    @order.setter
    def order(self, order):
        cv.check_type('Zernike order', order, Integral)
        cv.check_greater_than('Zernike order', order, 0, equality=True)
        self._order = order

    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, x):
        cv.check_type('x', x, Real)
        self._x = x

    @property
    def y(self):
        return self._y

    @y.setter
    def y(self, y):
        cv.check_type('y', y, Real)
        self._y = y

    @property
    def r(self):
        return self._r

    @r.setter
    def r(self, r):
        cv.check_type('r', r, Real)
        self._r = r

    @classmethod
    def from_hdf5(cls, group, **kwargs):
        if group['type'].value.decode() != cls.short_name.lower():
            raise ValueError("Expected HDF5 data for filter type '"
                             + cls.short_name.lower() + "' but got '"
                             + group['type'].value.decode() + " instead")

        filter_id = int(group.name.split('/')[-1].lstrip('filter '))
        order = group['order'].value
        x, y, r = group['x'].value, group['y'].value, group['r'].value

        return cls(order, x, y, r, filter_id)

    def to_xml_element(self):
        """Return XML Element representing the filter.

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing Zernike filter data

        """
        element = ET.Element('filter')
        element.set('id', str(self.id))
        element.set('type', self.short_name.lower())

        subelement = ET.SubElement(element, 'order')
        subelement.text = str(self.order)
        subelement = ET.SubElement(element, 'x')
        subelement.text = str(self.x)
        subelement = ET.SubElement(element, 'y')
        subelement.text = str(self.y)
        subelement = ET.SubElement(element, 'r')
        subelement.text = str(self.r)

        return element
