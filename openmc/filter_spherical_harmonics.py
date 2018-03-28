from numbers import Integral
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd

import openmc.checkvalue as cv
from . import Filter


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

    @property
    def num_bins(self):
        return (self._order + 1)**2

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
