from numbers import Integral, Real
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd

import openmc.checkvalue as cv
from . import Filter


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

    @property
    def num_bins(self):
        n = self._order
        return ((n + 1)*(n + 2))//2

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

    def get_pandas_dataframe(self, data_size, stride, **kwargs):
        """Builds a Pandas DataFrame for the Filter's bins.

        This method constructs a Pandas DataFrame object for the filter with
        columns annotated by filter bin information. This is a helper method for
        :meth:`Tally.get_pandas_dataframe`.

        Parameters
        ----------
        data_size : Integral
            The total number of bins in the tally corresponding to this filter

        Returns
        -------
        pandas.DataFrame
            A Pandas DataFrame with a column that is filled with strings
            indicating Zernike orders. The number of rows in the DataFrame is
            the same as the total number of bins in the corresponding tally.

        See also
        --------
        Tally.get_pandas_dataframe(), CrossFilter.get_pandas_dataframe()

        """
        # Initialize Pandas DataFrame
        df = pd.DataFrame()

        # Create list of strings for each order
        orders = []
        for n in range(self.order + 1):
            for m in range(-n, n + 1, 2):
                orders.append('Z{},{}'.format(n, m))

        bins = np.array(orders)
        filter_bins = np.repeat(bins, stride)
        tile_factor = data_size // len(filter_bins)
        filter_bins = np.tile(filter_bins, tile_factor)
        df = pd.concat([df, pd.DataFrame(
            {self.short_name.lower(): filter_bins})])

        return df

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
