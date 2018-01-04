from numbers import Integral, Real
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd

import openmc.checkvalue as cv
from . import Filter


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

    @property
    def num_bins(self):
        return self._order + 1

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
            indicating Legendre orders. The number of rows in the DataFrame is
            the same as the total number of bins in the corresponding tally.

        See also
        --------
        Tally.get_pandas_dataframe(), CrossFilter.get_pandas_dataframe()

        """
        # Initialize Pandas DataFrame
        df = pd.DataFrame()

        bins = np.array(['P{}'.format(i) for i in range(self.order + 1)])
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
