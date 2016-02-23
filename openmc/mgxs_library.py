from collections import Iterable
from numbers import Real, Integral
from xml.etree import ElementTree as ET
import warnings
import sys
if sys.version_info[0] >= 3:
    basestring = str

import numpy as np

import openmc
from openmc.mgxs import EnergyGroups
from openmc.checkvalue import check_type, check_value, check_greater_than, \
                              check_iterable_type
from openmc.clean_xml import *

# Supported incoming particle MGXS angular treatment representations
_REPRESENTATIONS = ['isotropic', 'angle']


def ndarray_to_string(arr):
    """Converts a numpy ndarray in to a join with spaces between entries
    similar to ' '.join(map(str,arr)) but applied to all sub-dimensions.

    Parameters
    ----------
    arr : ndarray
        Array to combine in to a string

    Returns
    -------
    text : str
        String representation of array in arr

    """

    shape = arr.shape
    ndim = arr.ndim
    tab = '    '
    indent = '\n' + tab + tab
    text = indent

    if ndim == 1:
        text += tab
        for i in range(shape[0]):
            text += '{:.7E} '.format(arr[i])
        text += indent
    elif ndim == 2:
        for i in range(shape[0]):
            text += tab
            for j in range(shape[1]):
                text += '{:.7E} '.format(arr[i, j])
            text += indent
    elif ndim == 3:
        for i in range(shape[0]):
            for j in range(shape[1]):
                text += tab
                for k in range(shape[2]):
                    text += '{:.7E} '.format(arr[i, j, k])
                text += indent
    elif ndim == 4:
        for i in range(shape[0]):
            for j in range(shape[1]):
                for k in range(shape[2]):
                    text += tab
                    for l in range(shape[3]):
                        text += '{:.7E} '.format(arr[i, j, k, l])
                    text += indent
    elif ndim == 5:
        for i in range(shape[0]):
            for j in range(shape[1]):
                for k in range(shape[2]):
                    for l in range(shape[3]):
                        text += tab
                        for m in range(shape[4]):
                            text += '{:.7E} '.format(arr[i, j, k, l, m])
                        text += indent

    return text


class XSdata(object):
    """A multi-group cross section data set providing all the
    multi-group data necessary for a multi-group OpenMC calculation.

    Parameters
    ----------
    name : str, optional
        Name of the mgxs data set.

    representation : {'isotropic', 'angle'}
        Method used in generating the MGXS (isotropic or angle-dependent flux
        weighting). Defaults to 'isotropic'

    Attributes
    ----------
    name : str
        Unique identifier for the xsdata object
    alias : str
        Separate unique identifier for the xsdata object
    kT : float
        Temperature (in units of MeV) of this data set.
    energy_groups : EnergyGroups
        Energy group structure
    fissionable : boolean
        Whether or not this is a fissionable data set.
    scatt_type : {'legendre', 'histogram', or 'tabular'}
        Angular distribution representation (legendre, histogram, or tabular)
    order : int
        Either the Legendre order, number of bins, or number of points used to
        describe the angular distribution associated with each group-to-group
        transfer probability.
    tabular_legendre : dict
        Set how to treat the Legendre scattering kernel (tabular or leave in
        Legendre polynomial form). Dict contains two keys: 'enable' and
        'num_points'.  'enable' is a boolean and 'num_points' is the
        number of points to use, if 'enable' is True.

    """
    def __init__(self, name, energy_groups, representation="isotropic"):
        # Initialize class attributes
        self._name = name
        self._energy_groups = energy_groups
        self._representation = representation
        self._alias = None
        self._kT = None
        self._fissionable = False
        self._scatt_type = 'legendre'
        self._order = None
        self._tabular_legendre = None
        self._num_polar = None
        self._num_azimuthal = None
        self._total = None
        self._absorption = None
        self._scatter = None
        self._multiplicity = None
        self._fission = None
        self._nu_fission = None
        self._k_fission = None
        self._chi = None
        self._use_chi = None

    @property
    def name(self):
        return self._name

    @property
    def energy_groups(self):
        return self._energy_groups

    @property
    def representation(self):
        return self._representation

    @property
    def alias(self):
        return self._alias

    @property
    def kT(self):
        return self._kT

    @property
    def scatt_type(self):
        return self._scatt_type

    @property
    def order(self):
        return self._order

    @property
    def tabular_legendre(self):
        return self._tabular_legendre

    @property
    def num_polar(self):
        return self._num_polar

    @property
    def num_azimuthal(self):
        return self._num_azimuthal

    @property
    def total(self):
        return self._total

    @property
    def absorption(self):
        return self._absorption

    @property
    def scatter(self):
        return self._scatter

    @property
    def multiplicity(self):
        return self._multiplicity

    @property
    def fission(self):
        return self._fission

    @property
    def nu_fission(self):
        return self._nu_fission

    @property
    def k_fission(self):
        return self._k_fission

    @property
    def chi(self):
        return self._chi

    @property
    def num_orders(self):
        if (self._order is not None) and (self._scatt_type is not None):
            if self._scatt_type is 'legendre':
                return self._order + 1
            else:
                return self._order

    @name.setter
    def name(self, name):
        check_type('name for XSdata', name, basestring)
        self._name = name

    @energy_groups.setter
    def energy_groups(self, energy_groups):
        # Check validity of energy_groups
        check_type("energy_groups", energy_groups, EnergyGroups)

        # Check that there is one or more groups
        if ((energy_groups.num_groups is None) or
            (energy_groups.num_groups < 1)):

            msg = 'energy_groups object incorrectly initialized.'
            raise ValueError(msg)

        self._energy_groups = energy_groups

    @representation.setter
    def representation(self, representation):
        # Check it is of valid type.
        check_value('representation', representation, _REPRESENTATIONS)
        self._representation = representation

    @alias.setter
    def alias(self, alias):
        if alias is not None:
            check_type('alias', alias, basestring)
            self._alias = alias
        else:
            self._alias = self._name

    @kT.setter
    def kT(self, kT):
        # Check validity of type and that the kT value is >= 0
        check_type("kT", kT, Real)
        check_greater_than("kT", kT, 0.0, equality=True)
        self._kT = kT

    @scatt_type.setter
    def scatt_type(self, scatt_type):
        # check to see it is of a valid type and value
        check_value("scatt_type", scatt_type, ['legendre', 'histogram',
                                               'tabular'])
        self._scatt_type = scatt_type

    @order.setter
    def order(self, order):
        # Check type and value
        check_type("order", order, Integral)
        check_greater_than("order", order, 0, equality=True)
        self._order = order

    @tabular_legendre.setter
    def tabular_legendre(self, tabular_legendre):
        # Check to make sure this is a dict and it has our keys with the
        # right values.
        check_type("tabular_legendre", tabular_legendre, dict)
        if 'enable' in tabular_legendre:
            enable = tabular_legendre['enable']
            check_type('enable', enable, bool)
        else:
            msg = "enable must be provided in tabular_legendre"
            raise ValueError(msg)
        if 'num_points' in tabular_legendre:
            num_points = tabular_legendre['num_points']
            check_value('num_points', num_points, Integral)
            check_greater_than('num_points', num_points, 0)
        else:
            if not enable:
                num_points = 1
            else:
                num_points = 33
        self._tabular_legendre = {'enable': enable, 'num_points': num_points}

    @num_polar.setter
    def num_polar(self, num_polar):
        # Make sure we have positive ints
        check_value("num_polar", num_polar, Integral)
        check_greater_than("num_polar", num_polar, 0)
        self._num_polar = num_polar

    @num_azimuthal.setter
    def num_azimuthal(self, num_azimuthal):
        check_value("num_azimuthal", num_azimuthal, Integral)
        check_greater_than("num_azimuthal", num_azimuthal, 0)
        self._num_azimuthal = num_azimuthal

    @total.setter
    def total(self, total):
        if self._representation is 'isotropic':
            shape = (self._energy_groups.num_groups,)
        elif self._representation is 'angle':
            shape = (self._num_polar, self._num_azimuthal,
                     self._energy_groups.num_groups)
        # check we have a numpy list
        check_type("total", total, np.ndarray, expected_iter_type=Real)
        if total.shape == shape:
            self._total = np.copy(total)
        else:
            msg = 'Shape of provided total "{0}" does not match shape ' \
                  'required, "{1}"'.format(total.shape, shape)
            raise ValueError(msg)

    @absorption.setter
    def absorption(self, absorption):
        if self._representation is 'isotropic':
            shape = (self._energy_groups.num_groups,)
        elif self._representation is 'angle':
            shape = (self._num_polar, self._num_azimuthal,
                     self._energy_groups.num_groups)
        # check we have a numpy list
        check_type("absorption", absorption, np.ndarray, expected_iter_type=Real)
        if absorption.shape == shape:
            self._absorption = np.copy(absorption)
        else:
            msg = 'Shape of provided absorption "{0}" does not match shape ' \
                  'required, "{1}"'.format(absorption.shape, shape)
            raise ValueError(msg)

    @fission.setter
    def fission(self, fission):
        if self._representation is 'isotropic':
            shape = (self._energy_groups.num_groups,)
        elif self._representation is 'angle':
            shape = (self._num_polar, self._num_azimuthal,
                     self._energy_groups.num_groups)
        # check we have a numpy list
        check_type("fission", fission, np.ndarray, expected_iter_type=Real)
        if fission.shape == shape:
            self._fission = np.copy(fission)
            if np.sum(self._fission) > 0.0:
                self._fissionable = True
        else:
            msg = 'Shape of provided fission "{0}" does not match shape ' \
                  'required, "{1}"'.format(fission.shape, shape)
            raise ValueError(msg)

    @k_fission.setter
    def k_fission(self, k_fission):
        if self._representation is 'isotropic':
            shape = (self._energy_groups.num_groups,)
        elif self._representation is 'angle':
            shape = (self._num_polar, self._num_azimuthal,
                     self._energy_groups.num_groups)
        # check we have a numpy list
        check_type("k_fission", k_fission, np.ndarray, expected_iter_type=Real)
        if k_fission.shape == shape:
            self._k_fission = np.copy(k_fission)
            if np.sum(self._k_fission) > 0.0:
                self._fissionable = True
        else:
            msg = 'Shape of provided k_fission "{0}" does not match shape ' \
                  'required, "{1}"'.format(k_fission.shape, shape)
            raise ValueError(msg)

    @chi.setter
    def chi(self, chi):
        if not self._use_chi:
            msg = 'Providing chi when nu_fission already provided as matrix!'
            raise ValueError(msg)
        if self._representation is 'isotropic':
            shape = (self._energy_groups.num_groups,)
        elif self._representation is 'angle':
            shape = (self._num_polar, self._num_azimuthal,
                     self._energy_groups.num_groups)
        # check we have a numpy list
        check_type("chi", chi, np.ndarray, expected_iter_type=Real)
        if chi.shape == shape:
            self._chi = np.copy(chi)
        else:
            msg = 'Shape of provided chi "{0}" does not match shape ' \
                  'required, "{1}"'.format(chi.shape, shape)
            raise ValueError(msg)
        if self._use_chi is not None:
            self._use_chi = True

    @scatter.setter
    def scatter(self, scatter):
        if self._representation is 'isotropic':
            shape = (self.num_orders, self._energy_groups.num_groups,
                     self._energy_groups.num_groups)
            max_depth = 3
        elif self._representation is 'angle':
            shape = (self._num_polar, self._num_azimuthal, self.num_orders,
                     self._energy_groups.num_groups,
                     self._energy_groups.num_groups)
            max_depth = 5
        # check we have a numpy list
        check_iterable_type("scatter", scatter, expected_type=Real,
                            max_depth=max_depth)
        if scatter.shape == shape:
            self._scatter = np.copy(scatter)
        else:
            msg = 'Shape of provided scatter "{0}" does not match shape ' \
                  'required, "{1}"'.format(scatter.shape, shape)
            raise ValueError(msg)

    @multiplicity.setter
    def multiplicity(self, multiplicity):
        if self._representation is 'isotropic':
            shape = (self._energy_groups.num_groups,
                     self._energy_groups.num_groups)
            max_depth = 2
        elif self._representation is 'angle':
            shape = (self._num_polar, self._num_azimuthal,
                     self._energy_groups.num_groups,
                     self._energy_groups.num_groups)
            max_depth = 4
        # check we have a numpy list
        check_iterable_type("multiplicity", multiplicity, expected_type=Real,
                            max_depth=max_depth)
        if multiplicity.shape == shape:
            self._multiplicity = np.copy(multiplicity)
        else:
            msg = 'Shape of provided multiplicity "{0}" does not match shape ' \
                  'required, "{1}"'.format(multiplicity.shape, shape)
            raise ValueError(msg)

    @nu_fission.setter
    def nu_fission(self, nu_fission):
        # nu_fission can be given as a vector or a matrix
        # Vector is used when chi also exists.
        # Matrix is used when chi does not exist.
        # We have to check that the correct form is given, but only if
        # chi already has been set.  If not, we just check that this is OK
        # and set the use_chi flag.

        # First lets set our dimensions here since they get used repeatedly
        # throughout this code.
        if self._representation is 'isotropic':
            shape_vec = (self._energy_groups.num_groups,)
            shape_mat = (self._energy_groups.num_groups,
                         self._energy_groups.num_groups)
        elif self._representation is 'angle':
            shape_vec = (self._num_polar, self._num_azimuthal,
                         self._energy_groups.num_groups)
            shape_mat = (self._num_polar, self._num_azimuthal,
                         self._energy_groups.num_groups,
                         self._energy_groups.num_groups)

        # Begin by checking the case when chi has already been given and thus
        # the rules for filling in nu_fission are set.
        if self._use_chi is not None:
            if self._use_chi:
                shape = shape_vec
            else:
                shape = shape_mat
            if nu_fission.shape != shape:
                msg = "Invalid Shape of Nu_fission!"
                raise ValueError(msg)
        else:
            # Get shape of nu_fission so we can figure if we need chi or not
            if nu_fission.shape == shape_vec:
                self._use_chi = True
                shape = shape_vec
            elif nu_fission.shape == shape_mat:
                self._use_chi = False
                shape = shape_mat
            else:
                msg = "Invalid Shape of Nu_fission!"
                raise ValueError(msg)

        # check we have a numpy list
        check_type("nu_fission", nu_fission, np.ndarray, expected_iter_type=Real)
        self._nu_fission = np.copy(nu_fission)

    def _get_xsdata_xml(self):
        element = ET.Element("xsdata")
        element.set("name", self._name)

        if self._alias is not None:
            subelement = ET.SubElement(element, 'alias')
            subelement.text = self.alias

        if self._kT is not None:
            subelement = ET.SubElement(element, 'kT')
            subelement.text = str(self._kT)

        if self._fissionable is not None:
            subelement = ET.SubElement(element, 'fissionable')
            subelement.text = str(self._fissionable)

        if self._representation is not None:
            subelement = ET.SubElement(element, 'representation')
            subelement.text = self._representation

        if self._representation == 'angle':
            if self._num_azimuthal is not None:
                subelement = ET.SubElement(element, 'num_azimuthal')
                subelement.text = str(self._num_azimuthal)
            if self._num_polar is not None:
                subelement = ET.SubElement(element, 'num_polar')
                subelement.text = str(self._num_polar)

        if self._scatt_type is not None:
            subelement = ET.SubElement(element, 'scatt_type')
            subelement.text = self._scatt_type

        if self._order is not None:
            subelement = ET.SubElement(element, 'order')
            subelement.text = str(self._order)

        if self._tabular_legendre is not None:
            subelement = ET.SubElement(element, 'tabular_legendre')
            subelement.set('enable', str(self._tabular_legendre['enable']))
            subelement.set('num_points', str(self._tabular_legendre['num_points']))

        if self._total is not None:
            subelement = ET.SubElement(element, 'total')
            subelement.text = ndarray_to_string(self._total)

        if self._absorption is not None:
            subelement = ET.SubElement(element, 'absorption')
            subelement.text = ndarray_to_string(self._absorption)

        if self._scatter is not None:
            subelement = ET.SubElement(element, 'scatter')
            subelement.text = ndarray_to_string(self._scatter)

        if self._multiplicity is not None:
            subelement = ET.SubElement(element, 'multiplicity')
            subelement.text = ndarray_to_string(self._multiplicity)

        if self._fissionable:
            if self._fission is not None:
                subelement = ET.SubElement(element, 'fission')
                subelement.text = ndarray_to_string(self._fission)

            if self._k_fission is not None:
                subelement = ET.SubElement(element, 'k_fission')
                subelement.text = ndarray_to_string(self._k_fission)

            if self._nu_fission is not None:
                subelement = ET.SubElement(element, 'nu_fission')
                subelement.text = ndarray_to_string(self._nu_fission)

            if self._chi is not None:
                subelement = ET.SubElement(element, 'chi')
                subelement.text = ndarray_to_string(self._chi)

        return element


class MGXSLibraryFile(object):
    """Multi-Group Cross Sections file used for an OpenMC simulation.
    Corresponds directly to the MG version of the cross_sections.xml input file.

    Attributes
    ----------
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure.
    inverse_velocities : Iterable of Real
        Inverse of velocities, units of sec/cm
    filename : str
        XML file to write to.
    xsdatas : Iterable of XSdata
        Iterable of multi-Group cross section data objects
    """

    def __init__(self, energy_groups):
        # Initialize MGXSLibraryFile class attributes
        self._xsdatas = []
        self._energy_groups = energy_groups
        self._inverse_velocities = None
        self._cross_sections_file = ET.Element("cross_sections")

    @property
    def inverse_velocities(self):
        return self._inverse_velocities

    @property
    def energy_groups(self):
        return self._energy_groups

    @inverse_velocities.setter
    def inverse_velocities(self, inverse_velocities):
        cv.check_type('inverse_velocities', inverse_velocities, Iterable, Real)
        cv.check_greater_than('number of inverse_velocities',
                              len(inverse_velocities), 0.0)
        self._inverse_velocities = np.array(inverse_velocities)

    @energy_groups.setter
    def energy_groups(self, energy_groups):
        check_type("energy groups", energy_groups, EnergyGroups)
        self._energy_groups = energy_groups

    def add_xsdata(self, xsdata):
        """Add an XSdata entry to the file.

        Parameters
        ----------
        xsdata : XSdata
            MGXS information to add

        """

        # Check the type
        if not isinstance(xsdata, XSdata):
            msg = 'Unable to add a non-XSdata "{0}" to the ' \
                  'MGXSLibraryFile'.format(xsdata)
            raise ValueError(msg)

        # Make sure energy groups match.
        if xsdata.energy_groups != self._energy_groups:
            msg = 'Energy groups of XSdata do not match that of MGXSLibraryFile!'
            raise ValueError(msg)

        self._xsdatas.append(xsdata)

    def add_xsdatas(self, xsdatas):
        """Add multiple xsdatas to the file.

        Parameters
        ----------
        xsdatas : tuple or list of XSdata
            XSdatas to add

        """

        if not isinstance(xsdatas, Iterable):
            msg = 'Unable to create OpenMC xsdatas.xml file from "{0}" which ' \
                  'is not iterable'.format(xsdatas)
            raise ValueError(msg)

        for xsdata in xsdatas:
            self.add_xsdata(xsdata)

    def remove_xsdata(self, xsdata):
        """Remove a xsdata from the file

        Parameters
        ----------
        xsdata : XSdata
            XSdata to remove

        """

        if not isinstance(xsdata, XSdata):
            msg = 'Unable to remove a non-XSdata "{0}" from the ' \
                  'XSdatasFile'.format(xsdata)
            raise ValueError(msg)

        self._xsdatas.remove(xsdata)

    def _create_groups_subelement(self):
        if self._energy_groups is not None:
            element = ET.SubElement(self._cross_sections_file, "groups")
            element.text = str(self._energy_groups.num_groups)

    def _create_group_structure_subelement(self):
        if self._energy_groups is not None:
            element = ET.SubElement(self._cross_sections_file,
                                    "group_structure")
            element.text = ' '.join(map(str, self._energy_groups.group_edges))

    def _create_inverse_velocities_subelement(self):
        if self._inverse_velocities is not None:
            element = ET.SubElement(self._cross_sections_file,
                                    "inverse_velocities")
            element.text = ' '.join(map(str, self._inverse_velocities))

    def _create_xsdata_subelements(self):
        for xsdata in self._xsdatas:
            xml_element = xsdata._get_xsdata_xml()
            self._cross_sections_file.append(xml_element)

    def export_to_xml(self, filename='mg_cross_sections.xml'):
        """Create an mg_cross_sections.xml file that can be used for a
        simulation.

        Parameters
        ----------
        filename : str, optional
            filename of file, default is mg_cross_sections.xml

        """

        # Reset xml element tree
        self._cross_sections_file.clear()

        self._create_groups_subelement()
        self._create_group_structure_subelement()
        self._create_inverse_velocities_subelement()
        self._create_xsdata_subelements()

        # Clean the indentation in the file to be user-readable
        sort_xml_elements(self._cross_sections_file)
        clean_xml_indentation(self._cross_sections_file)

        # Write the XML Tree to the xsdatas.xml file
        tree = ET.ElementTree(self._cross_sections_file)
        tree.write(filename, xml_declaration=True,
                             encoding='utf-8', method="xml")



