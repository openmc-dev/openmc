from collections import Iterable
from numbers import Real, Integral
from xml.etree import ElementTree as ET
import sys

import numpy as np

import openmc
import openmc.mgxs
from openmc.checkvalue import check_type, check_value, check_greater_than, \
    check_iterable_type
from openmc.clean_xml import *

if sys.version_info[0] >= 3:
    basestring = str

# Supported incoming particle MGXS angular treatment representations
_REPRESENTATIONS = ['isotropic', 'angle']


def ndarray_to_string(arr):
    """Converts a numpy ndarray in to a join with spaces between entries
    similar to ' '.join(map(str,arr)) but applied to all sub-dimensions.

    Parameters
    ----------
    arr : numpy.ndarray
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
    energy_groups : openmc.mgxs.EnergyGroups
        Energygroup structure
    representation : {'isotropic', 'angle'}, optional
        Method used in generating the MGXS (isotropic or angle-dependent flux
        weighting). Defaults to 'isotropic'

    Attributes
    ----------
    name : str
        Unique identifier for the xsdata object
    alias : str
        Separate unique identifier for the xsdata object
    zaid : int
        1000*(atomic number) + mass number. As an example, the zaid of U-235
        would be 92235.
    awr : float
        Atomic weight ratio of an isotope.  That is, the ratio of the mass
        of the isotope to the mass of a single neutron.
    kT : float
        Temperature (in units of MeV).
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure
    fissionable : bool
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
    representation : {'isotropic', 'angle'}
        Method used in generating the MGXS (isotropic or angle-dependent flux
        weighting).
    num_azimuthal : int
        Number of equal width angular bins that the azimuthal angular domain is
        subdivided into. This only applies when ``representation`` is "angle".
    num_polar : int
        Number of equal width angular bins that the polar angular domain is
        subdivided into. This only applies when ``representation`` is "angle".
    use_chi : bool
        Whether or not a chi vector or nu-fission matrix was used.
    vector_shape : iterable of int
        Dimensionality of vector multi-group cross sections (e.g., the total
        cross section).  The return result depends on the value of
        ``representation``.
    matrix_shape : iterable of int
        Dimensionality of matrix multi-group cross sections (e.g., the
        fission matrix cross section).  The return result depends on the
        value of ``representation``.
    pn_matrix_shape : iterable of int
        Dimensionality of scattering matrix data (e.g., the
        scattering matrix cross section).  The return result depends on the
        value of ``representation``.
    total : numpy.ndarray
        Group-wise total cross section ordered by increasing group index (i.e.,
        fast to thermal). If ``representation`` is "isotropic", then the length
        of this list should equal the number of groups described in the
        ``groups`` element.  If ``representation`` is "angle", then the length
        of this list should equal the number of groups times the number of
        azimuthal angles times the number of polar angles, with the
        inner-dimension being groups, intermediate-dimension being azimuthal
        angles and outer-dimension being the polar angles.
    absorption : numpy.ndarray
        Group-wise absorption cross section ordered by increasing group index
        (i.e., fast to thermal). If ``representation`` is "isotropic", then the
        length of this list should equal the number of groups described in the
        ``groups`` attribute. If ``representation`` is "angle", then the length
        of this list should equal the number of groups times the number of
        azimuthal angles times the number of polar angles, with the
        inner-dimension being groups, intermediate-dimension being azimuthal
        angles and outer-dimension being the polar angles.
    scatter : numpy.ndarray
        Scattering moment matrices presented with the columns representing
        incoming group and rows representing the outgoing group.  That is,
        down-scatter will be above the diagonal of the resultant matrix.  This
        matrix is repeated for every Legendre order (in order of increasing
        orders) if ``scatt_type`` is "legendre"; otherwise, this matrix is
        repeated for every bin of the histogram or tabular representation.
        Finally, if ``representation`` is "angle", the above is repeated for
        every azimuthal angle and every polar angle, in that order.
    multiplicity : numpy.ndarray
        Ratio of neutrons produced in scattering collisions to the neutrons
        which undergo scattering collisions; that is, the multiplicity provides
        the code with a scaling factor to account for neutrons produced in
        (n,xn) reactions.  This information is assumed isotropic and therefore
        does not need to be repeated for every Legendre moment or
        histogram/tabular bin.  This matrix follows the same arrangement as
        described for the ``scatter`` attribute, with the exception of the data
        needed to provide the scattering type information.
    fission : numpy.ndarray
        Group-wise fission cross section ordered by increasing group index
        (i.e., fast to thermal). If ``representation`` is "isotropic", then the
        length of this list should equal the number of groups described in the
        ``groups`` attribute. If ``representation`` is "angle", then the length
        of this list should equal the number of groups times the number of
        azimuthal angles times the number of polar angles, with the
        inner-dimension being groups, intermediate-dimension being azimuthal
        angles and outer-dimension being the polar angles.
    kappa_fission : numpy.ndarray
        Group-wise kappa-fission cross section ordered by increasing group
        index (i.e., fast to thermal).  If ``representation`` is "isotropic",
        then the length of this list should equal the number of groups in the
        ``groups`` attribute. If ``representation`` is "angle", then the length
        of this list should equal the number of groups times the number of
        azimuthal angles times the number of polar angles, with the
        inner-dimension being groups, intermediate-dimension being azimuthal
        angles and outer-dimension being the polar angles.
    chi : numpy.ndarray
        Group-wise fission spectra ordered by increasing group index (i.e.,
        fast to thermal).  This attribute should be used if making the common
        approximation that the fission spectra does not depend on incoming
        energy.  If the user does not wish to make this approximation, then
        this should not be provided and this information included in the
        ``nu_fission`` element instead.  If ``representation`` is "isotropic",
        then the length of this list should equal the number of groups
        in the ``groups`` element.  If ``representation`` is "angle", then the
        length of this list should equal the number of groups times the number
        of azimuthal angles times the number of polar angles, with the
        inner-dimension being groups, intermediate-dimension being azimuthal
        angles and outer-dimension being the polar angles.
    nu_fission : numpy.ndarray
        Group-wise fission production cross section vector (i.e., if ``chi`` is
        provided), or is the group-wise fission production matrix. If providing
        the vector, it should be ordered the same as the ``fission`` data.  If
        providing the matrix, it should be ordered the same as the
        ``multiplicity`` matrix.

    """

    def __init__(self, name, energy_groups, representation='isotropic'):
        # Initialize class attributes
        self._name = name
        self._energy_groups = energy_groups
        self._representation = representation
        self._alias = None
        self._zaid = None
        self._awr = None
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
        self._kappa_fission = None
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
    def zaid(self):
        return self._zaid

    @property
    def awr(self):
        return self._awr

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
    def use_chi(self):
        return self._use_chi

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
    def kappa_fission(self):
        return self._kappa_fission

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

    @property
    def vector_shape(self):
        if self.representation is 'isotropic':
            return (self.energy_groups.num_groups,)
        elif self.representation is 'angle':
            return (self.num_polar, self.num_azimuthal,
                    self.energy_groups.num_groups)

    @property
    def matrix_shape(self):
        if self.representation is 'isotropic':
            return (self.energy_groups.num_groups,
                    self.energy_groups.num_groups)
        elif self.representation is 'angle':
            return (self.num_polar, self.num_azimuthal,
                    self.energy_groups.num_groups,
                    self.energy_groups.num_groups)

    @property
    def pn_matrix_shape(self):
        if self.representation is 'isotropic':
            return (self.num_orders, self.energy_groups.num_groups,
                    self.energy_groups.num_groups)
        elif self.representation is 'angle':
            return (self.num_polar, self.num_azimuthal, self.num_orders,
                    self.energy_groups.num_groups,
                    self.energy_groups.num_groups)

    @name.setter
    def name(self, name):
        check_type('name for XSdata', name, basestring)
        self._name = name

    @energy_groups.setter
    def energy_groups(self, energy_groups):
        # Check validity of energy_groups
        check_type('energy_groups', energy_groups, openmc.mgxs.EnergyGroups)

        if energy_group.group_edges is None:
            msg = 'Unable to assign an EnergyGroups object ' \
                  'with uninitialized group edges'
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

    @zaid.setter
    def zaid(self, zaid):
        # Check type and value
        check_type('zaid', zaid, Integral)
        check_greater_than('zaid', zaid, 0)
        self._zaid = zaid

    @awr.setter
    def awr(self, awr):
        # Check validity of type and that the awr value is > 0
        check_type('awr', awr, Real)
        check_greater_than('awr', awr, 0.0)
        self._awr = awr

    @kT.setter
    def kT(self, kT):
        # Check validity of type and that the kT value is >= 0
        check_type('kT', kT, Real)
        check_greater_than('kT', kT, 0.0, equality=True)
        self._kT = kT

    @scatt_type.setter
    def scatt_type(self, scatt_type):
        # check to see it is of a valid type and value
        check_value('scatt_type', scatt_type, ['legendre', 'histogram',
                                               'tabular'])
        self._scatt_type = scatt_type

    @order.setter
    def order(self, order):
        # Check type and value
        check_type('order', order, Integral)
        check_greater_than('order', order, 0, equality=True)
        self._order = order

    @tabular_legendre.setter
    def tabular_legendre(self, tabular_legendre):
        # Check to make sure this is a dict and it has our keys with the
        # right values.
        check_type('tabular_legendre', tabular_legendre, dict)
        if 'enable' in tabular_legendre:
            enable = tabular_legendre['enable']
            check_type('enable', enable, bool)
        else:
            msg = 'The tabular_legendre dict must include a value keyed by ' \
                  '"enable"'
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
        check_value('num_polar', num_polar, Integral)
        check_greater_than('num_polar', num_polar, 0)
        self._num_polar = num_polar

    @num_azimuthal.setter
    def num_azimuthal(self, num_azimuthal):
        check_value('num_azimuthal', num_azimuthal, Integral)
        check_greater_than('num_azimuthal', num_azimuthal, 0)
        self._num_azimuthal = num_azimuthal

    @use_chi.setter
    def use_chi(self, use_chi):
        check_type('use_chi', use_chi, bool)
        self._use_chi = use_chi

    @total.setter
    def total(self, total):
        check_type('total', total, Iterable, expected_iter_type=Real)
        # Convert to a numpy array so we can easily get the shape for
        # checking
        nptotal = np.asarray(total)
        check_value('total shape', nptotal.shape, [self.vector_shape])

        self._total = nptotal

    @absorption.setter
    def absorption(self, absorption):
        check_type('absorption', absorption, Iterable, expected_iter_type=Real)
        # Convert to a numpy array so we can easily get the shape for
        # checking
        npabsorption = np.asarray(absorption)
        check_value('absorption shape', npabsorption.shape,
                    [self.vector_shape])

        self._absorption = npabsorption

    @fission.setter
    def fission(self, fission):
        check_type('fission', fission, Iterable, expected_iter_type=Real)
        # Convert to a numpy array so we can easily get the shape for
        # checking
        npfission = np.asarray(fission)
        check_value('fission shape', npfission.shape, [self.vector_shape])

        self._fission = npfission

        if np.sum(self._fission) > 0.0:
            self._fissionable = True

    @kappa_fission.setter
    def kappa_fission(self, kappa_fission):
        check_type('kappa_fission', kappa_fission, Iterable,
                   expected_iter_type=Real)
        # Convert to a numpy array so we can easily get the shape for
        # checking
        npkappa_fission = np.asarray(kappa_fission)
        check_value('kappa fission shape', npkappa_fission.shape,
                    [self.vector_shape])

        self._kappa_fission = npkappa_fission

        if np.sum(self._kappa_fission) > 0.0:
            self._fissionable = True

    @chi.setter
    def chi(self, chi):
        if self.use_chi is not None:
            if not self.use_chi:
                msg = 'Providing "chi" when "nu-fission" already provided as a' \
                      'matrix'
                raise ValueError(msg)

        check_type('chi', chi, Iterable, expected_iter_type=Real)
        # Convert to a numpy array so we can easily get the shape for
        # checking
        npchi = np.asarray(chi)
        # Check the shape
        if npchi.shape != self.vector_shape:
            msg = 'Provided chi iterable does not have the expected shape.'
            raise ValueError(msg)

        self._chi = npchi

        if self.use_chi is not None:
            self.use_chi = True

    @scatter.setter
    def scatter(self, scatter):
        # Convert to a numpy array so we can easily get the shape for
        # checking
        npscatter = np.asarray(scatter)
        check_iterable_type('scatter', npscatter, Real,
                            max_depth=len(npscatter.shape))
        check_value('scatter shape', npscatter.shape, [self.pn_matrix_shape])

        self._scatter = npscatter

    @multiplicity.setter
    def multiplicity(self, multiplicity):
        # Convert to a numpy array so we can easily get the shape for
        # checking
        npmultiplicity = np.asarray(multiplicity)
        check_iterable_type('multiplicity', npmultiplicity, Real,
                            max_depth=len(npmultiplicity.shape))
        check_value('multiplicity shape', npmultiplicity.shape,
                    [self.matrix_shape])

        self._multiplicity = npmultiplicity

    @nu_fission.setter
    def nu_fission(self, nu_fission):
        # The NuFissionXS class does not have the capability to produce
        # a fission matrix and therefore if this path is pursued, we know
        # chi must be used.
        # nu_fission can be given as a vector or a matrix
        # Vector is used when chi also exists.
        # Matrix is used when chi does not exist.
        # We have to check that the correct form is given, but only if
        # chi already has been set.  If not, we just check that this is OK
        # and set the use_chi flag accordingly

        # Convert to a numpy array so we can easily get the shape for
        # checking
        npnu_fission = np.asarray(nu_fission)

        check_iterable_type('nu_fission', npnu_fission, Real,
                            max_depth=len(npnu_fission.shape))

        if self.use_chi is not None:
            if self.use_chi:
                check_value('nu_fission shape', npnu_fission.shape,
                            [self.vector_shape])
            else:
                check_value('nu_fission shape', npnu_fission.shape,
                            [self.matrix_shape])
        else:
            check_value('nu_fission shape', npnu_fission.shape,
                        [self.vector_shape, self.matrix_shape])
            # Find out if we have a nu-fission matrix or vector
            # and set a flag to allow other methods to check this later.
            if npnu_fission.shape == self.vector_shape:
                self.use_chi = True
            else:
                self.use_chi = False

        self._nu_fission = npnu_fission
        if np.sum(self._nu_fission) > 0.0:
            self._fissionable = True

    def set_total_mgxs(self, total, nuclide='total', xs_type='macro'):
        """This method allows for an openmc.mgxs.TotalXS or
        openmc.mgxs.TransportXS to be used to set the total cross section
        for this XSdata object.

        Parameters
        ----------
        total: openmc.mgxs.TotalXS or openmc.mgxs.TransportXS
            MGXS Object containing the total or transport cross section
            for the domain of interest.
        nuclide : str
            Individual nuclide (or 'total' if obtaining material-wise data)
            to gather data for.  Defaults to 'total'.
        xs_type: {'macro', 'micro'}
            Provide the macro or micro cross section in units of cm^-1 or
            barns. Defaults to 'macro'.

        See also
        --------
        openmc.mgxs.Library.create_mg_library()
        openmc.mgxs.Library.get_xsdata

        """

        check_type('total', total, (openmc.mgxs.TotalXS,
                                    openmc.mgxs.TransportXS))
        check_value('energy_groups', total.energy_groups, [self.energy_groups])
        check_value('domain_type', total.domain_type,
                    ['universe', 'cell', 'material'])

        if self.representation is 'isotropic':
            self._total = total.get_xs(nuclides=nuclide, xs_type=xs_type)
        elif self.representation is 'angle':
            msg = 'Angular-Dependent MGXS have not yet been implemented'
            raise ValueError(msg)

    def set_absorption_mgxs(self, absorption, nuclide='total',
                            xs_type='macro'):
        """This method allows for an openmc.mgxs.AbsorptionXS
        to be used to set the absorption cross section for this XSdata object.

        Parameters
        ----------
        absorption: openmc.mgxs.AbsorptionXS
            MGXS Object containing the absorption cross section
            for the domain of interest.
        nuclide : str
            Individual nuclide (or 'total' if obtaining material-wise data)
            to gather data for.  Defaults to 'total'.
        xs_type: {'macro', 'micro'}
            Provide the macro or micro cross section in units of cm^-1 or
            barns. Defaults to 'macro'.

        See also
        --------
        openmc.mgxs.Library.create_mg_library()
        openmc.mgxs.Library.get_xsdata

        """

        check_type('absorption', absorption, openmc.mgxs.AbsorptionXS)
        check_value('energy_groups', absorption.energy_groups,
                    [self.energy_groups])
        check_value('domain_type', absorption.domain_type,
                    ['universe', 'cell', 'material'])

        if self.representation is 'isotropic':
            self._absorption = absorption.get_xs(nuclides=nuclide,
                                                 xs_type=xs_type)
        elif self.representation is 'angle':
            msg = 'Angular-Dependent MGXS have not yet been implemented'
            raise ValueError(msg)

    def set_fission_mgxs(self, fission, nuclide='total', xs_type='macro'):
        """This method allows for an openmc.mgxs.FissionXS
        to be used to set the fission cross section for this XSdata object.

        Parameters
        ----------
        fission: openmc.mgxs.FissionXS
            MGXS Object containing the fission cross section
            for the domain of interest.
        nuclide : str
            Individual nuclide (or 'total' if obtaining material-wise data)
            to gather data for.  Defaults to 'total'.
        xs_type: {'macro', 'micro'}
            Provide the macro or micro cross section in units of cm^-1 or
            barns. Defaults to 'macro'.

        See also
        --------
        openmc.mgxs.Library.create_mg_library()
        openmc.mgxs.Library.get_xsdata

        """

        check_type('fission', fission, openmc.mgxs.FissionXS)
        check_value('energy_groups', fission.energy_groups,
                    [self.energy_groups])
        check_value('domain_type', fission.domain_type,
                    ['universe', 'cell', 'material'])

        if self.representation is 'isotropic':
            self._fission = fission.get_xs(nuclides=nuclide,
                                           xs_type=xs_type)
        elif self.representation is 'angle':
            msg = 'Angular-Dependent MGXS have not yet been implemented'
            raise ValueError(msg)

    def set_nu_fission_mgxs(self, nu_fission, nuclide='total',
                            xs_type='macro'):
        """This method allows for an openmc.mgxs.NuFissionXS
        to be used to set the nu-fission cross section for this XSdata object.

        Parameters
        ----------
        nu_fission: openmc.mgxs.NuFissionXS
            MGXS Object containing the nu-fission cross section
            for the domain of interest.
        nuclide : str
            Individual nuclide (or 'total' if obtaining material-wise data)
            to gather data for.  Defaults to 'total'.
        xs_type: {'macro', 'micro'}
            Provide the macro or micro cross section in units of cm^-1 or
            barns. Defaults to 'macro'.

        See also
        --------
        openmc.mgxs.Library.create_mg_library()
        openmc.mgxs.Library.get_xsdata

        """

        check_type('nu_fission', nu_fission, (openmc.mgxs.NuFissionXS,
                                              openmc.mgxs.NuFissionMatrixXS))
        check_value('energy_groups', nu_fission.energy_groups,
                    [self.energy_groups])
        check_value('domain_type', nu_fission.domain_type,
                    ['universe', 'cell', 'material'])

        if self.representation is 'isotropic':
            self._nu_fission = nu_fission.get_xs(nuclides=nuclide,
                                                 xs_type=xs_type)
        elif self.representation is 'angle':
            msg = 'Angular-Dependent MGXS have not yet been implemented'
            raise ValueError(msg)

        if isinstance(nu_fission, openmc.mgxs.NuFissionMatrixXS):
            self.use_chi = False
        else:
            self.use_chi = True

        if np.sum(self._nu_fission) > 0.0:
            self._fissionable = True

    def set_kappa_fission_mgxs(self, k_fission, nuclide='total',
                               xs_type='macro'):
        """This method allows for an openmc.mgxs.KappaFissionXS
        to be used to set the kappa-fission cross section for this XSdata
        object.

        Parameters
        ----------
        kappa_fission: openmc.mgxs.KappaFissionXS
            MGXS Object containing the kappa-fission cross section
            for the domain of interest.
        nuclide : str
            Individual nuclide (or 'total' if obtaining material-wise data)
            to gather data for.  Defaults to 'total'.
        xs_type: {'macro', 'micro'}
            Provide the macro or micro cross section in units of cm^-1 or
            barns. Defaults to 'macro'.

        See also
        --------
        openmc.mgxs.Library.create_mg_library()
        openmc.mgxs.Library.get_xsdata

        """

        check_type('k_fission', k_fission, openmc.mgxs.KappaFissionXS)
        check_value('energy_groups', k_fission.energy_groups,
                    [self.energy_groups])
        check_value('domain_type', k_fission.domain_type,
                    ['universe', 'cell', 'material'])

        if self.representation is 'isotropic':
            self._kappa_fission = k_fission.get_xs(nuclides=nuclide,
                                                   xs_type=xs_type)
        elif self.representation is 'angle':
            msg = 'Angular-Dependent MGXS have not yet been implemented'
            raise ValueError(msg)

    def set_chi_mgxs(self, chi, nuclide='total', xs_type='macro'):
        """This method allows for an openmc.mgxs.Chi
        to be used to set chi for this XSdata object.

        Parameters
        ----------
        chi: openmc.mgxs.Chi
            MGXS Object containing chi for the domain of interest.
        nuclide : str
            Individual nuclide (or 'total' if obtaining material-wise data)
            to gather data for.  Defaults to 'total'.
        xs_type: {'macro', 'micro'}
            Provide the macro or micro cross section in units of cm^-1 or
            barns. Defaults to 'macro'.

        See also
        --------
        openmc.mgxs.Library.create_mg_library()
        openmc.mgxs.Library.get_xsdata

        """

        if self.use_chi is not None:
            if not self.use_chi:
                msg = 'Providing chi when nu_fission already provided as a ' \
                      'matrix!'
                raise ValueError(msg)

        check_type('chi', chi, openmc.mgxs.Chi)
        check_value('energy_groups', chi.energy_groups, [self.energy_groups])
        check_value('domain_type', chi.domain_type,
                    ['universe', 'cell', 'material'])

        if self.representation is 'isotropic':
            self._chi = chi.get_xs(nuclides=nuclide,
                                   xs_type=xs_type)
        elif self.representation is 'angle':
            msg = 'Angular-Dependent MGXS have not yet been implemented'
            raise ValueError(msg)

        if self.use_chi is not None:
            self.use_chi = True

    def set_scatter_mgxs(self, scatter, nuclide='total', xs_type='macro'):
        """This method allows for an openmc.mgxs.ScatterMatrixXS
        to be used to set the scatter matrix cross section for this XSdata
        object.  If the XSdata.order attribute has not yet been set, then
        it will be set based on the properties of scatter.

        Parameters
        ----------
        scatter: openmc.mgxs.ScatterMatrixXS
            MGXS Object containing the scatter matrix cross section
            for the domain of interest.
        nuclide : str
            Individual nuclide (or 'total' if obtaining material-wise data)
            to gather data for.  Defaults to 'total'.
        xs_type: {'macro', 'micro'}
            Provide the macro or micro cross section in units of cm^-1 or
            barns. Defaults to 'macro'.

        See also
        --------
        openmc.mgxs.Library.create_mg_library()
        openmc.mgxs.Library.get_xsdata

        """

        check_type('scatter', scatter, openmc.mgxs.ScatterMatrixXS)
        check_value('energy_groups', scatter.energy_groups,
                    [self.energy_groups])
        check_value('domain_type', scatter.domain_type,
                    ['universe', 'cell', 'material'])

        if (self.scatt_type != 'legendre'):
            msg = 'Anisotropic scattering representations other than ' \
                  'Legendre expansions have not yet been implemented in ' \
                  'openmc.mgxs.'
            raise ValueError(msg)

        # If the user has not defined XSdata.order, then we will set
        # the order based on the data within scatter.
        # Otherwise, we will check to see that XSdata.order to match
        # the order of scatter
        if self.order is None:
            self.order = scatter.legendre_order
        else:
            check_value('legendre_order', scatter.legendre_order,
                        [self.order])

        if self.representation is 'isotropic':
            # Get the scattering orders in the outermost dimension
            self._scatter = np.zeros((self.num_orders,
                                      self.energy_groups.num_groups,
                                      self.energy_groups.num_groups))
            for moment in range(self.num_orders):
                self._scatter[moment, :, :] = scatter.get_xs(nuclides=nuclide,
                                                             xs_type=xs_type,
                                                             moment=moment)

        elif self.representation is 'angle':
            msg = 'Angular-Dependent MGXS have not yet been implemented'
            raise ValueError(msg)

    def set_multiplicity_mgxs(self, nuscatter, scatter=None, nuclide='total',
                              xs_type='macro'):
        """This method allows for either the direct use of only an
        openmc.mgxs.MultiplicityMatrixXS OR an openmc.mgxs.NuScatterMatrixXS and
        openmc.mgxs.ScatterMatrixXS to be used to set the scattering
        multiplicity for this XSdata object. Multiplicity,
        in OpenMC parlance, is a factor used to account for the production
        of neutrons introduced by scattering multiplication reactions, i.e.,
        (n,xn) events. In this sense, the multiplication matrix is simply
        defined as the ratio of the nu-scatter and scatter matrices.

        Parameters
        ----------
        nuscatter: {openmc.mgxs.NuScatterMatrixXS,
                    openmc.mgxs.MultiplicityMatrixXS}
            MGXS Object containing the matrix cross section for the domain
            of interest.
        scatter: openmc.mgxs.ScatterMatrixXS
            MGXS Object containing the scattering matrix cross section
            for the domain of interest.
        nuclide : str
            Individual nuclide (or 'total' if obtaining material-wise data)
            to gather data for.  Defaults to 'total'.
        xs_type: {'macro', 'micro'}
            Provide the macro or micro cross section in units of cm^-1 or
            barns. Defaults to 'macro'.

        See also
        --------
        openmc.mgxs.Library.create_mg_library()
        openmc.mgxs.Library.get_xsdata

        """

        check_type('nuscatter', nuscatter, (openmc.mgxs.NuScatterMatrixXS,
                                            openmc.mgxs.MultiplicityMatrixXS))
        check_value('energy_groups', nuscatter.energy_groups,
                    [self.energy_groups])
        check_value('domain_type', nuscatter.domain_type,
                    ['universe', 'cell', 'material'])
        if scatter is not None:
            check_type('scatter', scatter, openmc.mgxs.ScatterMatrixXS)
            if isinstance(nuscatter, openmc.mgxs.MultiplicityMatrixXS):
                msg = 'Either an MultiplicityMatrixXS object must be passed ' \
                      'for "nuscatter" or the "scatter" argument must be ' \
                      'provided.'
                raise ValueError(msg)
            check_value('energy_groups', scatter.energy_groups,
                        [self.energy_groups])
            check_value('domain_type', scatter.domain_type,
                        ['universe', 'cell', 'material'])

        if self.representation is 'isotropic':
            nuscatt = nuscatter.get_xs(nuclides=nuclide,
                                       xs_type=xs_type, moment=0)
            if isinstance(nuscatter, openmc.mgxs.MultiplicityMatrixXS):
                self._multiplicity = nuscatt
            else:
                scatt = scatter.get_xs(nuclides=nuclide,
                                       xs_type=xs_type, moment=0)
                self._multiplicity = np.divide(nuscatt, scatt)
        elif self.representation is 'angle':
            msg = 'Angular-Dependent MGXS have not yet been implemented'
            raise ValueError(msg)
        self._multiplicity = np.nan_to_num(self._multiplicity)

    def _get_xsdata_xml(self):
        element = ET.Element('xsdata')
        element.set('name', self._name)

        if self._alias is not None:
            subelement = ET.SubElement(element, 'alias')
            subelement.text = self.alias

        if self._kT is not None:
            subelement = ET.SubElement(element, 'kT')
            subelement.text = str(self._kT)

        if self._zaid is not None:
            subelement = ET.SubElement(element, 'zaid')
            subelement.text = str(self._zaid)

        if self._awr is not None:
            subelement = ET.SubElement(element, 'awr')
            subelement.text = str(self._awr)

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
            subelement.set('num_points',
                           str(self._tabular_legendre['num_points']))

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

            if self._kappa_fission is not None:
                subelement = ET.SubElement(element, 'k_fission')
                subelement.text = ndarray_to_string(self._kappa_fission)

            if self._nu_fission is not None:
                subelement = ET.SubElement(element, 'nu_fission')
                subelement.text = ndarray_to_string(self._nu_fission)

            if self._chi is not None:
                subelement = ET.SubElement(element, 'chi')
                subelement.text = ndarray_to_string(self._chi)

        return element


class MGXSLibrary(object):
    """Multi-Group Cross Sections file used for an OpenMC simulation.
    Corresponds directly to the MG version of the cross_sections.xml input
    file.

    Attributes
    ----------
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure.
    inverse_velocities : Iterable of Real
        Inverse of velocities, units of sec/cm
    xsdatas : Iterable of openmc.XSdata
        Iterable of multi-Group cross section data objects
    """

    def __init__(self, energy_groups):
        self._xsdatas = []
        self._energy_groups = energy_groups
        self._inverse_velocities = None
        self._cross_sections_file = ET.Element('cross_sections')

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
        check_type('energy groups', energy_groups, openmc.mgxs.EnergyGroups)
        self._energy_groups = energy_groups

    def add_xsdata(self, xsdata):
        """Add an XSdata entry to the file.

        Parameters
        ----------
        xsdata : openmc.XSdata
            MGXS information to add

        """
        if not isinstance(xsdata, XSdata):
            msg = 'Unable to add a non-XSdata "{0}" to the ' \
                  'MGXSLibrary instance'.format(xsdata)
            raise ValueError(msg)
        if xsdata.energy_groups != self._energy_groups:
            msg = 'Energy groups of XSdata do not match that of MGXSLibrary.'
            raise ValueError(msg)

        self._xsdatas.append(xsdata)

    def add_xsdatas(self, xsdatas):
        """Add multiple xsdatas to the file.

        Parameters
        ----------
        xsdatas : tuple or list of openmc.XSdata
            XSdatas to add

        """
        check_iterable_type('xsdatas', xsdatas, XSdata)

        for xsdata in xsdatas:
            self.add_xsdata(xsdata)

    def remove_xsdata(self, xsdata):
        """Remove a xsdata from the file

        Parameters
        ----------
        xsdata : openmc.XSdata
            XSdata to remove

        """

        if not isinstance(xsdata, XSdata):
            msg = 'Unable to remove a non-XSdata "{0}" from the ' \
                  'MGXSLibrary instance'.format(xsdata)
            raise ValueError(msg)

        self._xsdatas.remove(xsdata)

    def _create_groups_subelement(self):
        if self._energy_groups is not None:
            element = ET.SubElement(self._cross_sections_file, 'groups')
            element.text = str(self._energy_groups.num_groups)

    def _create_group_structure_subelement(self):
        if self._energy_groups is not None:
            element = ET.SubElement(self._cross_sections_file,
                                    'group_structure')
            element.text = ' '.join(map(str, self._energy_groups.group_edges))

    def _create_inverse_velocities_subelement(self):
        if self._inverse_velocities is not None:
            element = ET.SubElement(self._cross_sections_file,
                                    'inverse_velocities')
            element.text = ' '.join(map(str, self._inverse_velocities))

    def _create_xsdata_subelements(self):
        for xsdata in self._xsdatas:
            xml_element = xsdata._get_xsdata_xml()
            self._cross_sections_file.append(xml_element)

    def export_to_xml(self, filename='mgxs.xml'):
        """Create an mgxs.xml file that can be used for a
        simulation.

        Parameters
        ----------
        filename : str, optional
            filename of file, default is mgxs.xml

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
                   encoding='utf-8', method='xml')
