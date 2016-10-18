from collections import Iterable
from numbers import Real, Integral
import sys

import numpy as np
import h5py

import openmc
import openmc.mgxs
from openmc.checkvalue import check_type, check_value, check_greater_than, \
    check_iterable_type

if sys.version_info[0] >= 3:
    basestring = str

# Supported incoming particle MGXS angular treatment representations
_REPRESENTATIONS = ['isotropic', 'angle']
_SCATTER_TYPES = ['tabular', 'legendre', 'histogram']
_SCATTER_SHAPES = ["[Order][G][G']"]


class XSdata(object):
    """A multi-group cross section data set providing all the
    multi-group data necessary for a multi-group OpenMC calculation.

    Parameters
    ----------
    name : str
        Name of the mgxs data set.
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure
    representation : {'isotropic', 'angle'}, optional
        Method used in generating the MGXS (isotropic or angle-dependent flux
        weighting). Defaults to 'isotropic'
    temperatures : Iterable of float
        Temperatures (in units of Kelvin) of the provided datasets.  Defaults
        to a single temperature at 294K.

    Attributes
    ----------
    name : str
        Unique identifier for the xsdata object
    aromic_weight_ratio : float
        Atomic weight ratio of an isotope.  That is, the ratio of the mass
        of the isotope to the mass of a single neutron.
    temperatures : numpy.ndarray
        Temperatures (in units of Kelvin) of the provided datasets.  Defaults
        to a single temperature at 294K.
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure
    fissionable : bool
        Whether or not this is a fissionable data set.
    scatter_format : {'legendre', 'histogram', or 'tabular'}
        Angular distribution representation (legendre, histogram, or tabular)
    scatter_shapes : {"[Order][G][G']"}
        Dimensionality of the scattering and multiplicity matrices
    order : int
        Either the Legendre order, number of bins, or number of points used to
        describe the angular distribution associated with each group-to-group
        transfer probability.
    representation : {'isotropic', 'angle'}
        Method used in generating the MGXS (isotropic or angle-dependent flux
        weighting).
    num_azimuthal : int
        Number of equal width angular bins that the azimuthal angular domain is
        subdivided into. This only applies when :attr:`XSdata.representation`
        is "angle".
    num_polar : int
        Number of equal width angular bins that the polar angular domain is
        subdivided into. This only applies when :attr:`XSdata.representation`
        is "angle".
    use_chi : bool
        Whether or not a chi vector or nu-fission matrix was used.
    vector_shape : iterable of int
        Dimensionality of vector multi-group cross sections (e.g., the total
        cross section).  The return result depends on the value of
        :attr:`XSdata.representation`.
    matrix_shape : iterable of int
        Dimensionality of matrix multi-group cross sections (e.g., the
        fission matrix cross section).  The return result depends on the
        value of :attr:`XSdata.representation`.
    pn_matrix_shape : iterable of int
        Dimensionality of scattering matrix data (e.g., the
        scattering matrix cross section).  The return result depends on the
        value of :attr:`XSdata.representation`.
    total : dict of numpy.ndarray
        Group-wise total cross section.
    absorption : dict of numpy.ndarray
        Group-wise absorption cross section.
    scatter_matrix : dict of numpy.ndarray
        Scattering moment matrices presented with the columns representing
        incoming group and rows representing the outgoing group.  That is,
        down-scatter will be above the diagonal of the resultant matrix.
    multiplicity_matrix : dict of numpy.ndarray
        Ratio of neutrons produced in scattering collisions to the neutrons
        which undergo scattering collisions; that is, the multiplicity provides
        the code with a scaling factor to account for neutrons produced in
        (n,xn) reactions.
    fission : dict of numpy.ndarray
        Group-wise fission cross section.
    kappa_fission : dict of numpy.ndarray
        Group-wise kappa-fission cross section.
    chi : dict of numpy.ndarray
        Group-wise fission spectra ordered by increasing group index (i.e.,
        fast to thermal). This attribute should be used if making the common
        approximation that the fission spectra does not depend on incoming
        energy. If the user does not wish to make this approximation, then
        this should not be provided and this information included in the
        :attr:`XSdata.nu_fission` attribute instead.
    nu_fission : dict of numpy.ndarray
        Group-wise fission production cross section vector (i.e., if
        :attr:`XSdata.chi` is provided), or is the group-wise fission
        production matrix.
    inverse_velocities : dict of numpy.ndarray
        Inverse of velocities, in units of sec/cm.

    Notes
    -----
    The parameters containing cross section data have dimensionalities which
    depend upon the value of :attr:`XSdata.representation` as well as the
    number of Legendre or other angular dimensions as described by
    :attr:`XSdata.order`. The :attr:`XSdata.vector_shape`,
    :attr:`XSdata.matrix_shape`, and :attr:`XSdata.pn_matrix_shape` properties
    are provided to obtain the dimensionality of the data for each temperature.

    The following are cross sections which should use each of these properties:

    vector_shape :
        total, absorption, fission, kappa_fission, chi, nu_fission (if chi is
        provided), inverse_velocities
    matrix_shape :
        multiplicity_matrix, nu_fission (if chi is not provided)
    pn_matrix_shape :
        scatter_matrix

    """

    def __init__(self, name, energy_groups, temperatures=[294.],
                 representation='isotropic'):
        # Initialize class attributes
        self.name = name
        self.energy_groups = energy_groups
        self.temperatures = temperatures
        self.representation = representation
        self._atomic_weight_ratio = None
        self._fissionable = False
        self._scatter_format = 'legendre'
        self._scatter_shape = "[Order][G][G']"
        self._order = None
        self._num_polar = None
        self._num_azimuthal = None
        self._use_chi = None
        self._total = len(temperatures) * [None]
        self._absorption = len(temperatures) * [None]
        self._scatter_matrix = len(temperatures) * [None]
        self._multiplicity_matrix = len(temperatures) * [None]
        self._fission = len(temperatures) * [None]
        self._nu_fission = len(temperatures) * [None]
        self._kappa_fission = len(temperatures) * [None]
        self._chi = len(temperatures) * [None]
        self._inverse_velocities = len(temperatures) * [None]

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
    def atomic_weight_ratio(self):
        return self._atomic_weight_ratio

    @property
    def fissionable(self):
        return self._fissionable

    @property
    def temperatures(self):
        return self._temperatures

    @property
    def scatter_format(self):
        return self._scatter_format

    @property
    def scatter_shape(self):
        return self._scatter_shape

    @property
    def order(self):
        return self._order

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
    def scatter_matrix(self):
        return self._scatter_matrix

    @property
    def multiplicity_matrix(self):
        return self._multiplicity_matrix

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
        if self._order is not None:
            if self._scatter_format in (None, 'legendre'):
                return self._order + 1
            else:
                return self._order

    @property
    def vector_shape(self):
        if self.representation == 'isotropic':
            return (self.energy_groups.num_groups,)
        elif self.representation == 'angle':
            return (self.num_polar, self.num_azimuthal,
                    self.energy_groups.num_groups)

    @property
    def matrix_shape(self):
        if self.representation == 'isotropic':
            return (self.energy_groups.num_groups,
                    self.energy_groups.num_groups)
        elif self.representation == 'angle':
            return (self.num_polar, self.num_azimuthal,
                    self.energy_groups.num_groups,
                    self.energy_groups.num_groups)

    @property
    def pn_matrix_shape(self):
        if self.representation == 'isotropic':
            return (self.num_orders, self.energy_groups.num_groups,
                    self.energy_groups.num_groups)
        elif self.representation == 'angle':
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

        if energy_groups.group_edges is None:
            msg = 'Unable to assign an EnergyGroups object ' \
                  'with uninitialized group edges'
            raise ValueError(msg)
        self._energy_groups = energy_groups

    @representation.setter
    def representation(self, representation):
        # Check it is of valid type.
        check_value('representation', representation, _REPRESENTATIONS)
        self._representation = representation

    @atomic_weight_ratio.setter
    def atomic_weight_ratio(self, atomic_weight_ratio):
        # Check validity of type and that the atomic_weight_ratio value is > 0
        check_type('atomic_weight_ratio', atomic_weight_ratio, Real)
        check_greater_than('atomic_weight_ratio', atomic_weight_ratio, 0.0)
        self._atomic_weight_ratio = atomic_weight_ratio

    @temperatures.setter
    def temperatures(self, temperatures):
        check_iterable_type('temperatures', temperatures, Real)

        self._temperatures = np.array(temperatures)

    @scatter_format.setter
    def scatter_format(self, scatter_format):
        # check to see it is of a valid type and value
        check_value('scatter_format', scatter_format, _SCATTER_TYPES)
        self._scatter_format = scatter_format

    @scatter_shape.setter
    def scatter_shape(self, scatter_shape):
        # check to see it is of a valid type and value
        check_value('scatter_shape', scatter_shape, _SCATTER_SHAPES)
        self._scatter_shape = scatter_shape

    @order.setter
    def order(self, order):
        # Check type and value
        check_type('order', order, Integral)
        check_greater_than('order', order, 0, equality=True)
        self._order = order

    @num_polar.setter
    def num_polar(self, num_polar):
        # Make sure we have positive ints
        check_type('num_polar', num_polar, Integral)
        check_greater_than('num_polar', num_polar, 0)
        self._num_polar = num_polar

    @num_azimuthal.setter
    def num_azimuthal(self, num_azimuthal):
        check_type('num_azimuthal', num_azimuthal, Integral)
        check_greater_than('num_azimuthal', num_azimuthal, 0)
        self._num_azimuthal = num_azimuthal

    @use_chi.setter
    def use_chi(self, use_chi):
        check_type('use_chi', use_chi, bool)
        self._use_chi = use_chi

    def add_temperature(self, temperature):
        """This method re-sizes the attributes of this XSdata object so that it
        can accomodate an additional temperature.  Note that the set_* methods
        will still need to be executed.

        Parameters
        ----------
        temperature : float
            Temperature (in units of Kelvin) of the provided dataset.

        """
        check_type('temperature', temperature, Real)

        temp_store = self.temperatures.tolist().append(temperature)
        self.temperatures = temp_store

        self._total.append(None)
        self._absorption.append(None)
        self._scatter_matrix.append(None)
        self._multiplicity_matrix.append(None)
        self._fission.append(None)
        self._nu_fission.append(None)
        self._kappa_fission.append(None)
        self._chi.append(None)
        self._inverse_velocities.append(None)

    def set_total(self, total, temperature=294.):
        """This method sets the cross section for this XSdata object at the
        provided temperature.

        Parameters
        ----------
        total: np.ndarray
            Total Cross Section
        temperature : float
            Temperature (in Kelvin) of the data. Defaults to room temperature
            (294K).

        See also
        --------
        openmc.mgxs_library.set_total_mgxs()

        """
        check_type('total', total, Iterable, expected_iter_type=Real)
        # Convert to a numpy array so we can easily get the shape for checking
        total = np.asarray(total)
        check_value('total shape', total.shape, [self.vector_shape])
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        self._total[i] = total

    def set_absorption(self, absorption, temperature=294.):
        """This method sets the cross section for this XSdata object at the
        provided temperature.

        Parameters
        ----------
        absorption: np.ndarray
            Absorption Cross Section
        temperature : float
            Temperature (in Kelvin) of the data. Defaults to room temperature
            (294K).

        See also
        --------
        openmc.mgxs_library.set_absorption_mgxs()

        """
        check_type('absorption', absorption, Iterable, expected_iter_type=Real)
        # Convert to a numpy array so we can easily get the shape for checking
        absorption = np.asarray(absorption)
        check_value('absorption shape', absorption.shape, [self.vector_shape])
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        self._absorption[i] = absorption

    def set_fission(self, fission, temperature=294.):
        """This method sets the cross section for this XSdata object at the
        provided temperature.

        Parameters
        ----------
        fission: np.ndarray
            Fission Cross Section
        temperature : float
            Temperature (in Kelvin) of the data. Defaults to room temperature
            (294K).

        See also
        --------
        openmc.mgxs_library.set_fission_mgxs()

        """
        check_type('fission', fission, Iterable, expected_iter_type=Real)
        # Convert to a numpy array so we can easily get the shape for checking
        fission = np.asarray(fission)
        check_value('fission shape', fission.shape, [self.vector_shape])
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        self._fission[i] = fission

        if np.sum(fission) > 0.0:
            self._fissionable = True

    def set_kappa_fission(self, kappa_fission, temperature=294.):
        """This method sets the cross section for this XSdata object at the
        provided temperature.

        Parameters
        ----------
        kappa_fission: np.ndarray
            Kappa-Fission Cross Section
        temperature : float
            Temperature (in Kelvin) of the data. Defaults to room temperature
            (294K).

        See also
        --------
        openmc.mgxs_library.set_kappa_fission_mgxs()

        """
        check_type('kappa_fission', kappa_fission, Iterable,
                   expected_iter_type=Real)
        # Convert to a numpy array so we can easily get the shape for checking
        kappa_fission = np.asarray(kappa_fission)
        check_value('kappa fission shape', kappa_fission.shape,
                    [self.vector_shape])
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        self._kappa_fission[i] = kappa_fission

        if np.sum(kappa_fission) > 0.0:
            self._fissionable = True

    def set_chi(self, chi, temperature=294.):
        """This method sets the cross section for this XSdata object at the
        provided temperature.

        Parameters
        ----------
        chi: np.ndarray
            Fission Spectrum
        temperature : float
            Temperature (in Kelvin) of the data. Defaults to room temperature
            (294K).

        See also
        --------
        openmc.mgxs_library.set_chi_mgxs()

        """
        if self.use_chi is not None:
            if not self.use_chi:
                msg = 'Providing "chi" when "nu-fission" already provided ' \
                      'as a matrix'
                raise ValueError(msg)

        check_type('chi', chi, Iterable, expected_iter_type=Real)
        # Convert to a numpy array so we can easily get the shape for checking
        chi = np.asarray(chi)
        # Check the shape
        if chi.shape != self.vector_shape:
            msg = 'Provided chi iterable does not have the expected shape.'
            raise ValueError(msg)
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        self._chi[i] = chi

        if self.use_chi is not None:
            self.use_chi = True

    def set_scatter_matrix(self, scatter, temperature=294.):
        """This method sets the cross section for this XSdata object at the
        provided temperature.

        Parameters
        ----------
        scatter: np.ndarray
            Scattering Matrix Cross Section
        temperature : float
            Temperature (in Kelvin) of the data. Defaults to room temperature
            (294K).

        See also
        --------
        openmc.mgxs_library.set_scatter_matrix_mgxs()

        """
        # Convert to a numpy array so we can easily get the shape for checking
        scatter = np.asarray(scatter)
        check_iterable_type('scatter', scatter, Real,
                            max_depth=len(scatter.shape))
        check_value('scatter shape', scatter.shape, [self.pn_matrix_shape])
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        self._scatter_matrix[i] = scatter

    def set_multiplicity_matrix(self, multiplicity, temperature=294.):
        """This method sets the cross section for this XSdata object at the
        provided temperature.

        Parameters
        ----------
        multiplicity: np.ndarray
            Multiplicity Matrix Cross Section
        temperature : float
            Temperature (in Kelvin) of the data. Defaults to room temperature
            (294K).

        See also
        --------
        openmc.mgxs_library.set_multiplicity_matrix_mgxs()

        """
        # Convert to a numpy array so we can easily get the shape for checking
        multiplicity = np.asarray(multiplicity)
        check_iterable_type('multiplicity', multiplicity, Real,
                            max_depth=len(multiplicity.shape))
        check_value('multiplicity shape', multiplicity.shape,
                    [self.matrix_shape])
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        self._multiplicity_matrix[i] = multiplicity

    def set_nu_fission(self, nu_fission, temperature=294.):
        """This method sets the cross section for this XSdata object at the
        provided temperature.

        Parameters
        ----------
        nu_fission: np.ndarray
            Nu-fission Cross Section
        temperature : float
            Temperature (in Kelvin) of the data. Defaults to room temperature
            (294K).

        See also
        --------
        openmc.mgxs_library.set_nu_fission_mgxs()

        """
        # The NuFissionXS class does not have the capability to produce
        # a fission matrix and therefore if this path is pursued, we know
        # chi must be used.
        # nu_fission can be given as a vector or a matrix
        # Vector is used when chi also exists.
        # Matrix is used when chi does not exist.
        # We have to check that the correct form is given, but only if
        # chi already has been set.  If not, we just check that this is OK
        # and set the use_chi flag accordingly

        # Convert to a numpy array so we can easily get the shape for checking
        nu_fission = np.asarray(nu_fission)

        check_iterable_type('nu_fission', nu_fission, Real,
                            max_depth=len(nu_fission.shape))
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        if self.use_chi is not None:
            if self.use_chi:
                check_value('nu_fission shape', nu_fission.shape,
                            [self.vector_shape])
            else:
                check_value('nu_fission shape', nu_fission.shape,
                            [self.matrix_shape])
        else:
            check_value('nu_fission shape', nu_fission.shape,
                        [self.vector_shape, self.matrix_shape])
            # Find out if we have a nu-fission matrix or vector
            # and set a flag to allow other methods to check this later.
            if nu_fission.shape == self.vector_shape:
                self.use_chi = True
            else:
                self.use_chi = False

        i = np.where(self.temperatures == temperature)[0][0]
        self._nu_fission[i] = nu_fission
        if np.sum(nu_fission) > 0.0:
            self._fissionable = True

    def set_inverse_velocities(self, inv_vel, temperature=294.):
        """This method sets the inverse velocities for this XSdata object at the
        provided temperature.

        Parameters
        ----------
        inv_vel: np.ndarray
            Inverse velocities in units of sec/cm.
        temperature : float
            Temperature (in Kelvin) of the data. Defaults to room temperature
            (294K).

        """
        check_type('inverse velocities', inv_vel, Iterable,
                   expected_iter_type=Real)
        # Convert to a numpy array so we can easily get the shape for checking
        inv_vel = np.asarray(inv_vel)
        check_value('inverse velocities shape', inv_vel.shape,
                    [self.vector_shape])
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        self._inverse_velocities[i] = inv_vel

    def set_total_mgxs(self, total, temperature=294., nuclide='total',
                       xs_type='macro', subdomain=None):
        """This method allows for an openmc.mgxs.TotalXS or
        openmc.mgxs.TransportXS to be used to set the total cross section
        for this XSdata object.

        Parameters
        ----------
        total: openmc.mgxs.TotalXS or openmc.mgxs.TransportXS
            MGXS Object containing the total or transport cross section
            for the domain of interest.
        temperature : float
            Temperature (in Kelvin) of the data. Defaults to room temperature
            (294K).
        nuclide : str
            Individual nuclide (or 'total' if obtaining material-wise data)
            to gather data for.  Defaults to 'total'.
        xs_type: {'macro', 'micro'}
            Provide the macro or micro cross section in units of cm^-1 or
            barns. Defaults to 'macro'.
        subdomain : iterable of int
            If the MGXS contains a mesh domain type, the subdomain parameter
            specifies which mesh cell (i.e., [i, j, k] index) to use.

        See also
        --------
        openmc.mgxs.Library.create_mg_library()
        openmc.mgxs.Library.get_xsdata

        """

        check_type('total', total, (openmc.mgxs.TotalXS,
                                    openmc.mgxs.TransportXS))
        check_value('energy_groups', total.energy_groups, [self.energy_groups])
        check_value('domain_type', total.domain_type,
                    ['universe', 'cell', 'material', 'mesh'])
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        if self.representation == 'isotropic':
            self._total[i] = total.get_xs(nuclides=nuclide, xs_type=xs_type,
                                          subdomains=subdomain)
        elif self.representation == 'angle':
            msg = 'Angular-Dependent MGXS have not yet been implemented'
            raise ValueError(msg)

    def set_absorption_mgxs(self, absorption, temperature=294.,
                            nuclide='total', xs_type='macro', subdomain=None):
        """This method allows for an openmc.mgxs.AbsorptionXS
        to be used to set the absorption cross section for this XSdata object.

        Parameters
        ----------
        absorption: openmc.mgxs.AbsorptionXS
            MGXS Object containing the absorption cross section
            for the domain of interest.
        temperature : float
            Temperature (in Kelvin) of the data. Defaults to room temperature
            (294K).
        nuclide : str
            Individual nuclide (or 'total' if obtaining material-wise data)
            to gather data for.  Defaults to 'total'.
        xs_type: {'macro', 'micro'}
            Provide the macro or micro cross section in units of cm^-1 or
            barns. Defaults to 'macro'.
        subdomain : iterable of int
            If the MGXS contains a mesh domain type, the subdomain parameter
            specifies which mesh cell (i.e., [i, j, k] index) to use.

        See also
        --------
        openmc.mgxs.Library.create_mg_library()
        openmc.mgxs.Library.get_xsdata

        """

        check_type('absorption', absorption, openmc.mgxs.AbsorptionXS)
        check_value('energy_groups', absorption.energy_groups,
                    [self.energy_groups])
        check_value('domain_type', absorption.domain_type,
                    ['universe', 'cell', 'material', 'mesh'])
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        if self.representation == 'isotropic':
            self._absorption[i] = absorption.get_xs(nuclides=nuclide,
                                                    xs_type=xs_type,
                                                    subdomains=subdomain)
        elif self.representation == 'angle':
            msg = 'Angular-Dependent MGXS have not yet been implemented'
            raise ValueError(msg)

    def set_fission_mgxs(self, fission, temperature=294., nuclide='total',
                         xs_type='macro', subdomain=None):
        """This method allows for an openmc.mgxs.FissionXS
        to be used to set the fission cross section for this XSdata object.

        Parameters
        ----------
        fission: openmc.mgxs.FissionXS
            MGXS Object containing the fission cross section
            for the domain of interest.
        temperature : float
            Temperature (in Kelvin) of the data. Defaults to room temperature
            (294K).
        nuclide : str
            Individual nuclide (or 'total' if obtaining material-wise data)
            to gather data for.  Defaults to 'total'.
        xs_type: {'macro', 'micro'}
            Provide the macro or micro cross section in units of cm^-1 or
            barns. Defaults to 'macro'.
        subdomain : iterable of int
            If the MGXS contains a mesh domain type, the subdomain parameter
            specifies which mesh cell (i.e., [i, j, k] index) to use.

        See also
        --------
        openmc.mgxs.Library.create_mg_library()
        openmc.mgxs.Library.get_xsdata

        """

        check_type('fission', fission, openmc.mgxs.FissionXS)
        check_value('energy_groups', fission.energy_groups,
                    [self.energy_groups])
        check_value('domain_type', fission.domain_type,
                    ['universe', 'cell', 'material', 'mesh'])
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        if self.representation == 'isotropic':
            self._fission[i] = fission.get_xs(nuclides=nuclide,
                                              xs_type=xs_type,
                                              subdomains=subdomain)
        elif self.representation == 'angle':
            msg = 'Angular-Dependent MGXS have not yet been implemented'
            raise ValueError(msg)

    def set_nu_fission_mgxs(self, nu_fission, temperature=294.,
                            nuclide='total', xs_type='macro', subdomain=None):
        """This method allows for an openmc.mgxs.NuFissionXS
        to be used to set the nu-fission cross section for this XSdata object.

        Parameters
        ----------
        nu_fission: openmc.mgxs.NuFissionXS
            MGXS Object containing the nu-fission cross section
            for the domain of interest.
        temperature : float
            Temperature (in Kelvin) of the data. Defaults to room temperature
            (294K).
        nuclide : str
            Individual nuclide (or 'total' if obtaining material-wise data)
            to gather data for.  Defaults to 'total'.
        xs_type: {'macro', 'micro'}
            Provide the macro or micro cross section in units of cm^-1 or
            barns. Defaults to 'macro'.
        subdomain : iterable of int
            If the MGXS contains a mesh domain type, the subdomain parameter
            specifies which mesh cell (i.e., [i, j, k] index) to use.

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
                    ['universe', 'cell', 'material', 'mesh'])
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        if self.representation == 'isotropic':
            self._nu_fission[i] = nu_fission.get_xs(nuclides=nuclide,
                                                    xs_type=xs_type,
                                                    subdomains=subdomain)
        elif self.representation == 'angle':
            msg = 'Angular-Dependent MGXS have not yet been implemented'
            raise ValueError(msg)

        if isinstance(nu_fission, openmc.mgxs.NuFissionMatrixXS):
            self.use_chi = False
        else:
            self.use_chi = True

        if np.sum(self._nu_fission) > 0.0:
            self._fissionable = True

    def set_kappa_fission_mgxs(self, k_fission, temperature=294.,
                               nuclide='total', xs_type='macro',
                               subdomain=None):
        """This method allows for an openmc.mgxs.KappaFissionXS
        to be used to set the kappa-fission cross section for this XSdata
        object.

        Parameters
        ----------
        kappa_fission: openmc.mgxs.KappaFissionXS
            MGXS Object containing the kappa-fission cross section
            for the domain of interest.
        temperature : float
            Temperature (in Kelvin) of the data. Defaults to room temperature
            (294K).
        nuclide : str
            Individual nuclide (or 'total' if obtaining material-wise data)
            to gather data for.  Defaults to 'total'.
        xs_type: {'macro', 'micro'}
            Provide the macro or micro cross section in units of cm^-1 or
            barns. Defaults to 'macro'.
        subdomain : iterable of int
            If the MGXS contains a mesh domain type, the subdomain parameter
            specifies which mesh cell (i.e., [i, j, k] index) to use.

        See also
        --------
        openmc.mgxs.Library.create_mg_library()
        openmc.mgxs.Library.get_xsdata

        """

        check_type('k_fission', k_fission, openmc.mgxs.KappaFissionXS)
        check_value('energy_groups', k_fission.energy_groups,
                    [self.energy_groups])
        check_value('domain_type', k_fission.domain_type,
                    ['universe', 'cell', 'material', 'mesh'])
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        if self.representation == 'isotropic':
            self._kappa_fission[i] = k_fission.get_xs(nuclides=nuclide,
                                                      xs_type=xs_type,
                                                      subdomains=subdomain)
        elif self.representation == 'angle':
            msg = 'Angular-Dependent MGXS have not yet been implemented'
            raise ValueError(msg)

    def set_chi_mgxs(self, chi, temperature=294., nuclide='total',
                     xs_type='macro', subdomain=None):
        """This method allows for an openmc.mgxs.Chi
        to be used to set chi for this XSdata object.

        Parameters
        ----------
        chi: openmc.mgxs.Chi
            MGXS Object containing chi for the domain of interest.
        temperature : float
            Temperature (in Kelvin) of the data. Defaults to room temperature
            (294K).
        nuclide : str
            Individual nuclide (or 'total' if obtaining material-wise data)
            to gather data for.  Defaults to 'total'.
        xs_type: {'macro', 'micro'}
            Provide the macro or micro cross section in units of cm^-1 or
            barns. Defaults to 'macro'.
        subdomain : iterable of int
            If the MGXS contains a mesh domain type, the subdomain parameter
            specifies which mesh cell (i.e., [i, j, k] index) to use.

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
                    ['universe', 'cell', 'material', 'mesh'])
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        if self.representation == 'isotropic':
            self._chi[i] = chi.get_xs(nuclides=nuclide,
                                      xs_type=xs_type, subdomains=subdomain)
        elif self.representation == 'angle':
            msg = 'Angular-Dependent MGXS have not yet been implemented'
            raise ValueError(msg)

        if self.use_chi is not None:
            self.use_chi = True

    def set_scatter_matrix_mgxs(self, scatter, temperature=294.,
                                nuclide='total', xs_type='macro',
                                subdomain=None):
        """This method allows for an openmc.mgxs.ScatterMatrixXS
        to be used to set the scatter matrix cross section for this XSdata
        object.  If the XSdata.order attribute has not yet been set, then
        it will be set based on the properties of scatter.

        Parameters
        ----------
        scatter: openmc.mgxs.ScatterMatrixXS
            MGXS Object containing the scatter matrix cross section
            for the domain of interest.
        temperature : float
            Temperature (in Kelvin) of the data. Defaults to room temperature
            (294K).
        nuclide : str
            Individual nuclide (or 'total' if obtaining material-wise data)
            to gather data for.  Defaults to 'total'.
        xs_type: {'macro', 'micro'}
            Provide the macro or micro cross section in units of cm^-1 or
            barns. Defaults to 'macro'.
        subdomain : iterable of int
            If the MGXS contains a mesh domain type, the subdomain parameter
            specifies which mesh cell (i.e., [i, j, k] index) to use.

        See also
        --------
        openmc.mgxs.Library.create_mg_library()
        openmc.mgxs.Library.get_xsdata

        """

        check_type('scatter', scatter, openmc.mgxs.ScatterMatrixXS)
        check_value('energy_groups', scatter.energy_groups,
                    [self.energy_groups])
        check_value('domain_type', scatter.domain_type,
                    ['universe', 'cell', 'material', 'mesh'])
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        if self.scatter_format != 'legendre':
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

        i = np.where(self.temperatures == temperature)[0][0]
        if self.representation == 'isotropic':
            # Get the scattering orders in the outermost dimension
            self._scatter_matrix[i] = np.zeros((self.num_orders,
                                                self.energy_groups.num_groups,
                                                self.energy_groups.num_groups))
            for moment in range(self.num_orders):
                self._scatter_matrix[i][moment, :, :] = \
                    scatter.get_xs(nuclides=nuclide, xs_type=xs_type,
                                   moment=moment, subdomains=subdomain)

        elif self.representation == 'angle':
            msg = 'Angular-Dependent MGXS have not yet been implemented'
            raise ValueError(msg)

    def set_multiplicity_matrix_mgxs(self, nuscatter, scatter=None,
                                     temperature=294., nuclide='total',
                                     xs_type='macro', subdomain=None):
        """This method allows for either the direct use of only an
        openmc.mgxs.MultiplicityMatrixXS OR
        an openmc.mgxs.NuScatterMatrixXS and
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
        temperature : float
            Temperature (in Kelvin) of the data. Defaults to room temperature
            (294K).
        nuclide : str
            Individual nuclide (or 'total' if obtaining material-wise data)
            to gather data for.  Defaults to 'total'.
        xs_type: {'macro', 'micro'}
            Provide the macro or micro cross section in units of cm^-1 or
            barns. Defaults to 'macro'.
        subdomain : iterable of int
            If the MGXS contains a mesh domain type, the subdomain parameter
            specifies which mesh cell (i.e., [i, j, k] index) to use.

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
                    ['universe', 'cell', 'material', 'mesh'])
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

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
                        ['universe', 'cell', 'material', 'mesh'])

        i = np.where(self.temperatures == temperature)[0][0]
        if self.representation == 'isotropic':
            nuscatt = nuscatter.get_xs(nuclides=nuclide,
                                       xs_type=xs_type, moment=0,
                                       subdomains=subdomain)
            if isinstance(nuscatter, openmc.mgxs.MultiplicityMatrixXS):
                self._multiplicity_matrix[i] = nuscatt
            else:
                scatt = scatter.get_xs(nuclides=nuclide,
                                       xs_type=xs_type, moment=0,
                                       subdomains=subdomain)
                self._multiplicity_matrix[i] = np.divide(nuscatt, scatt)
        elif self.representation == 'angle':
            msg = 'Angular-Dependent MGXS have not yet been implemented'
            raise ValueError(msg)
        self._multiplicity_matrix[i] = \
            np.nan_to_num(self._multiplicity_matrix[i])

    def to_hdf5(self, file):
        """Write XSdata to an HDF5 file

        Parameters
        ----------
        file : h5py.File
            HDF5 File (a root Group) to write to

        """
        grp = file.create_group(self.name)
        if self.atomic_weight_ratio is not None:
            grp.attrs['atomic_weight_ratio'] = self.atomic_weight_ratio
        if self.fissionable is not None:
            grp.attrs['fissionable'] = self.fissionable
        if self.representation is not None:
            grp.attrs['representation'] = np.string_(self.representation)
            if self.representation == 'angle':
                if self.num_azimuthal is not None:
                    grp.attrs['num_azimuthal'] = self.num_azimuthal
                if self.num_polar is not None:
                    grp.attrs['num_polar'] = self.num_polar
        if self.scatter_format is not None:
            grp.attrs['scatter_format'] = np.string_(self.scatter_format)
        if self.scatter_shape is not None:
            grp.attrs['scatter_shape'] = np.string_(self.scatter_shape)
        if self.order is not None:
            grp.attrs['order'] = self.order

        ktg = grp.create_group('kTs')
        for temperature in self.temperatures:
            temp_label = str(int(np.round(temperature))) + "K"
            kT = temperature * openmc.data.K_BOLTZMANN
            ktg.create_dataset(temp_label, data=kT)

        # Create the temperature datasets
        for i, temperature in enumerate(self.temperatures):
            xs_grp = grp.create_group(str(int(np.round(temperature))) + "K")
            if self._total[i] is None:
                raise ValueError('total data must be provided when writing '
                                 'the HDF5 library')
            xs_grp.create_dataset("total", data=self._total[i])
            if self._absorption[i] is None:
                raise ValueError('absorption data must be provided when '
                                 'writing the HDF5 library')
            xs_grp.create_dataset("absorption", data=self._absorption[i])
            if self.fissionable:
                if self._fission[i] is not None:
                    xs_grp.create_dataset("fission", data=self._fission[i])
                if self._kappa_fission[i] is not None:
                    xs_grp.create_dataset("kappa-fission",
                                          data=self._kappa_fission[i])
                if self._chi[i] is not None:
                    xs_grp.create_dataset("chi", data=self._chi[i])
                if self._nu_fission[i] is None:
                    raise ValueError('nu-fission data must be provided when '
                                     'writing the HDF5 library')
                xs_grp.create_dataset("nu-fission", data=self._nu_fission[i])

            if self._scatter_matrix[i] is None:
                raise ValueError('Scatter matrix must be provided when '
                                 'writing the HDF5 library')

            # Get the sparse scattering data to print to the library
            G = self.energy_groups.num_groups
            if self.representation == 'isotropic':
                g_out_bounds = np.zeros((G, 2), dtype=np.int)
                for g_in in range(G):
                    nz = np.nonzero(self._scatter_matrix[i][0, g_in, :])
                    g_out_bounds[g_in, 0] = nz[0][0]
                    g_out_bounds[g_in, 1] = nz[0][-1]
                # Now create the flattened scatter matrix array
                matrix = self._scatter_matrix[i]
                flat_scatt = []
                for g_in in range(G):
                    for g_out in range(g_out_bounds[g_in, 0],
                                       g_out_bounds[g_in, 1] + 1):
                        for l in range(len(matrix[:, g_in, g_out])):
                            flat_scatt.append(matrix[l, g_in, g_out])
                # And write it.
                scatt_grp = xs_grp.create_group('scatter_data')
                scatt_grp.create_dataset("scatter_matrix",
                                         data=np.array(flat_scatt))
                # Repeat for multiplicity
                if self._multiplicity_matrix[i] is not None:
                    # Now create the flattened scatter matrix array
                    matrix = self._multiplicity_matrix[i][:, :]
                    flat_mult = []
                    for g_in in range(G):
                        for g_out in range(g_out_bounds[g_in, 0],
                                           g_out_bounds[g_in, 1] + 1):
                            flat_mult.append(matrix[g_in, g_out])
                    scatt_grp.create_dataset("multiplicity_matrix",
                                             data=np.array(flat_mult))

                # And finally, adjust g_out_bounds for 1-based group counting
                # and write it.
                g_out_bounds[:, :] += 1
                scatt_grp.create_dataset("g_min", data=g_out_bounds[:, 0])
                scatt_grp.create_dataset("g_max", data=g_out_bounds[:, 1])

            elif self.representation == 'angle':
                Np = self.num_polar
                Na = self.num_azimuthal
                g_out_bounds = np.zeros((Np, Na, G, 2), dtype=np.int)
                for p in range(Np):
                    for a in range(Na):
                        for g_in in range(G):
                            matrix = self._scatter_matrix[i][p, a, 0, g_in, :]
                            nz = np.nonzero(matrix)
                            g_out_bounds[p, a, g_in, 0] = nz[0][0]
                            g_out_bounds[p, a, g_in, 1] = nz[0][-1]
                # Now create the flattened scatter matrix array
                flat_scatt = []
                for p in range(Np):
                    for a in range(Na):
                        matrix = self._scatter_matrix[i][p, a, :, :, :]
                        for g_in in range(G):
                            for g_out in range(g_out_bounds[p, a, g_in, 0],
                                               g_out_bounds[p, a, g_in, 1] + 1):
                                for l in range(len(matrix[:, g_in, g_out])):
                                    flat_scatt.append(matrix[l, g_in, g_out])
                # And write it.
                scatt_grp = xs_grp.create_group('scatter_data')
                scatt_grp.create_dataset("scatter_matrix",
                                         data=np.array(flat_scatt))
                # Repeat for multiplicity
                if self._multiplicity_matrix[i] is not None:
                    # Now create the flattened scatter matrix array
                    flat_mult = []
                    for p in range(Np):
                        for a in range(Na):
                            matrix = self._multiplicity_matrix[i][p, a, :, :]
                            for g_in in range(G):
                                for g_out in range(g_out_bounds[p, a, g_in, 0],
                                                   g_out_bounds[p, a, g_in, 1] + 1):
                                    flat_mult.append(matrix[g_in, g_out])
                # And write it.
                scatt_grp.create_dataset("multiplicity_matrix",
                                         data=np.array(flat_mult))

                # And finally, adjust g_out_bounds for 1-based group counting
                # and write it.
                g_out_bounds[:, :, :, :] += 1
                scatt_grp.create_dataset("g_min", data=g_out_bounds[:, :, :, 0])
                scatt_grp.create_dataset("g_max", data=g_out_bounds[:, :, :, 1])

            # Add the kinetics data
            if self._inverse_velocities[i] is not None:
                xs_grp.create_dataset("inverse_velocities",
                                      data=self._inverse_velocities[i])


class MGXSLibrary(object):
    """Multi-Group Cross Sections file used for an OpenMC simulation.
    Corresponds directly to the MG version of the cross_sections.xml input
    file.

    Parameters
    ----------
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure

    Attributes
    ----------
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure.
    xsdatas : Iterable of openmc.XSdata
        Iterable of multi-Group cross section data objects
    """

    def __init__(self, energy_groups):
        self.energy_groups = energy_groups
        self._xsdatas = []

    @property
    def energy_groups(self):
        return self._energy_groups

    @property
    def temperatures(self):
        return self._temperatures

    @property
    def xsdatas(self):
        return self._xsdatas

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

    def export_to_hdf5(self, filename='mgxs.h5'):
        """Create an hdf5 file that can be used for a simulation.

        Parameters
        ----------
        filename : str
            Filename of file, default is mgxs.h5.

        """

        check_type('filename', filename, basestring)

        # Create and write to the HDF5 file
        file = h5py.File(filename, "w")
        file.attrs['groups'] = self.energy_groups.num_groups
        file.attrs['group structure'] = self.energy_groups.group_edges

        for xsdata in self._xsdatas:
            xsdata.to_hdf5(file)

        file.close()
