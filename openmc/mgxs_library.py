import copy
from numbers import Real, Integral
import os

import numpy as np
import h5py
from scipy.interpolate import interp1d
from scipy.integrate import simps
from scipy.special import eval_legendre

import openmc
import openmc.mgxs
from openmc.checkvalue import check_type, check_value, check_greater_than, \
    check_iterable_type, check_less_than, check_filetype_version

ROOM_TEMPERATURE_KELVIN = 294.0

# Supported incoming particle MGXS angular treatment representations
REPRESENTATION_ISOTROPIC = 'isotropic'
REPRESENTATION_ANGLE = 'angle'
_REPRESENTATIONS = [
    REPRESENTATION_ISOTROPIC,
    REPRESENTATION_ANGLE
]

# Supported scattering angular distribution representations
SCATTER_TABULAR = 'tabular'
SCATTER_LEGENDRE = 'legendre'
SCATTER_HISTOGRAM = 'histogram'
_SCATTER_TYPES = [
    SCATTER_TABULAR,
    SCATTER_LEGENDRE,
    SCATTER_HISTOGRAM
]

# List of MGXS indexing schemes
_XS_SHAPES = ["[G][G'][Order]", "[G]", "[G']", "[G][G']", "[DG]", "[DG][G]",
              "[DG][G']", "[DG][G][G']"]

# Number of mu points for conversion between scattering formats
_NMU = 257

# Filetype name of the MGXS Library
_FILETYPE_MGXS_LIBRARY = 'mgxs'

# Current version of the MGXS Library Format
_VERSION_MGXS_LIBRARY = 1


class XSdata:
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
    num_delayed_groups : int
        Number of delayed groups

    Attributes
    ----------
    name : str
        Unique identifier for the xsdata object
    atomic_weight_ratio : float
        Atomic weight ratio of an isotope.  That is, the ratio of the mass
        of the isotope to the mass of a single neutron.
    temperatures : numpy.ndarray
        Temperatures (in units of Kelvin) of the provided datasets.  Defaults
        to a single temperature at 294K.
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure
    num_delayed_groups : int
        Num delayed groups
    fissionable : bool
        Whether or not this is a fissionable data set.
    scatter_format : {'legendre', 'histogram', or 'tabular'}
        Angular distribution representation (legendre, histogram, or tabular)
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
    total : list of numpy.ndarray
        Group-wise total cross section.
    absorption : list of numpy.ndarray
        Group-wise absorption cross section.
    scatter_matrix : list of numpy.ndarray
        Scattering moment matrices presented with the columns representing
        incoming group and rows representing the outgoing group.  That is,
        down-scatter will be above the diagonal of the resultant matrix.
    multiplicity_matrix : list of numpy.ndarray
        Ratio of neutrons produced in scattering collisions to the neutrons
        which undergo scattering collisions; that is, the multiplicity provides
        the code with a scaling factor to account for neutrons produced in
        (n,xn) reactions.
    fission : list of numpy.ndarray
        Group-wise fission cross section.
    kappa_fission : list of numpy.ndarray
        Group-wise kappa_fission cross section.
    chi : list of numpy.ndarray
        Group-wise fission spectra ordered by increasing group index (i.e.,
        fast to thermal). This attribute should be used if making the common
        approximation that the fission spectra does not depend on incoming
        energy. If the user does not wish to make this approximation, then
        this should not be provided and this information included in the
        :attr:`XSdata.nu_fission` attribute instead.
    chi_prompt : list of numpy.ndarray
        Group-wise prompt fission spectra ordered by increasing group index
        (i.e., fast to thermal). This attribute should be used if chi from
        prompt and delayed neutrons is being set separately.
    chi_delayed : list of numpy.ndarray
        Group-wise delayed fission spectra ordered by increasing group index
        (i.e., fast to thermal). This attribute should be used if chi from
        prompt and delayed neutrons is being set separately.
    nu_fission : list of numpy.ndarray
        Group-wise fission production cross section vector (i.e., if ``chi`` is
        provided), or is the group-wise fission production matrix.
    prompt_nu_fission : list of numpy.ndarray
        Group-wise prompt fission production cross section vector.
    delayed_nu_fission : list of numpy.ndarray
        Group-wise delayed fission production cross section vector.
    beta : list of numpy.ndarray
        Delayed-group-wise delayed neutron fraction cross section vector.
    decay_rate : list of numpy.ndarray
        Delayed-group-wise decay rate vector.
    inverse_velocity : list of numpy.ndarray
        Inverse of velocity, in units of sec/cm.
    xs_shapes : dict of iterable of int
        Dictionary with keys of _XS_SHAPES and iterable of int values with the
        corresponding shapes where "Order" corresponds to the pn scattering
        order, "G" corresponds to incoming energy group, "G'" corresponds to
        outgoing energy group, and "DG" corresponds to delayed group.

    Notes
    -----
    The parameters containing cross section data have dimensionalities which
    depend upon the value of :attr:`XSdata.representation` as well as the
    number of Legendre or other angular dimensions as described by
    :attr:`XSdata.order`. The :attr:`XSdata.xs_shapes` are provided to obtain
    the dimensionality of the data for each temperature.

    The following are cross sections which should use each of the properties.
    Note that some cross sections can be input in more than one shape so they
    are listed multiple times:

    [G][G'][Order]: scatter_matrix

    [G]: total, absorption, fission, kappa_fission, nu_fission,
         prompt_nu_fission, delayed_nu_fission, inverse_velocity

    [G']: chi, chi_prompt, chi_delayed

    [G][G']: multiplicity_matrix, nu_fission, prompt_nu_fission

    [DG]: beta, decay_rate

    [DG][G]: delayed_nu_fission, beta, decay_rate

    [DG][G']: chi_delayed

    [DG][G][G']: delayed_nu_fission

    """

    def __init__(self, name, energy_groups, temperatures=[ROOM_TEMPERATURE_KELVIN],
                 representation=REPRESENTATION_ISOTROPIC, num_delayed_groups=0):

        # Initialize class attributes
        self.name = name
        self.energy_groups = energy_groups
        self.num_delayed_groups = num_delayed_groups
        self.temperatures = temperatures
        self.representation = representation
        self._atomic_weight_ratio = None
        self._fissionable = False
        self._scatter_format = SCATTER_LEGENDRE
        self._order = None
        self._num_polar = None
        self._num_azimuthal = None
        self._total = len(temperatures) * [None]
        self._absorption = len(temperatures) * [None]
        self._scatter_matrix = len(temperatures) * [None]
        self._multiplicity_matrix = len(temperatures) * [None]
        self._fission = len(temperatures) * [None]
        self._nu_fission = len(temperatures) * [None]
        self._prompt_nu_fission = len(temperatures) * [None]
        self._delayed_nu_fission = len(temperatures) * [None]
        self._kappa_fission = len(temperatures) * [None]
        self._chi = len(temperatures) * [None]
        self._chi_prompt = len(temperatures) * [None]
        self._chi_delayed = len(temperatures) * [None]
        self._beta = len(temperatures) * [None]
        self._decay_rate = len(temperatures) * [None]
        self._inverse_velocity = len(temperatures) * [None]
        self._xs_shapes = None

    def __deepcopy__(self, memo):
        existing = memo.get(id(self))

        # If this is the first time we have tried to copy this object, copy it
        if existing is None:
            clone = type(self).__new__(type(self))
            clone._name = self.name
            clone._energy_groups = copy.deepcopy(self.energy_groups, memo)
            clone._num_delayed_groups = self.num_delayed_groups
            clone._temperatures = copy.deepcopy(self.temperatures, memo)
            clone._representation = self.representation
            clone._atomic_weight_ratio = self._atomic_weight_ratio
            clone._fissionable = self._fissionable
            clone._scatter_format = self._scatter_format
            clone._order = self._order
            clone._num_polar = self._num_polar
            clone._num_azimuthal = self._num_azimuthal
            clone._total = copy.deepcopy(self._total, memo)
            clone._absorption = copy.deepcopy(self._absorption, memo)
            clone._scatter_matrix = copy.deepcopy(self._scatter_matrix, memo)
            clone._multiplicity_matrix = \
                copy.deepcopy(self._multiplicity_matrix, memo)
            clone._fission = copy.deepcopy(self._fission, memo)
            clone._nu_fission = copy.deepcopy(self._nu_fission, memo)
            clone._prompt_nu_fission = \
                copy.deepcopy(self._prompt_nu_fission, memo)
            clone._delayed_nu_fission = \
                copy.deepcopy(self._delayed_nu_fission, memo)
            clone._kappa_fission = copy.deepcopy(self._kappa_fission, memo)
            clone._chi = copy.deepcopy(self._chi, memo)
            clone._chi_prompt = copy.deepcopy(self._chi_prompt, memo)
            clone._chi_delayed = copy.deepcopy(self._chi_delayed, memo)
            clone._beta = copy.deepcopy(self._beta, memo)
            clone._decay_rate = copy.deepcopy(self._decay_rate, memo)
            clone._inverse_velocity = \
                copy.deepcopy(self._inverse_velocity, memo)
            clone._xs_shapes = copy.deepcopy(self._xs_shapes, memo)

            memo[id(self)] = clone

            return clone

        # If this object has been copied before, return the first copy made
        else:
            return existing

    @property
    def name(self):
        return self._name

    @property
    def energy_groups(self):
        return self._energy_groups

    @property
    def num_delayed_groups(self):
        return self._num_delayed_groups

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
    def order(self):
        return self._order

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
    def prompt_nu_fission(self):
        return self._prompt_nu_fission

    @property
    def delayed_nu_fission(self):
        return self._delayed_nu_fission

    @property
    def kappa_fission(self):
        return self._kappa_fission

    @property
    def chi(self):
        return self._chi

    @property
    def chi_prompt(self):
        return self._chi_prompt

    @property
    def chi_delayed(self):
        return self._chi_delayed

    @property
    def num_orders(self):
        if self._order is None:
            raise ValueError('Order has not been set.')

        if self._scatter_format in (None, SCATTER_LEGENDRE):
            return self._order + 1
        else:
            return self._order

    @property
    def xs_shapes(self):

        if self._xs_shapes is None:

            self._xs_shapes = {}
            self._xs_shapes["[G]"] = (self.energy_groups.num_groups,)
            self._xs_shapes["[G']"] = (self.energy_groups.num_groups,)
            self._xs_shapes["[G][G']"] = (self.energy_groups.num_groups,
                                          self.energy_groups.num_groups)
            self._xs_shapes["[DG]"] = (self.num_delayed_groups,)
            self._xs_shapes["[DG][G]"] = (self.num_delayed_groups,
                                          self.energy_groups.num_groups)
            self._xs_shapes["[DG][G']"] = (self.num_delayed_groups,
                                           self.energy_groups.num_groups)
            self._xs_shapes["[DG][G][G']"] = (self.num_delayed_groups,
                                              self.energy_groups.num_groups,
                                              self.energy_groups.num_groups)

            self._xs_shapes["[G][G'][Order]"] \
                = (self.energy_groups.num_groups,
                   self.energy_groups.num_groups, self.num_orders)

            # If representation is by angle prepend num polar and num azim
            if self.representation == REPRESENTATION_ANGLE:
                for key, shapes in self._xs_shapes.items():
                    self._xs_shapes[key] \
                        = (self.num_polar, self.num_azimuthal) + shapes

        return self._xs_shapes

    @name.setter
    def name(self, name):

        check_type('name for XSdata', name, str)
        self._name = name

    @energy_groups.setter
    def energy_groups(self, energy_groups):

        check_type('energy_groups', energy_groups, openmc.mgxs.EnergyGroups)
        if energy_groups.group_edges is None:
            msg = 'Unable to assign an EnergyGroups object ' \
                  'with uninitialized group edges'
            raise ValueError(msg)

        self._energy_groups = energy_groups

    @num_delayed_groups.setter
    def num_delayed_groups(self, num_delayed_groups):

        check_type('num_delayed_groups', num_delayed_groups, Integral)
        check_less_than('num_delayed_groups', num_delayed_groups,
                        openmc.mgxs.MAX_DELAYED_GROUPS, equality=True)
        check_greater_than('num_delayed_groups', num_delayed_groups, 0,
                           equality=True)
        self._num_delayed_groups = num_delayed_groups

    @representation.setter
    def representation(self, representation):

        check_value('representation', representation, _REPRESENTATIONS)
        self._representation = representation

    @atomic_weight_ratio.setter
    def atomic_weight_ratio(self, atomic_weight_ratio):

        check_type('atomic_weight_ratio', atomic_weight_ratio, Real)
        check_greater_than('atomic_weight_ratio', atomic_weight_ratio, 0.0)
        self._atomic_weight_ratio = atomic_weight_ratio

    @temperatures.setter
    def temperatures(self, temperatures):

        check_iterable_type('temperatures', temperatures, Real)
        self._temperatures = np.array(temperatures)

    @scatter_format.setter
    def scatter_format(self, scatter_format):

        check_value('scatter_format', scatter_format, _SCATTER_TYPES)
        self._scatter_format = scatter_format

    @order.setter
    def order(self, order):

        check_type('order', order, Integral)
        check_greater_than('order', order, 0, equality=True)
        self._order = order

    @num_polar.setter
    def num_polar(self, num_polar):

        check_type('num_polar', num_polar, Integral)
        check_greater_than('num_polar', num_polar, 0)
        self._num_polar = num_polar

    @num_azimuthal.setter
    def num_azimuthal(self, num_azimuthal):

        check_type('num_azimuthal', num_azimuthal, Integral)
        check_greater_than('num_azimuthal', num_azimuthal, 0)
        self._num_azimuthal = num_azimuthal

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
        self._prompt_nu_fission.append(None)
        self._delayed_nu_fission.append(None)
        self._kappa_fission.append(None)
        self._chi.append(None)
        self._chi_prompt.append(None)
        self._chi_delayed.append(None)
        self._beta.append(None)
        self._decay_rate.append(None)
        self._inverse_velocity.append(None)

    def _check_temperature(self, temperature):
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

    def _temperature_index(self, temperature):
        return np.where(self.temperatures == temperature)[0][0]

    def _set_fissionable(self, array):
        if np.sum(array) > 0:
            self._fissionable = True

    def set_total(self, total, temperature=ROOM_TEMPERATURE_KELVIN):
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

        # Get the accepted shapes for this xs
        shapes = [self.xs_shapes["[G]"]]

        # Convert to a numpy array so we can easily get the shape for checking
        total = np.asarray(total)
        check_value('total shape', total.shape, shapes)
        self._check_temperature(temperature)

        i = self._temperature_index(temperature)
        self._total[i] = total

    def set_absorption(self, absorption, temperature=ROOM_TEMPERATURE_KELVIN):
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

        # Get the accepted shapes for this xs
        shapes = [self.xs_shapes["[G]"]]

        # Convert to a numpy array so we can easily get the shape for checking
        absorption = np.asarray(absorption)
        check_value('absorption shape', absorption.shape, shapes)
        self._check_temperature(temperature)

        i = self._temperature_index(temperature)
        self._absorption[i] = absorption

    def set_fission(self, fission, temperature=ROOM_TEMPERATURE_KELVIN):
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

        # Get the accepted shapes for this xs
        shapes = [self.xs_shapes["[G]"]]

        # Convert to a numpy array so we can easily get the shape for checking
        fission = np.asarray(fission)
        check_value('fission shape', fission.shape, shapes)
        self._check_temperature(temperature)

        i = self._temperature_index(temperature)
        self._fission[i] = fission

        self._set_fissionable(fission)

    def set_kappa_fission(self, kappa_fission, temperature=ROOM_TEMPERATURE_KELVIN):
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

        # Get the accepted shapes for this xs
        shapes = [self.xs_shapes["[G]"]]

        # Convert to a numpy array so we can easily get the shape for checking
        kappa_fission = np.asarray(kappa_fission)
        check_value('kappa fission shape', kappa_fission.shape, shapes)
        self._check_temperature(temperature)

        i = self._temperature_index(temperature)
        self._kappa_fission[i] = kappa_fission

        self._set_fissionable(kappa_fission)

    def set_chi(self, chi, temperature=ROOM_TEMPERATURE_KELVIN):
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

        # Get the accepted shapes for this xs
        shapes = [self.xs_shapes["[G']"]]

        # Convert to a numpy array so we can easily get the shape for checking
        chi = np.asarray(chi)
        check_value('chi shape', chi.shape, shapes)
        self._check_temperature(temperature)

        i = self._temperature_index(temperature)
        self._chi[i] = chi

    def set_chi_prompt(self, chi_prompt, temperature=ROOM_TEMPERATURE_KELVIN):
        """This method sets the cross section for this XSdata object at the
        provided temperature.

        Parameters
        ----------
        chi_prompt : np.ndarray
            Prompt fission Spectrum
        temperature : float
            Temperature (in units of Kelvin) of the provided dataset. Defaults
            to 294K

        See also
        --------
        openmc.mgxs_library.set_chi_prompt_mgxs()

        """

        # Get the accepted shapes for this xs
        shapes = [self.xs_shapes["[G']"]]

        # Convert to a numpy array so we can easily get the shape for checking
        chi_prompt = np.asarray(chi_prompt)
        check_value('chi prompt shape', chi_prompt.shape, shapes)
        self._check_temperature(temperature)

        i = self._temperature_index(temperature)
        self._chi_prompt[i] = chi_prompt

    def set_chi_delayed(self, chi_delayed, temperature=ROOM_TEMPERATURE_KELVIN):
        """This method sets the cross section for this XSdata object at the
        provided temperature.

        Parameters
        ----------
        chi_delayed : np.ndarray
            Delayed fission Spectrum
        temperature : float
            Temperature (in units of Kelvin) of the provided dataset. Defaults
            to 294K

        See also
        --------
        openmc.mgxs_library.set_chi_delayed_mgxs()

        """

        # Get the accepted shapes for this xs
        shapes = [self.xs_shapes["[G']"], self.xs_shapes["[DG][G']"]]

        # Convert to a numpy array so we can easily get the shape for checking
        chi_delayed = np.asarray(chi_delayed)
        check_value('chi delayed shape', chi_delayed.shape, shapes)
        self._check_temperature(temperature)
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = self._temperature_index(temperature)
        self._chi_delayed[i] = chi_delayed

    def set_beta(self, beta, temperature=ROOM_TEMPERATURE_KELVIN):
        """This method sets the cross section for this XSdata object at the
        provided temperature.

        Parameters
        ----------
        beta : np.ndarray
            Delayed fission spectrum
        temperature : float
            Temperature (in units of Kelvin) of the provided dataset. Defaults
            to 294K

        See also
        --------
        openmc.mgxs_library.set_beta_mgxs()

        """

        # Get the accepted shapes for this xs
        shapes = [self.xs_shapes["[DG]"], self.xs_shapes["[DG][G]"]]

        # Convert to a numpy array so we can easily get the shape for checking
        beta = np.asarray(beta)
        check_value('beta shape', beta.shape, shapes)
        self._check_temperature(temperature)

        i = self._temperature_index(temperature)
        self._beta[i] = beta

    def set_decay_rate(self, decay_rate, temperature=ROOM_TEMPERATURE_KELVIN):
        """This method sets the cross section for this XSdata object at the
        provided temperature.

        Parameters
        ----------
        decay_rate : np.ndarray
            Delayed neutron precursor decay rate
        temperature : float
            Temperature (in units of Kelvin) of the provided dataset. Defaults
            to 294K

        See also
        --------
        openmc.mgxs_library.set_decay_rate_mgxs()

        """

        # Get the accepted shapes for this xs
        shapes = [self.xs_shapes["[DG]"], self.xs_shapes["[DG][G]"]]

        # Convert to a numpy array so we can easily get the shape for checking
        decay_rate = np.asarray(decay_rate)
        check_value('decay rate shape', decay_rate.shape, shapes)
        self._check_temperature(temperature)

        i = self._temperature_index(temperature)
        self._decay_rate[i] = decay_rate

    def set_scatter_matrix(self, scatter, temperature=ROOM_TEMPERATURE_KELVIN):
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

        # Get the accepted shapes for this xs
        shapes = [self.xs_shapes["[G][G'][Order]"]]

        # Convert to a numpy array so we can easily get the shape for checking
        scatter = np.asarray(scatter)
        check_iterable_type('scatter', scatter, Real,
                            max_depth=len(scatter.shape))
        check_value('scatter shape', scatter.shape, shapes)
        self._check_temperature(temperature)

        i = self._temperature_index(temperature)
        self._scatter_matrix[i] = scatter

    def set_multiplicity_matrix(self, multiplicity, temperature=ROOM_TEMPERATURE_KELVIN):
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

        # Get the accepted shapes for this xs
        shapes = [self.xs_shapes["[G][G']"]]

        # Convert to a numpy array so we can easily get the shape for checking
        multiplicity = np.asarray(multiplicity)
        check_iterable_type('multiplicity', multiplicity, Real,
                            max_depth=len(multiplicity.shape))
        check_value('multiplicity shape', multiplicity.shape, shapes)
        self._check_temperature(temperature)

        i = self._temperature_index(temperature)
        self._multiplicity_matrix[i] = multiplicity

    def set_nu_fission(self, nu_fission, temperature=ROOM_TEMPERATURE_KELVIN):
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

        # Get the accepted shapes for this xs
        shapes = [self.xs_shapes["[G]"], self.xs_shapes["[G][G']"]]

        # Convert to a numpy array so we can easily get the shape for checking
        nu_fission = np.asarray(nu_fission)
        check_value('nu_fission shape', nu_fission.shape, shapes)
        check_iterable_type('nu_fission', nu_fission, Real,
                            max_depth=len(nu_fission.shape))
        self._check_temperature(temperature)

        i = self._temperature_index(temperature)
        self._nu_fission[i] = nu_fission
        self._set_fissionable(nu_fission)

    def set_prompt_nu_fission(self, prompt_nu_fission, temperature=ROOM_TEMPERATURE_KELVIN):
        """This method sets the cross section for this XSdata object at the
        provided temperature.

        Parameters
        ----------
        prompt_nu_fission: np.ndarray
            Prompt-nu-fission Cross Section
        temperature : float
            Temperature (in units of Kelvin) of the provided dataset. Defaults
            to 294K

        See also
        --------
        openmc.mgxs_library.set_prompt_nu_fission_mgxs()

        """

        # Get the accepted shapes for this xs
        shapes = [self.xs_shapes["[G]"], self.xs_shapes["[G][G']"]]

        # Convert to a numpy array so we can easily get the shape for checking
        prompt_nu_fission = np.asarray(prompt_nu_fission)
        check_value('prompt_nu_fission shape', prompt_nu_fission.shape, shapes)
        check_iterable_type('prompt_nu_fission', prompt_nu_fission, Real,
                            max_depth=len(prompt_nu_fission.shape))
        self._check_temperature(temperature)

        i = self._temperature_index(temperature)
        self._prompt_nu_fission[i] = prompt_nu_fission
        self._set_fissionable(prompt_nu_fission)

    def set_delayed_nu_fission(self, delayed_nu_fission, temperature=ROOM_TEMPERATURE_KELVIN):
        """This method sets the cross section for this XSdata object at the
        provided temperature.

        Parameters
        ----------
        delayed_nu_fission: np.ndarray
            Delayed-nu-fission Cross Section
        temperature : float
            Temperature (in units of Kelvin) of the provided dataset. Defaults
            to 294K

        See also
        --------
        openmc.mgxs_library.set_delayed_nu_fission_mgxs()

        """

        # Get the accepted shapes for this xs
        shapes = [self.xs_shapes["[DG][G]"], self.xs_shapes["[DG][G][G']"]]

        # Convert to a numpy array so we can easily get the shape for checking
        delayed_nu_fission = np.asarray(delayed_nu_fission)
        check_value('delayed_nu_fission shape', delayed_nu_fission.shape,
                    shapes)
        check_iterable_type('delayed_nu_fission', delayed_nu_fission, Real,
                            max_depth=len(delayed_nu_fission.shape))
        self._check_temperature(temperature)

        i = self._temperature_index(temperature)
        self._delayed_nu_fission[i] = delayed_nu_fission
        self._set_fissionable(delayed_nu_fission)

    def set_inverse_velocity(self, inv_vel, temperature=ROOM_TEMPERATURE_KELVIN):
        """This method sets the inverse velocity for this XSdata object at the
        provided temperature.

        Parameters
        ----------
        inv_vel: np.ndarray
            Inverse velocity in units of sec/cm.
        temperature : float
            Temperature (in Kelvin) of the data. Defaults to room temperature
            (294K).

        """

        # Get the accepted shapes for this xs
        shapes = [self.xs_shapes["[G]"]]

        # Convert to a numpy array so we can easily get the shape for checking
        inv_vel = np.asarray(inv_vel)
        check_value('inverse_velocity shape', inv_vel.shape, shapes)
        self._check_temperature(temperature)

        i = self._temperature_index(temperature)
        self._inverse_velocity[i] = inv_vel

    def set_total_mgxs(self, total, temperature=ROOM_TEMPERATURE_KELVIN, nuclide='total',
                       xs_type='macro', subdomain=None):
        """This method allows for an openmc.mgxs.TotalXS or
        openmc.mgxs.TransportXS to be used to set the total cross section for
        this XSdata object.

        Parameters
        ----------
        total: openmc.mgxs.TotalXS or openmc.mgxs.TransportXS
            MGXS Object containing the total, transport or nu-transport cross
            section for the domain of interest.
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
        openmc.mgxs.Library.get_xsdata()

        """

        check_type('total', total, (openmc.mgxs.TotalXS,
                                    openmc.mgxs.TransportXS))
        check_value('energy_groups', total.energy_groups, [self.energy_groups])
        check_value('domain_type', total.domain_type, openmc.mgxs.DOMAIN_TYPES)
        self._check_temperature(temperature)

        i = self._temperature_index(temperature)
        self._total[i] = total.get_xs(nuclides=nuclide, xs_type=xs_type,
                                      subdomains=subdomain)

    def set_absorption_mgxs(self, absorption, temperature=ROOM_TEMPERATURE_KELVIN,
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
        openmc.mgxs.Library.get_xsdata()

        """

        check_type('absorption', absorption, openmc.mgxs.AbsorptionXS)
        check_value('energy_groups', absorption.energy_groups,
                    [self.energy_groups])
        check_value('domain_type', absorption.domain_type,
                    openmc.mgxs.DOMAIN_TYPES)
        self._check_temperature(temperature)

        i = self._temperature_index(temperature)
        self._absorption[i] = absorption.get_xs(nuclides=nuclide,
                                                xs_type=xs_type,
                                                subdomains=subdomain)

    def set_fission_mgxs(self, fission, temperature=ROOM_TEMPERATURE_KELVIN, nuclide='total',
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
        openmc.mgxs.Library.get_xsdata()

        """

        check_type('fission', fission, openmc.mgxs.FissionXS)
        check_value('energy_groups', fission.energy_groups,
                    [self.energy_groups])
        check_value('domain_type', fission.domain_type,
                    openmc.mgxs.DOMAIN_TYPES)
        self._check_temperature(temperature)

        i = self._temperature_index(temperature)
        self._fission[i] = fission.get_xs(nuclides=nuclide,
                                          xs_type=xs_type,
                                          subdomains=subdomain)

    def set_nu_fission_mgxs(self, nu_fission, temperature=ROOM_TEMPERATURE_KELVIN,
                            nuclide='total', xs_type='macro', subdomain=None):
        """This method allows for an openmc.mgxs.FissionXS
        to be used to set the nu-fission cross section for this XSdata object.

        Parameters
        ----------
        nu_fission: openmc.mgxs.FissionXS
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
        openmc.mgxs.Library.get_xsdata()

        """

        check_type('nu_fission', nu_fission, (openmc.mgxs.FissionXS,
                                              openmc.mgxs.NuFissionMatrixXS))
        if isinstance(nu_fission, openmc.mgxs.FissionXS):
            check_value('nu', nu_fission.nu, [True])
        check_value('energy_groups', nu_fission.energy_groups,
                    [self.energy_groups])
        check_value('domain_type', nu_fission.domain_type,
                    openmc.mgxs.DOMAIN_TYPES)
        self._check_temperature(temperature)

        i = self._temperature_index(temperature)
        self._nu_fission[i] = nu_fission.get_xs(nuclides=nuclide,
                                                xs_type=xs_type,
                                                subdomains=subdomain)

        self._set_fissionable(self._nu_fission)

    def set_prompt_nu_fission_mgxs(self, prompt_nu_fission, temperature=ROOM_TEMPERATURE_KELVIN,
                                   nuclide='total', xs_type='macro',
                                   subdomain=None):
        """Sets the prompt-nu-fission cross section.

        This method allows for an openmc.mgxs.FissionXS or
        openmc.mgxs.NuFissionMatrixXS to be used to set the prompt-nu-fission
        cross section for this XSdata object.

        Parameters
        ----------
        prompt_nu_fission: openmc.mgxs.FissionXS or openmc.mgxs.NuFissionMatrixXS
            MGXS Object containing the prompt-nu-fission cross section
            for the domain of interest.
        temperature : float
            Temperature (in units of Kelvin) of the provided dataset. Defaults
            to 294K
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
        openmc.mgxs.Library.get_xsdata()

        """

        check_type('prompt_nu_fission', prompt_nu_fission,
                   (openmc.mgxs.FissionXS, openmc.mgxs.NuFissionMatrixXS))
        check_value('prompt', prompt_nu_fission.prompt, [True])
        check_value('energy_groups', prompt_nu_fission.energy_groups,
                    [self.energy_groups])
        check_value('domain_type', prompt_nu_fission.domain_type,
                    openmc.mgxs.DOMAIN_TYPES)
        self._check_temperature(temperature)

        i = self._temperature_index(temperature)
        self._prompt_nu_fission[i] = prompt_nu_fission.get_xs(
            nuclides=nuclide, xs_type=xs_type, subdomains=subdomain)

        self._set_fissionable(self._prompt_nu_fission)

    def set_delayed_nu_fission_mgxs(self, delayed_nu_fission, temperature=ROOM_TEMPERATURE_KELVIN,
                                    nuclide='total', xs_type='macro',
                                    subdomain=None):
        """This method allows for an openmc.mgxs.DelayedNuFissionXS or
        openmc.mgxs.DelayedNuFissionMatrixXS to be used to set the
        delayed-nu-fission cross section for this XSdata object.

        Parameters
        ----------
        delayed_nu_fission: openmc.mgxs.DelayedNuFissionXS or openmc.mgxs.DelayedNuFissionMatrixXS
            MGXS Object containing the delayed-nu-fission cross section
            for the domain of interest.
        temperature : float
            Temperature (in units of Kelvin) of the provided dataset. Defaults
            to 294K
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
        openmc.mgxs.Library.get_xsdata()

        """

        check_type('delayed_nu_fission', delayed_nu_fission,
                   (openmc.mgxs.DelayedNuFissionXS,
                    openmc.mgxs.DelayedNuFissionMatrixXS))
        check_value('energy_groups', delayed_nu_fission.energy_groups,
                    [self.energy_groups])
        check_value('num_delayed_groups', delayed_nu_fission.num_delayed_groups,
                    [self.num_delayed_groups])
        check_value('domain_type', delayed_nu_fission.domain_type,
                    openmc.mgxs.DOMAIN_TYPES)
        self._check_temperature(temperature)

        i = self._temperature_index(temperature)
        self._delayed_nu_fission[i] = delayed_nu_fission.get_xs(
            nuclides=nuclide, xs_type=xs_type, subdomains=subdomain)

        self._set_fissionable(self._delayed_nu_fission)

    def set_kappa_fission_mgxs(self, k_fission, temperature=ROOM_TEMPERATURE_KELVIN,
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
        openmc.mgxs.Library.get_xsdata()

        """

        check_type('kappa_fission', k_fission, openmc.mgxs.KappaFissionXS)
        check_value('energy_groups', k_fission.energy_groups,
                    [self.energy_groups])
        check_value('domain_type', k_fission.domain_type,
                    openmc.mgxs.DOMAIN_TYPES)
        self._check_temperature(temperature)

        i = self._temperature_index(temperature)
        self._kappa_fission[i] = k_fission.get_xs(nuclides=nuclide,
                                                  xs_type=xs_type,
                                                  subdomains=subdomain)

    def set_chi_mgxs(self, chi, temperature=ROOM_TEMPERATURE_KELVIN, nuclide='total',
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
        openmc.mgxs.Library.get_xsdata()

        """

        check_type('chi', chi, openmc.mgxs.Chi)
        check_value('energy_groups', chi.energy_groups, [self.energy_groups])
        check_value('domain_type', chi.domain_type, openmc.mgxs.DOMAIN_TYPES)
        self._check_temperature(temperature)

        i = self._temperature_index(temperature)
        self._chi[i] = chi.get_xs(nuclides=nuclide, xs_type=xs_type,
                                  subdomains=subdomain)

    def set_chi_prompt_mgxs(self, chi_prompt, temperature=ROOM_TEMPERATURE_KELVIN,
                            nuclide='total', xs_type='macro', subdomain=None):
        """This method allows for an openmc.mgxs.Chi to be used to set
        chi-prompt for this XSdata object.

        Parameters
        ----------
        chi_prompt: openmc.mgxs.Chi
            MGXS Object containing chi-prompt for the domain of interest.
        temperature : float
            Temperature (in units of Kelvin) of the provided dataset. Defaults
            to 294K
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
        openmc.mgxs.Library.get_xsdata()

        """

        check_type('chi_prompt', chi_prompt, openmc.mgxs.Chi)
        check_value('prompt', chi_prompt.prompt, [True])
        check_value('energy_groups', chi_prompt.energy_groups,
                    [self.energy_groups])
        check_value('domain_type', chi_prompt.domain_type,
                    openmc.mgxs.DOMAIN_TYPES)
        self._check_temperature(temperature)

        i = self._temperature_index(temperature)
        self._chi_prompt[i] = chi_prompt.get_xs(nuclides=nuclide,
                                                xs_type=xs_type,
                                                subdomains=subdomain)

    def set_chi_delayed_mgxs(self, chi_delayed, temperature=ROOM_TEMPERATURE_KELVIN,
                             nuclide='total', xs_type='macro', subdomain=None):
        """This method allows for an openmc.mgxs.ChiDelayed
        to be used to set chi-delayed for this XSdata object.

        Parameters
        ----------
        chi_delayed: openmc.mgxs.ChiDelayed
            MGXS Object containing chi-delayed for the domain of interest.
        temperature : float
            Temperature (in units of Kelvin) of the provided dataset. Defaults
            to 294K
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
        openmc.mgxs.Library.get_xsdata()

        """

        check_type('chi_delayed', chi_delayed, openmc.mgxs.ChiDelayed)
        check_value('energy_groups', chi_delayed.energy_groups,
                    [self.energy_groups])
        check_value('num_delayed_groups', chi_delayed.num_delayed_groups,
                    [self.num_delayed_groups])
        check_value('domain_type', chi_delayed.domain_type,
                    openmc.mgxs.DOMAIN_TYPES)
        self._check_temperature(temperature)

        i = self._temperature_index(temperature)
        self._chi_delayed[i] = chi_delayed.get_xs(nuclides=nuclide,
                                                  xs_type=xs_type,
                                                  subdomains=subdomain)

    def set_beta_mgxs(self, beta, temperature=ROOM_TEMPERATURE_KELVIN,
                      nuclide='total', xs_type='macro', subdomain=None):
        """This method allows for an openmc.mgxs.Beta
        to be used to set beta for this XSdata object.

        Parameters
        ----------
        beta : openmc.mgxs.Beta
            MGXS Object containing beta for the domain of interest.
        temperature : float
            Temperature (in units of Kelvin) of the provided dataset. Defaults
            to 294K
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
        openmc.mgxs.Library.get_xsdata()

        """

        check_type('beta', beta, openmc.mgxs.Beta)
        check_value('num_delayed_groups', beta.num_delayed_groups,
                    [self.num_delayed_groups])
        check_value('domain_type', beta.domain_type, openmc.mgxs.DOMAIN_TYPES)
        self._check_temperature(temperature)

        i = self._temperature_index(temperature)
        self._beta[i] = beta.get_xs(nuclides=nuclide,
                                    xs_type=xs_type,
                                    subdomains=subdomain)

    def set_decay_rate_mgxs(self, decay_rate, temperature=ROOM_TEMPERATURE_KELVIN,
                            nuclide='total', xs_type='macro', subdomain=None):
        """This method allows for an openmc.mgxs.DecayRate
        to be used to set decay rate for this XSdata object.

        Parameters
        ----------
        decay_rate : openmc.mgxs.DecayRate
            MGXS Object containing decay rate for the domain of interest.
        temperature : float
            Temperature (in units of Kelvin) of the provided dataset. Defaults
            to 294K
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
        openmc.mgxs.Library.get_xsdata()

        """

        check_type('decay_rate', decay_rate, openmc.mgxs.DecayRate)
        check_value('num_delayed_groups', decay_rate.num_delayed_groups,
                    [self.num_delayed_groups])
        check_value('domain_type', decay_rate.domain_type,
                    openmc.mgxs.DOMAIN_TYPES)
        self._check_temperature(temperature)

        i = self._temperature_index(temperature)
        self._decay_rate[i] = decay_rate.get_xs(nuclides=nuclide,
                                                xs_type=xs_type,
                                                subdomains=subdomain)

    def set_scatter_matrix_mgxs(self, scatter, temperature=ROOM_TEMPERATURE_KELVIN,
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
        openmc.mgxs.Library.get_xsdata()

        """

        check_type('scatter', scatter, openmc.mgxs.ScatterMatrixXS)
        check_value('energy_groups', scatter.energy_groups,
                    [self.energy_groups])
        check_value('domain_type', scatter.domain_type,
                    openmc.mgxs.DOMAIN_TYPES)
        self._check_temperature(temperature)

        # Set the value of scatter_format based on the same value within
        # scatter
        self.scatter_format = scatter.scatter_format

        # If the user has not defined XSdata.order, then we will set
        # the order based on the data within scatter.
        # Otherwise, we will check to see that XSdata.order matches
        # the order of scatter
        if self.scatter_format == SCATTER_LEGENDRE:
            if self.order is None:
                self.order = scatter.legendre_order
            else:
                check_value('legendre_order', scatter.legendre_order,
                            [self.order])
        elif self.scatter_format == SCATTER_HISTOGRAM:
            if self.order is None:
                self.order = scatter.histogram_bins
            else:
                check_value('histogram_bins', scatter.histogram_bins,
                            [self.order])

        i = self._temperature_index(temperature)
        if self.scatter_format == SCATTER_LEGENDRE:
            self._scatter_matrix[i] = \
                np.zeros(self.xs_shapes["[G][G'][Order]"])
            # Get the scattering orders in the outermost dimension
            if self.representation == REPRESENTATION_ISOTROPIC:
                for moment in range(self.num_orders):
                    self._scatter_matrix[i][:, :, moment] = \
                        scatter.get_xs(nuclides=nuclide, xs_type=xs_type,
                                       moment=moment, subdomains=subdomain)
            elif self.representation == REPRESENTATION_ANGLE:
                for moment in range(self.num_orders):
                    self._scatter_matrix[i][:, :, :, :, moment] = \
                        scatter.get_xs(nuclides=nuclide, xs_type=xs_type,
                                       moment=moment, subdomains=subdomain)
        else:
            self._scatter_matrix[i] = \
                scatter.get_xs(nuclides=nuclide, xs_type=xs_type,
                               subdomains=subdomain)

    def set_multiplicity_matrix_mgxs(self, nuscatter, scatter=None,
                                     temperature=ROOM_TEMPERATURE_KELVIN, nuclide='total',
                                     xs_type='macro', subdomain=None):
        """This method allows for either the direct use of only an
        openmc.mgxs.MultiplicityMatrixXS or an openmc.mgxs.ScatterMatrixXS and
        openmc.mgxs.ScatterMatrixXS to be used to set the scattering
        multiplicity for this XSdata object. Multiplicity, in OpenMC parlance,
        is a factor used to account for the production of neutrons introduced by
        scattering multiplication reactions, i.e., (n,xn) events. In this sense,
        the multiplication matrix is simply defined as the ratio of the
        nu-scatter and scatter matrices.

        Parameters
        ----------
        nuscatter: openmc.mgxs.ScatterMatrixXS or openmc.mgxs.MultiplicityMatrixXS
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
        openmc.mgxs.Library.get_xsdata()

        """

        check_type('nuscatter', nuscatter, (openmc.mgxs.ScatterMatrixXS,
                                            openmc.mgxs.MultiplicityMatrixXS))
        check_value('energy_groups', nuscatter.energy_groups,
                    [self.energy_groups])
        check_value('domain_type', nuscatter.domain_type,
                    openmc.mgxs.DOMAIN_TYPES)
        self._check_temperature(temperature)

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
                        openmc.mgxs.DOMAIN_TYPES)
        i = self._temperature_index(temperature)
        nuscatt = nuscatter.get_xs(nuclides=nuclide,
                                   xs_type=xs_type, moment=0,
                                   subdomains=subdomain)
        if isinstance(nuscatter, openmc.mgxs.MultiplicityMatrixXS):
            self._multiplicity_matrix[i] = nuscatt
        else:
            scatt = scatter.get_xs(nuclides=nuclide,
                                   xs_type=xs_type, moment=0,
                                   subdomains=subdomain)
            if scatter.scatter_format == SCATTER_HISTOGRAM:
                scatt = np.sum(scatt, axis=2)
            if nuscatter.scatter_format == SCATTER_HISTOGRAM:
                nuscatt = np.sum(nuscatt, axis=2)
            self._multiplicity_matrix[i] = np.divide(nuscatt, scatt)

        self._multiplicity_matrix[i] = \
            np.nan_to_num(self._multiplicity_matrix[i])

    def set_inverse_velocity_mgxs(self, inverse_velocity, temperature=ROOM_TEMPERATURE_KELVIN,
                                  nuclide='total', xs_type='macro',
                                  subdomain=None):
        """This method allows for an openmc.mgxs.InverseVelocity
        to be used to set the inverse velocity for this XSdata object.

        Parameters
        ----------
        inverse_velocity : openmc.mgxs.InverseVelocity
            MGXS object containing the inverse velocity for the domain of
            interest.
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
        openmc.mgxs.Library.get_xsdata()

        """

        check_type('inverse_velocity', inverse_velocity,
                   openmc.mgxs.InverseVelocity)
        check_value('energy_groups', inverse_velocity.energy_groups,
                    [self.energy_groups])
        check_value('domain_type', inverse_velocity.domain_type,
                    openmc.mgxs.DOMAIN_TYPES)
        self._check_temperature(temperature)

        i = self._temperature_index(temperature)
        self._inverse_velocity[i] = inverse_velocity.get_xs(
            nuclides=nuclide, xs_type=xs_type, subdomains=subdomain)

    def convert_representation(self, target_representation, num_polar=None,
                               num_azimuthal=None):
        """Produce a new XSdata object with the same data, but converted to the
        new representation (isotropic or angle-dependent).

        This method cannot be used to change the number of polar or
        azimuthal bins of an XSdata object that already uses an angular
        representation. Finally, this method simply uses an arithmetic mean to
        convert from an angular to isotropic representation; no flux-weighting
        is applied and therefore reaction rates will not be preserved.

        Parameters
        ----------
        target_representation : {'isotropic', 'angle'}
            Representation of the MGXS (isotropic or angle-dependent flux
            weighting).
        num_polar : int, optional
            Number of equal width angular bins that the polar angular domain is
            subdivided into. This is required when `target_representation` is
            "angle".

        num_azimuthal : int, optional
            Number of equal width angular bins that the azimuthal angular domain
            is subdivided into. This is required when `target_representation` is
            "angle".

        Returns
        -------
        openmc.XSdata

            Multi-group cross section data with the same data as self, but
            represented as specified in `target_representation`.

        """

        check_value('target_representation', target_representation,
                    _REPRESENTATIONS)
        if target_representation == REPRESENTATION_ANGLE:
            check_type('num_polar', num_polar, Integral)
            check_type('num_azimuthal', num_azimuthal, Integral)
            check_greater_than('num_polar', num_polar, 0)
            check_greater_than('num_azimuthal', num_azimuthal, 0)

        xsdata = copy.deepcopy(self)

        # First handle the case where the current and requested
        # representations are the same
        if target_representation == self.representation:
            # Check to make sure the num_polar and num_azimuthal values match
            if target_representation == REPRESENTATION_ANGLE:
                if num_polar != self.num_polar or num_azimuthal != self.num_azimuthal:
                    raise ValueError("Cannot translate between `angle`"
                                     " representations with different angle"
                                     " bin structures")
            # Nothing to do as the same structure was requested
            return xsdata

        xsdata.representation = target_representation
        # We have different actions depending on the representation conversion
        if target_representation == REPRESENTATION_ISOTROPIC:
            # This is not needed for the correct functionality, but these
            # values are changed back to None for clarity
            xsdata._num_polar = None
            xsdata._num_azimuthal = None

        elif target_representation == REPRESENTATION_ANGLE:
            xsdata.num_polar = num_polar
            xsdata.num_azimuthal = num_azimuthal

        # Reset xs_shapes so it is recalculated the next time it is needed
        xsdata._xs_shapes = None

        for i, temp in enumerate(xsdata.temperatures):
            for xs in ['total', 'absorption', 'fission', 'nu_fission',
                       'scatter_matrix', 'multiplicity_matrix',
                       'prompt_nu_fission', 'delayed_nu_fission',
                       'kappa_fission', 'chi', 'chi_prompt', 'chi_delayed',
                       'beta', 'decay_rate', 'inverse_velocity']:
                # Get the original data
                orig_data = getattr(self, '_' + xs)[i]
                if orig_data is not None:

                    if target_representation == 'isotropic':
                        # Since we are going from angle to isotropic, the
                        # current data is just the average over the angle bins
                        new_data = orig_data.mean(axis=(0, 1))

                    elif target_representation == REPRESENTATION_ANGLE:
                        # Since we are going from isotropic to angle, the
                        # current data is just copied for every angle bin
                        new_shape = (num_polar, num_azimuthal) + \
                            orig_data.shape
                        new_data = np.resize(orig_data, new_shape)

                    setter = getattr(xsdata, 'set_' + xs)
                    setter(new_data, temp)

        return xsdata

    def convert_scatter_format(self, target_format, target_order=None):
        """Produce a new MGXSLibrary object with the same data, but converted
        to the new scatter format and order

        Parameters
        ----------
        target_format : {'tabular', 'legendre', 'histogram'}
            Representation of the scattering angle distribution
        target_order : int
            Either the Legendre target_order, number of bins, or number of
            points used to describe the angular distribution associated with
            each group-to-group transfer probability

        Returns
        -------
        openmc.XSdata
            Multi-group cross section data with the same data as in self, but
            represented as specified in `target_format`.

        """

        check_value('target_format', target_format, _SCATTER_TYPES)
        check_type('target_order', target_order, Integral)
        if target_format == SCATTER_LEGENDRE:
            check_greater_than('target_order', target_order, 0, equality=True)
        else:
            check_greater_than('target_order', target_order, 0)

        xsdata = copy.deepcopy(self)
        xsdata.scatter_format = target_format
        xsdata.order = target_order

        # Reset and re-generate XSdata.xs_shapes with the new scattering format
        xsdata._xs_shapes = None

        for i, temp in enumerate(xsdata.temperatures):
            orig_data = self._scatter_matrix[i]
            new_shape = orig_data.shape[:-1] + (xsdata.num_orders,)
            new_data = np.zeros(new_shape)

            if self.scatter_format == SCATTER_LEGENDRE:
                if target_format == SCATTER_LEGENDRE:
                    # Then we are changing orders and only need to change
                    # dimensionality of the mu data and pad/truncate as needed
                    order = min(xsdata.num_orders, self.num_orders)
                    new_data[..., :order] = orig_data[..., :order]

                elif target_format == SCATTER_TABULAR:
                    mu = np.linspace(-1, 1, xsdata.num_orders)
                    # Evaluate the legendre on the mu grid
                    for imu in range(len(mu)):
                        for l in range(self.num_orders):
                            new_data[..., imu] += (
                                 (l + 0.5) * eval_legendre(l, mu[imu]) *
                                 orig_data[..., l])

                elif target_format == SCATTER_HISTOGRAM:
                    # This code uses the vectorized integration capabilities
                    # instead of having an isotropic and angle representation
                    # path.
                    # Set the histogram mu grid
                    mu = np.linspace(-1, 1, xsdata.num_orders + 1)
                    # For every bin perform simpson integration of a finely
                    # sampled orig_data
                    for h_bin in range(xsdata.num_orders):
                        mu_fine = np.linspace(mu[h_bin], mu[h_bin + 1], _NMU)
                        table_fine = np.zeros(new_data.shape[:-1] + (_NMU,))
                        for imu in range(len(mu_fine)):
                            for l in range(self.num_orders):
                                table_fine[..., imu] += ((l + 0.5)
                                     * eval_legendre(l, mu_fine[imu]) *
                                     orig_data[..., l])
                        new_data[..., h_bin] = simps(table_fine, mu_fine)

            elif self.scatter_format == SCATTER_TABULAR:
                # Calculate the mu points of the current data
                mu_self = np.linspace(-1, 1, self.num_orders)

                if target_format == SCATTER_LEGENDRE:
                    # Find the Legendre coefficients via integration. To best
                    # use the vectorized integration capabilities of scipy,
                    # this is done with fixed sample integration routines.
                    mu_fine = np.linspace(-1, 1, _NMU)
                    for l in range(xsdata.num_orders):
                        y = (interp1d(mu_self, orig_data)(mu_fine) *
                             eval_legendre(l, mu_fine))
                        new_data[..., l] = simps(y, mu_fine)

                elif target_format == SCATTER_TABULAR:
                    # Simply use an interpolating function to get the new data
                    mu = np.linspace(-1, 1, xsdata.num_orders)
                    new_data[..., :] = interp1d(mu_self, orig_data)(mu)

                elif target_format == SCATTER_HISTOGRAM:
                    # Use an interpolating function to do the bin-wise
                    # integrals
                    mu = np.linspace(-1, 1, xsdata.num_orders + 1)

                    # Like the tabular -> legendre path above, this code will
                    # be written to utilize the vectorized integration
                    # capabilities instead of having an isotropic and
                    # angle representation path.
                    interp = interp1d(mu_self, orig_data)
                    for h_bin in range(xsdata.num_orders):
                        mu_fine = np.linspace(mu[h_bin], mu[h_bin + 1], _NMU)
                        new_data[..., h_bin] = simps(interp(mu_fine), mu_fine)

            elif self.scatter_format == SCATTER_HISTOGRAM:
                # The histogram format does not have enough information to
                # convert to the other forms without inducing some amount of
                # error. We will make the assumption that the center of the bin
                # has the value of the bin. The mu=-1 and 1 points will be
                # extrapolated from the shape.
                mu_midpoint = np.linspace(-1, 1, self.num_orders,
                                          endpoint=False)
                mu_midpoint += (mu_midpoint[1] - mu_midpoint[0]) * 0.5
                interp = interp1d(mu_midpoint, orig_data,
                                  fill_value='extrapolate')
                # Now get the distribution normalization factor to take from
                # an integral quantity to a point-wise quantity
                norm = float(self.num_orders) / 2.0

                # We now have a tabular distribution in tab_data on mu_self.
                # We now proceed just like the tabular branch above.
                if target_format == SCATTER_LEGENDRE:
                    # find the legendre coefficients via integration. To best
                    # use the vectorized integration capabilities of scipy,
                    # this will be done with fixed sample integration routines.
                    mu_fine = np.linspace(-1, 1, _NMU)
                    for l in range(xsdata.num_orders):
                        y = interp(mu_fine) * norm * eval_legendre(l, mu_fine)
                        new_data[..., l] = simps(y, mu_fine)

                elif target_format == SCATTER_TABULAR:
                    # Simply use an interpolating function to get the new data
                    mu = np.linspace(-1, 1, xsdata.num_orders)
                    new_data[..., :] = interp(mu) * norm

                elif target_format == SCATTER_HISTOGRAM:
                    # Use an interpolating function to do the bin-wise
                    # integrals
                    mu = np.linspace(-1, 1, xsdata.num_orders + 1)

                    # Like the tabular -> legendre path above, this code will
                    # be written to utilize the vectorized integration
                    # capabilities instead of having an isotropic and
                    # angle representation path.
                    for h_bin in range(xsdata.num_orders):
                        mu_fine = np.linspace(mu[h_bin], mu[h_bin + 1], _NMU)
                        new_data[..., h_bin] = \
                            norm * simps(interp(mu_fine), mu_fine)

            # Remove small values resulting from numerical precision issues
            new_data[..., np.abs(new_data) < 1.E-10] = 0.

            xsdata.set_scatter_matrix(new_data, temp)

        return xsdata

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
            if self.representation == REPRESENTATION_ANGLE:
                if self.num_azimuthal is not None:
                    grp.attrs['num_azimuthal'] = self.num_azimuthal

                if self.num_polar is not None:
                    grp.attrs['num_polar'] = self.num_polar

        grp.attrs['scatter_shape'] = np.string_("[G][G'][Order]")
        if self.scatter_format is not None:
            grp.attrs['scatter_format'] = np.string_(self.scatter_format)
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

                if self._chi_prompt[i] is not None:
                    xs_grp.create_dataset("chi-prompt",
                                          data=self._chi_prompt[i])

                if self._chi_delayed[i] is not None:
                    xs_grp.create_dataset("chi-delayed",
                                          data=self._chi_delayed[i])

                if self._nu_fission[i] is None and \
                   (self._delayed_nu_fission[i] is None or \
                    self._prompt_nu_fission[i] is None):
                    raise ValueError('nu-fission or prompt-nu-fission and '
                                     'delayed-nu-fission data must be '
                                     'provided when writing the HDF5 library')

                if self._nu_fission[i] is not None:
                    xs_grp.create_dataset("nu-fission",
                                          data=self._nu_fission[i])

                if self._prompt_nu_fission[i] is not None:
                    xs_grp.create_dataset("prompt-nu-fission",
                                          data=self._prompt_nu_fission[i])

                if self._delayed_nu_fission[i] is not None:
                    xs_grp.create_dataset("delayed-nu-fission",
                                          data=self._delayed_nu_fission[i])

                if self._beta[i] is not None:
                    xs_grp.create_dataset("beta", data=self._beta[i])

                if self._decay_rate[i] is not None:
                    xs_grp.create_dataset("decay rate",
                                          data=self._decay_rate[i])

            if self._scatter_matrix[i] is None:
                raise ValueError('Scatter matrix must be provided when '
                                 'writing the HDF5 library')

            # Get the sparse scattering data to print to the library
            G = self.energy_groups.num_groups
            if self.representation == REPRESENTATION_ISOTROPIC:
                Np = 1
                Na = 1
            elif self.representation == REPRESENTATION_ANGLE:
                Np = self.num_polar
                Na = self.num_azimuthal

            g_out_bounds = np.zeros((Np, Na, G, 2), dtype=np.int)
            for p in range(Np):
                for a in range(Na):
                    for g_in in range(G):
                        if self.scatter_format == SCATTER_LEGENDRE:
                            if self.representation == REPRESENTATION_ISOTROPIC:
                                matrix = \
                                    self._scatter_matrix[i][g_in, :, 0]
                            elif self.representation == REPRESENTATION_ANGLE:
                                matrix = \
                                    self._scatter_matrix[i][p, a, g_in, :, 0]
                        else:
                            if self.representation == REPRESENTATION_ISOTROPIC:
                                matrix = \
                                    np.sum(self._scatter_matrix[i][g_in, :, :],
                                           axis=1)
                            elif self.representation == REPRESENTATION_ANGLE:
                                matrix = \
                                    np.sum(self._scatter_matrix[i][p, a, g_in, :, :],
                                           axis=1)
                        nz = np.nonzero(matrix)
                        # It is possible that there only zeros in matrix
                        # and therefore nz will be empty, in that case set
                        # g_out_bounds to 0s
                        if len(nz[0]) == 0:
                            g_out_bounds[p, a, g_in, :] = 0
                        else:
                            g_out_bounds[p, a, g_in, 0] = nz[0][0]
                            g_out_bounds[p, a, g_in, 1] = nz[0][-1]

            # Now create the flattened scatter matrix array
            flat_scatt = []
            for p in range(Np):
                for a in range(Na):
                    if self.representation == REPRESENTATION_ISOTROPIC:
                        matrix = self._scatter_matrix[i][:, :, :]
                    elif self.representation == REPRESENTATION_ANGLE:
                        matrix = self._scatter_matrix[i][p, a, :, :, :]
                    for g_in in range(G):
                        for g_out in range(g_out_bounds[p, a, g_in, 0],
                                           g_out_bounds[p, a, g_in, 1] + 1):
                            for l in range(len(matrix[g_in, g_out, :])):
                                flat_scatt.append(matrix[g_in, g_out, l])

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
                        if self.representation == REPRESENTATION_ISOTROPIC:
                            matrix = self._multiplicity_matrix[i][:, :]
                        elif self.representation == REPRESENTATION_ANGLE:
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
            if self.representation == REPRESENTATION_ISOTROPIC:
                scatt_grp.create_dataset("g_min", data=g_out_bounds[0, 0, :, 0])
                scatt_grp.create_dataset("g_max", data=g_out_bounds[0, 0, :, 1])
            elif self.representation == REPRESENTATION_ANGLE:
                scatt_grp.create_dataset("g_min", data=g_out_bounds[:, :, :, 0])
                scatt_grp.create_dataset("g_max", data=g_out_bounds[:, :, :, 1])

            # Add the kinetics data
            if self._inverse_velocity[i] is not None:
                xs_grp.create_dataset("inverse-velocity",
                                      data=self._inverse_velocity[i])

    @classmethod
    def from_hdf5(cls, group, name, energy_groups, num_delayed_groups):
        """Generate XSdata object from an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from
        name : str
            Name of the mgxs data set.
        energy_groups : openmc.mgxs.EnergyGroups
            Energy group structure
        num_delayed_groups : int
            Number of delayed groups

        Returns
        -------
        openmc.XSdata
            Multi-group cross section data

        """

        # Get a list of all the subgroups which will contain our temperature
        # strings
        subgroups = group.keys()
        temperatures = []
        for subgroup in subgroups:
            if subgroup != 'kTs':
                temperatures.append(subgroup)

        # To ensure the actual floating point temperature used when creating
        # the new library is consistent with that used when originally creating
        # the file, get the floating point temperatures straight from the kTs
        # group.
        kTs_group = group['kTs']
        float_temperatures = []
        for temperature in temperatures:
            kT = kTs_group[temperature][()]
            float_temperatures.append(kT / openmc.data.K_BOLTZMANN)

        attrs = group.attrs.keys()
        if 'representation' in attrs:
            representation = group.attrs['representation'].decode()
        else:
            representation = REPRESENTATION_ISOTROPIC

        data = cls(name, energy_groups, float_temperatures, representation,
                   num_delayed_groups)

        if 'scatter_format' in attrs:
            data.scatter_format = group.attrs['scatter_format'].decode()

        # Get the remaining optional attributes
        if 'atomic_weight_ratio' in attrs:
            data.atomic_weight_ratio = group.attrs['atomic_weight_ratio']
        if 'order' in attrs:
            data.order = group.attrs['order']
        if data.representation == REPRESENTATION_ANGLE:
            data.num_azimuthal = group.attrs['num_azimuthal']
            data.num_polar = group.attrs['num_polar']

        # Read the temperature-dependent datasets
        for temp, float_temp in zip(temperatures, float_temperatures):
            xs_types = ['total', 'absorption', 'fission', 'kappa-fission',
                        'chi', 'chi-prompt', 'chi-delayed', 'nu-fission',
                        'prompt-nu-fission', 'delayed-nu-fission', 'beta',
                        'decay rate', 'inverse-velocity']

            temperature_group = group[temp]

            for xs_type in xs_types:
                set_func = 'set_' + xs_type.replace(' ', '_').replace('-', '_')
                if xs_type in temperature_group:
                    getattr(data, set_func)(temperature_group[xs_type][()],
                                            float_temp)

            scatt_group = temperature_group['scatter_data']

            # Get scatter matrix and 'un-flatten' it
            g_max = scatt_group['g_max']
            g_min = scatt_group['g_min']
            flat_scatter = scatt_group['scatter_matrix'][()]
            scatter_matrix = np.zeros(data.xs_shapes["[G][G'][Order]"])
            G = data.energy_groups.num_groups
            if data.representation == REPRESENTATION_ISOTROPIC:
                Np = 1
                Na = 1
            elif data.representation == REPRESENTATION_ANGLE:
                Np = data.num_polar
                Na = data.num_azimuthal
            flat_index = 0
            for p in range(Np):
                for a in range(Na):
                    for g_in in range(G):
                        if data.representation == REPRESENTATION_ISOTROPIC:
                            g_mins = g_min[g_in]
                            g_maxs = g_max[g_in]
                        elif data.representation == REPRESENTATION_ANGLE:
                            g_mins = g_min[p, a, g_in]
                            g_maxs = g_max[p, a, g_in]
                        for g_out in range(g_mins - 1, g_maxs):
                            for ang in range(data.num_orders):
                                if data.representation == REPRESENTATION_ISOTROPIC:
                                    scatter_matrix[g_in, g_out, ang] = \
                                        flat_scatter[flat_index]
                                elif data.representation == REPRESENTATION_ANGLE:
                                    scatter_matrix[p, a, g_in, g_out, ang] = \
                                        flat_scatter[flat_index]
                                flat_index += 1
            data.set_scatter_matrix(scatter_matrix, float_temp)

            # Repeat for multiplicity
            if 'multiplicity_matrix' in scatt_group:
                flat_mult = scatt_group['multiplicity_matrix'][()]
                mult_matrix = np.zeros(data.xs_shapes["[G][G']"])
                flat_index = 0
                for p in range(Np):
                    for a in range(Na):
                        for g_in in range(G):
                            if data.representation == REPRESENTATION_ISOTROPIC:
                                g_mins = g_min[g_in]
                                g_maxs = g_max[g_in]
                            elif data.representation == REPRESENTATION_ANGLE:
                                g_mins = g_min[p, a, g_in]
                                g_maxs = g_max[p, a, g_in]
                            for g_out in range(g_mins - 1, g_maxs):
                                if data.representation == REPRESENTATION_ISOTROPIC:
                                    mult_matrix[g_in, g_out] = \
                                        flat_mult[flat_index]
                                elif data.representation == REPRESENTATION_ANGLE:
                                    mult_matrix[p, a, g_in, g_out] = \
                                        flat_mult[flat_index]
                                flat_index += 1
                data.set_multiplicity_matrix(mult_matrix, float_temp)

        return data


class MGXSLibrary:
    """Multi-Group Cross Sections file used for an OpenMC simulation.
    Corresponds directly to the MG version of the cross_sections.xml input
    file.

    Parameters
    ----------
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure
    num_delayed_groups : int
        Num delayed groups

    Attributes
    ----------
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure.
    num_delayed_groups : int
        Num delayed groups
    xsdatas : Iterable of openmc.XSdata
        Iterable of multi-Group cross section data objects
    """

    def __init__(self, energy_groups, num_delayed_groups=0):
        self.energy_groups = energy_groups
        self.num_delayed_groups = num_delayed_groups
        self._xsdatas = []

    def __deepcopy__(self, memo):
        existing = memo.get(id(self))

        # If this is the first time we have tried to copy this object, copy it
        if existing is None:
            clone = type(self).__new__(type(self))
            clone._energy_groups = copy.deepcopy(self.energy_groups, memo)
            clone._num_delayed_groups = self.num_delayed_groups
            clone._xsdatas = copy.deepcopy(self.xsdatas, memo)

            memo[id(self)] = clone

            return clone

        # If this object has been copied before, return the first copy made
        else:
            return existing

    @property
    def energy_groups(self):
        return self._energy_groups

    @property
    def num_delayed_groups(self):
        return self._num_delayed_groups

    @property
    def xsdatas(self):
        return self._xsdatas

    @property
    def names(self):
        return [xsdata.name for xsdata in self.xsdatas]

    @energy_groups.setter
    def energy_groups(self, energy_groups):
        check_type('energy groups', energy_groups, openmc.mgxs.EnergyGroups)
        self._energy_groups = energy_groups

    @num_delayed_groups.setter
    def num_delayed_groups(self, num_delayed_groups):
        check_type('num_delayed_groups', num_delayed_groups, Integral)
        check_greater_than('num_delayed_groups', num_delayed_groups, 0,
                           equality=True)
        check_less_than('num_delayed_groups', num_delayed_groups,
                        openmc.mgxs.MAX_DELAYED_GROUPS, equality=True)
        self._num_delayed_groups = num_delayed_groups

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
        """Add multiple XSdatas to the file.

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

    def get_by_name(self, name):
        """Access the XSdata objects by name

        Parameters
        ----------
        name : str
            Name of openmc.XSdata object to obtain

        Returns
        -------
        result : openmc.XSdata or None
            Provides the matching XSdata object or None, if not found

        """
        check_type("name", name, str)
        result = None
        for xsdata in self.xsdatas:
            if name == xsdata.name:
                result = xsdata
        return result

    def convert_representation(self, target_representation, num_polar=None,
                               num_azimuthal=None):
        """Produce a new XSdata object with the same data, but converted to the
        new representation (isotropic or angle-dependent).

        This method cannot be used to change the number of polar or
        azimuthal bins of an XSdata object that already uses an angular
        representation. Finally, this method simply uses an arithmetic mean to
        convert from an angular to isotropic representation; no flux-weighting
        is applied and therefore the reaction rates will not be preserved.

        Parameters
        ----------
        target_representation : {'isotropic', 'angle'}
            Representation of the MGXS (isotropic or angle-dependent flux
            weighting).
        num_polar : int, optional
            Number of equal width angular bins that the polar angular domain is
            subdivided into. This is required when `target_representation` is
            "angle".
        num_azimuthal : int, optional
            Number of equal width angular bins that the azimuthal angular domain
            is subdivided into. This is required when `target_representation` is
            "angle".

        Returns
        -------
        openmc.MGXSLibrary
            Multi-group Library with the same data as self, but represented as
            specified in `target_representation`.

        """

        library = copy.deepcopy(self)
        for i, xsdata in enumerate(self.xsdatas):
            library.xsdatas[i] = \
                xsdata.convert_representation(target_representation,
                                              num_polar, num_azimuthal)
        return library

    def convert_scatter_format(self, target_format, target_order):
        """Produce a new MGXSLibrary object with the same data, but converted
        to the new scatter format and order

        Parameters
        ----------
        target_format : {'tabular', 'legendre', 'histogram'}
            Representation of the scattering angle distribution
        target_order : int
            Either the Legendre target_order, number of bins, or number of
            points used to describe the angular distribution associated with
            each group-to-group transfer probability

        Returns
        -------
        openmc.MGXSLibrary
            Multi-group Library with the same data as self, but with the scatter
            format represented as specified in `target_format` and
            `target_order`.

        """

        library = copy.deepcopy(self)
        for i, xsdata in enumerate(self.xsdatas):
            library.xsdatas[i] = \
                xsdata.convert_scatter_format(target_format, target_order)

        return library

    def export_to_hdf5(self, filename='mgxs.h5', libver='earliest'):
        """Create an hdf5 file that can be used for a simulation.

        Parameters
        ----------
        filename : str
            Filename of file, default is mgxs.h5.
        libver : {'earliest', 'latest'}
            Compatibility mode for the HDF5 file. 'latest' will produce files
            that are less backwards compatible but have performance benefits.

        """

        check_type('filename', filename, str)

        # Create and write to the HDF5 file
        file = h5py.File(filename, "w", libver=libver)
        file.attrs['filetype'] = np.string_(_FILETYPE_MGXS_LIBRARY)
        file.attrs['version'] = [_VERSION_MGXS_LIBRARY, 0]
        file.attrs['energy_groups'] = self.energy_groups.num_groups
        file.attrs['delayed_groups'] = self.num_delayed_groups
        file.attrs['group structure'] = self.energy_groups.group_edges

        for xsdata in self._xsdatas:
            xsdata.to_hdf5(file)

        file.close()

    @classmethod
    def from_hdf5(cls, filename=None):
        """Generate an MGXS Library from an HDF5 group or file

        Parameters
        ----------
        filename : str, optional
            Name of HDF5 file containing MGXS data. Default is None.
            If not provided, the value of the OPENMC_MG_CROSS_SECTIONS
            environmental variable will be used

        Returns
        -------
        openmc.MGXSLibrary
            Multi-group cross section data object.

        """
        # If filename is None, get the cross sections from the
        # OPENMC_CROSS_SECTIONS environment variable
        if filename is None:
            filename = os.environ.get('OPENMC_MG_CROSS_SECTIONS')

        # Check to make sure there was an environmental variable.
        if filename is None:
            raise ValueError("Either path or OPENMC_MG_CROSS_SECTIONS "
                             "environmental variable must be set")

        check_type('filename', filename, str)
        file = h5py.File(filename, 'r')

        # Check filetype and version
        check_filetype_version(file, _FILETYPE_MGXS_LIBRARY,
                               _VERSION_MGXS_LIBRARY)

        group_structure = file.attrs['group structure']
        num_delayed_groups = file.attrs['delayed_groups']
        energy_groups = openmc.mgxs.EnergyGroups(group_structure)
        data = cls(energy_groups, num_delayed_groups)

        for group_name, group in file.items():
            data.add_xsdata(openmc.XSdata.from_hdf5(group, group_name,
                                                    energy_groups,
                                                    num_delayed_groups))

        return data
