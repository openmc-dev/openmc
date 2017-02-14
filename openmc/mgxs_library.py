import copy
from numbers import Real, Integral
import os

from six import string_types
import numpy as np
import h5py

import openmc
import openmc.mgxs
from openmc.checkvalue import check_type, check_value, check_greater_than, \
    check_iterable_type, check_less_than


# Supported incoming particle MGXS angular treatment representations
_REPRESENTATIONS = ['isotropic', 'angle']
# Supported scattering angular distribution representations
_SCATTER_TYPES = ['tabular', 'legendre', 'histogram']
# List of MGXS dimension types
_XS_SHAPES = ["[G][G'][Order]", "[G]", "[G']", "[G][G']", "[DG]", "[DG][G]",
              "[DG][G']", "[DG][G][G']"]
# Number of mu points for conversion between scattering formats
_NMU = 257


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

    def __init__(self, name, energy_groups, temperatures=[294.],
                 representation='isotropic', num_delayed_groups=0):

        # Initialize class attributes
        self.name = name
        self.energy_groups = energy_groups
        self.num_delayed_groups = num_delayed_groups
        self.temperatures = temperatures
        self.representation = representation
        self._atomic_weight_ratio = None
        self._fissionable = False
        self._scatter_format = 'legendre'
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
        if self._order is not None:
            if self._scatter_format in (None, 'legendre'):
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
            self._xs_shapes["[DG'][G']"] = (self.num_delayed_groups,
                                            self.energy_groups.num_groups)
            self._xs_shapes["[DG][G][G']"] = (self.num_delayed_groups,
                                              self.energy_groups.num_groups,
                                              self.energy_groups.num_groups)

            self._xs_shapes["[G][G'][Order]"] \
                = (self.energy_groups.num_groups,
                   self.energy_groups.num_groups, self.num_orders)

            # If representation is by angle prepend num polar and num azim
            if self.representation == 'angle':
                for key, shapes in self._xs_shapes.items():
                    self._xs_shapes[key] \
                        = (self.num_polar, self.num_azimuthal) + shapes

        return self._xs_shapes

    @name.setter
    def name(self, name):
        check_type('name for XSdata', name, string_types)
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

    @num_delayed_groups.setter
    def num_delayed_groups(self, num_delayed_groups):

        # Check validity of num_delayed_groups
        check_type('num_delayed_groups', num_delayed_groups, Integral)
        check_less_than('num_delayed_groups', num_delayed_groups,
                        openmc.mgxs.MAX_DELAYED_GROUPS, equality=True)
        check_greater_than('num_delayed_groups', num_delayed_groups, 0,
                           equality=True)
        self._num_delayed_groups = num_delayed_groups

    @representation.setter
    def representation(self, representation):

        # Check it is of valid value.
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

        # Get the accepted shapes for this xs
        shapes = [self.xs_shapes["[G]"]]

        # Convert to a numpy array so we can easily get the shape for checking
        total = np.asarray(total)
        check_value('total shape', total.shape, shapes)
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

        # Get the accepted shapes for this xs
        shapes = [self.xs_shapes["[G]"]]

        # Convert to a numpy array so we can easily get the shape for checking
        absorption = np.asarray(absorption)
        check_value('absorption shape', absorption.shape, shapes)
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

        # Get the accepted shapes for this xs
        shapes = [self.xs_shapes["[G]"]]

        # Convert to a numpy array so we can easily get the shape for checking
        fission = np.asarray(fission)
        check_value('fission shape', fission.shape, shapes)
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

        # Get the accepted shapes for this xs
        shapes = [self.xs_shapes["[G]"]]

        # Convert to a numpy array so we can easily get the shape for checking
        kappa_fission = np.asarray(kappa_fission)
        check_value('kappa fission shape', kappa_fission.shape, shapes)
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

        # Get the accepted shapes for this xs
        shapes = [self.xs_shapes["[G']"]]

        # Convert to a numpy array so we can easily get the shape for checking
        chi = np.asarray(chi)
        check_value('chi shape', chi.shape, shapes)
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        self._chi[i] = chi

    def set_chi_prompt(self, chi_prompt, temperature=294.):
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
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        self._chi_prompt[i] = chi_prompt

    def set_chi_delayed(self, chi_delayed, temperature=294.):
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
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        self._chi_delayed[i] = chi_delayed

    def set_beta(self, beta, temperature=294.):
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
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        self._beta[i] = beta

    def set_decay_rate(self, decay_rate, temperature=294.):
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
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        self._decay_rate[i] = decay_rate

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

        # Get the accepted shapes for this xs
        shapes = [self.xs_shapes["[G][G'][Order]"]]

        # Convert to a numpy array so we can easily get the shape for checking
        scatter = np.asarray(scatter)
        check_iterable_type('scatter', scatter, Real,
                            max_depth=len(scatter.shape))
        check_value('scatter shape', scatter.shape, shapes)
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

        # Get the accepted shapes for this xs
        shapes = [self.xs_shapes["[G][G']"]]

        # Convert to a numpy array so we can easily get the shape for checking
        multiplicity = np.asarray(multiplicity)
        check_iterable_type('multiplicity', multiplicity, Real,
                            max_depth=len(multiplicity.shape))
        check_value('multiplicity shape', multiplicity.shape, shapes)
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

        # Get the accepted shapes for this xs
        shapes = [self.xs_shapes["[G]"], self.xs_shapes["[G][G']"]]

        # Convert to a numpy array so we can easily get the shape for checking
        nu_fission = np.asarray(nu_fission)
        check_value('nu_fission shape', nu_fission.shape, shapes)
        check_iterable_type('nu_fission', nu_fission, Real,
                            max_depth=len(nu_fission.shape))
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        self._nu_fission[i] = nu_fission
        if np.sum(nu_fission) > 0.0:
            self._fissionable = True

    def set_prompt_nu_fission(self, prompt_nu_fission, temperature=294.):
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
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        self._prompt_nu_fission[i] = prompt_nu_fission
        if np.sum(prompt_nu_fission) > 0.0:
            self._fissionable = True

    def set_delayed_nu_fission(self, delayed_nu_fission, temperature=294.):
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
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        self._delayed_nu_fission[i] = delayed_nu_fission
        if np.sum(delayed_nu_fission) > 0.0:
            self._fissionable = True

    def set_inverse_velocity(self, inv_vel, temperature=294.):
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
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        self._inverse_velocity[i] = inv_vel

    def set_total_mgxs(self, total, temperature=294., nuclide='total',
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
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        self._total[i] = total.get_xs(nuclides=nuclide, xs_type=xs_type,
                                      subdomains=subdomain)

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
        openmc.mgxs.Library.get_xsdata()

        """

        check_type('absorption', absorption, openmc.mgxs.AbsorptionXS)
        check_value('energy_groups', absorption.energy_groups,
                    [self.energy_groups])
        check_value('domain_type', absorption.domain_type,
                    openmc.mgxs.DOMAIN_TYPES)
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        self._absorption[i] = absorption.get_xs(nuclides=nuclide,
                                                xs_type=xs_type,
                                                subdomains=subdomain)

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
        openmc.mgxs.Library.get_xsdata()

        """

        check_type('fission', fission, openmc.mgxs.FissionXS)
        check_value('energy_groups', fission.energy_groups,
                    [self.energy_groups])
        check_value('domain_type', fission.domain_type,
                    openmc.mgxs.DOMAIN_TYPES)
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        self._fission[i] = fission.get_xs(nuclides=nuclide,
                                          xs_type=xs_type,
                                          subdomains=subdomain)

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
        openmc.mgxs.Library.get_xsdata()

        """

        check_type('nu_fission', nu_fission, (openmc.mgxs.NuFissionXS,
                                              openmc.mgxs.NuFissionMatrixXS))
        check_value('energy_groups', nu_fission.energy_groups,
                    [self.energy_groups])
        check_value('domain_type', nu_fission.domain_type,
                    openmc.mgxs.DOMAIN_TYPES)
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        self._nu_fission[i] = nu_fission.get_xs(nuclides=nuclide,
                                                xs_type=xs_type,
                                                subdomains=subdomain)

        if np.sum(self._nu_fission) > 0.0:
            self._fissionable = True

    def set_prompt_nu_fission_mgxs(self, prompt_nu_fission, temperature=294.,
                                   nuclide='total', xs_type='macro',
                                   subdomain=None):
        """Sets the prompt-nu-fission cross section.

        This method allows for an openmc.mgxs.PromptNuFissionXS or
        openmc.mgxs.PromptNuFissionMatrixXS to be used to set the
        prompt-nu-fission cross section for this XSdata object.

        Parameters
        ----------
        prompt_nu_fission: openmc.mgxs.PromptNuFissionXS or openmc.mgxs.PromptNuFissionMatrixXS
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
                   (openmc.mgxs.PromptNuFissionXS,
                    openmc.mgxs.PromptNuFissionMatrixXS))
        check_value('energy_groups', prompt_nu_fission.energy_groups,
                    [self.energy_groups])
        check_value('domain_type', prompt_nu_fission.domain_type,
                    openmc.mgxs.DOMAIN_TYPES)
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        self._prompt_nu_fission[i] = prompt_nu_fission.get_xs(
            nuclides=nuclide, xs_type=xs_type, subdomains=subdomain)

        if np.sum(self._prompt_nu_fission) > 0.0:
            self._fissionable = True

    def set_delayed_nu_fission_mgxs(self, delayed_nu_fission, temperature=294.,
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
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        self._delayed_nu_fission[i] = delayed_nu_fission.get_xs(
            nuclides=nuclide, xs_type=xs_type, subdomains=subdomain)

        if np.sum(self._delayed_nu_fission) > 0.0:
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
        openmc.mgxs.Library.get_xsdata()

        """

        check_type('kappa_fission', k_fission, openmc.mgxs.KappaFissionXS)
        check_value('energy_groups', k_fission.energy_groups,
                    [self.energy_groups])
        check_value('domain_type', k_fission.domain_type,
                    openmc.mgxs.DOMAIN_TYPES)
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        self._kappa_fission[i] = k_fission.get_xs(nuclides=nuclide,
                                                  xs_type=xs_type,
                                                  subdomains=subdomain)

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
        openmc.mgxs.Library.get_xsdata()

        """

        check_type('chi', chi, openmc.mgxs.Chi)
        check_value('energy_groups', chi.energy_groups, [self.energy_groups])
        check_value('domain_type', chi.domain_type, openmc.mgxs.DOMAIN_TYPES)
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        self._chi[i] = chi.get_xs(nuclides=nuclide, xs_type=xs_type,
                                  subdomains=subdomain)

    def set_chi_prompt_mgxs(self, chi_prompt, temperature=294.,
                            nuclide='total', xs_type='macro', subdomain=None):
        """This method allows for an openmc.mgxs.ChiPrompt
        to be used to set chi-prompt for this XSdata object.

        Parameters
        ----------
        chi_prompt: openmc.mgxs.ChiPrompt
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

        check_type('chi_prompt', chi_prompt, openmc.mgxs.ChiPrompt)
        check_value('energy_groups', chi_prompt.energy_groups,
                    [self.energy_groups])
        check_value('domain_type', chi_prompt.domain_type,
                    openmc.mgxs.DOMAIN_TYPES)
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        self._chi_prompt[i] = chi_prompt.get_xs(nuclides=nuclide,
                                                xs_type=xs_type,
                                                subdomains=subdomain)

    def set_chi_delayed_mgxs(self, chi_delayed, temperature=294.,
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
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        self._chi_delayed[i] = chi_delayed.get_xs(nuclides=nuclide,
                                                  xs_type=xs_type,
                                                  subdomains=subdomain)

    def set_beta_mgxs(self, beta, temperature=294.,
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
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        self._beta[i] = beta.get_xs(nuclides=nuclide,
                                    xs_type=xs_type,
                                    subdomains=subdomain)

    def set_decay_rate_mgxs(self, decay_rate, temperature=294.,
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
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        self._decay_rate[i] = decay_rate.get_xs(nuclides=nuclide,
                                                xs_type=xs_type,
                                                subdomains=subdomain)

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
        openmc.mgxs.Library.get_xsdata()

        """

        check_type('scatter', scatter, openmc.mgxs.ScatterMatrixXS)
        check_value('energy_groups', scatter.energy_groups,
                    [self.energy_groups])
        check_value('domain_type', scatter.domain_type,
                    openmc.mgxs.DOMAIN_TYPES)
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        # Set the value of scatter_format based on the same value within
        # scatter
        self.scatter_format = scatter.scatter_format

        # If the user has not defined XSdata.order, then we will set
        # the order based on the data within scatter.
        # Otherwise, we will check to see that XSdata.order matches
        # the order of scatter
        if self.scatter_format == 'legendre':
            if self.order is None:
                self.order = scatter.legendre_order
            else:
                check_value('legendre_order', scatter.legendre_order,
                            [self.order])
        elif self.scatter_format == 'histogram':
            if self.order is None:
                self.order = scatter.histogram_bins
            else:
                check_value('histogram_bins', scatter.histogram_bins,
                            [self.order])

        i = np.where(self.temperatures == temperature)[0][0]
        if self.scatter_format == 'legendre':
            self._scatter_matrix[i] = \
                np.zeros(self.xs_shapes["[G][G'][Order]"])
            # Get the scattering orders in the outermost dimension
            if self.representation == 'isotropic':
                for moment in range(self.num_orders):
                    self._scatter_matrix[i][:, :, moment] = \
                        scatter.get_xs(nuclides=nuclide, xs_type=xs_type,
                                       moment=moment, subdomains=subdomain)
            elif self.representation == 'angle':
                for moment in range(self.num_orders):
                    self._scatter_matrix[i][:, :, :, :, moment] = \
                        scatter.get_xs(nuclides=nuclide, xs_type=xs_type,
                                       moment=moment, subdomains=subdomain)
        else:
            self._scatter_matrix[i] = \
                scatter.get_xs(nuclides=nuclide, xs_type=xs_type,
                               subdomains=subdomain)

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
        openmc.mgxs.Library.get_xsdata()

        """

        check_type('nuscatter', nuscatter, (openmc.mgxs.NuScatterMatrixXS,
                                            openmc.mgxs.MultiplicityMatrixXS))
        check_value('energy_groups', nuscatter.energy_groups,
                    [self.energy_groups])
        check_value('domain_type', nuscatter.domain_type,
                    openmc.mgxs.DOMAIN_TYPES)
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
                        openmc.mgxs.DOMAIN_TYPES)
        i = np.where(self.temperatures == temperature)[0][0]
        nuscatt = nuscatter.get_xs(nuclides=nuclide,
                                   xs_type=xs_type, moment=0,
                                   subdomains=subdomain)
        if isinstance(nuscatter, openmc.mgxs.MultiplicityMatrixXS):
            self._multiplicity_matrix[i] = nuscatt
        else:
            scatt = scatter.get_xs(nuclides=nuclide,
                                   xs_type=xs_type, moment=0,
                                   subdomains=subdomain)
            if scatter.scatter_format == 'histogram':
                scatt = np.sum(scatt, axis=0)
            if nuscatter.scatter_format == 'histogram':
                nuscatt = np.sum(nuscatt, axis=0)
            self._multiplicity_matrix[i] = np.divide(nuscatt, scatt)

        self._multiplicity_matrix[i] = \
            np.nan_to_num(self._multiplicity_matrix[i])

    def set_inverse_velocity_mgxs(self, inverse_velocity, temperature=294.,
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
        check_type('temperature', temperature, Real)
        check_value('temperature', temperature, self.temperatures)

        i = np.where(self.temperatures == temperature)[0][0]
        self._inverse_velocity[i] = inverse_velocity.get_xs(
            nuclides=nuclide, xs_type=xs_type, subdomains=subdomain)

    def convert_representation(self, target_representation, num_polar=None,
                               num_azimuthal=None):
        """Produce a new XSdata object with the same data, but converted to the
        new representation

        Parameters
        ----------
        target_representation : {'isotropic', 'angle'}
            Representation of the MGXS (isotropic or angle-dependent flux
            weighting).
        num_polar : int, optional
            Number of equal width angular bins that the polar angular
            domain is subdivided into. This is required when
            :param:`target_representation` is "angle".
        num_azimuthal : int, optional
            Number of equal width angular bins that the azimuthal angular
            domain is subdivided into. This is required when
            :param:`target_representation` is "angle".

        Returns
        -------
        openmc.XSdata
            Multi-group cross section data with the same data as self, but
            represented as specified in :param:`target_representation`.

        """

        check_value('target_representation', target_representation,
                    _REPRESENTATIONS)
        if target_representation == 'angle':
            check_type('num_polar', num_polar, Integral)
            check_type('num_azimuthal', num_azimuthal, Integral)
            check_greater_than('num_polar', num_polar, 0)
            check_greater_than('num_azimuthal', num_azimuthal, 0)

        xsdata = copy.deepcopy(self)
        if target_representation == self.representation:
            # Check to make sure the num_polar and num_azimuthal values match
            if target_representation == 'angle':
                if num_polar != self.num_polar or num_azimuthal != self.num_azimuthal:
                    raise NotImplementedError("XSdata.convert_representation "
                                              "cannot translate between "
                                              "`angle` representations with "
                                              "different angle bin structures")
            # Nothing to do
            return xsdata

        types = ['total', 'absorption', 'fission', 'nu_fission',
                 'scatter_matrix', 'multiplicity_matrix', 'prompt_nu_fission',
                 'delayed_nu_fission', 'kappa_fission', 'chi', 'chi_prompt',
                 'chi_delayed', 'beta', 'decay_rate', 'inverse_velocity']

        xsdata.representation = target_representation
        # We have different actions depending on the representation conversion
        if target_representation == 'isotropic':
            # This is not needed for the correct functionality, but these
            # values are changed back to None for clarity
            xsdata._num_polar = None
            xsdata._num_azimuthal = None
        elif target_representation == 'angle':
            xsdata.num_polar = num_polar
            xsdata.num_azimuthal = num_azimuthal

        # Reset XSdata.xs_shapes to accomodate the new shape
        xsdata._xs_shapes = None

        for i, temp in enumerate(xsdata.temperatures):
            for xs in types:
                # Get the original data
                orig_data = getattr(self, '_' + xs)[i]
                if orig_data is not None:
                    if target_representation == 'isotropic':
                        # Since we are going from angle to isotropic, the
                        # current data is just averaged over the angle bins
                        new_data = orig_data.mean(axis=(0, 1))
                    elif target_representation == 'angle':
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
            Multi-group cross section data with the same data as self, but
            represented as specified in :param:`target_format`.

        """

        from scipy.interpolate import interp1d
        from scipy.integrate import simps
        from scipy.special import eval_legendre

        check_value('target_format', target_format, _SCATTER_TYPES)
        check_type('target_order', target_order, Integral)
        if target_format == 'legendre':
            check_greater_than('target_order', target_order, 0, equality=True)
        else:
            check_greater_than('target_order', target_order, 0)

        xsdata = copy.deepcopy(self)
        xsdata.scatter_format = target_format
        xsdata.order = target_order
        # Reset XSdata.xs_shapes to accomodate the new shape
        xsdata._xs_shapes = None
        xsdata.xs_shapes

        # We have to accomodate the following possibilities:
        # histogram -> tabular w/ same or diff order
        # histogram -> legendre
        # histogram -> histogram w/ same or diff order

        for i, temp in enumerate(xsdata.temperatures):
            orig_data = self._scatter_matrix[i]
            new_shape = orig_data.shape[:-1] + (xsdata.num_orders,)
            new_data = np.zeros(new_shape)
            if self.scatter_format == 'legendre':
                if target_format == 'legendre':
                    # Then we are changing orders and only need to change
                    # dimensionality of the mu data and pad/truncate as needed
                    order = min(xsdata.num_orders, self.num_orders)
                    new_data[..., :order] = orig_data[..., :order]
                elif target_format == 'tabular':
                    mu = np.linspace(-1, 1, xsdata.num_orders)
                    # Evaluate the legendre on the tabular grid within mu
                    for imu in range(len(mu)):
                        new_data[..., imu] = \
                            np.sum((l + 0.5) * eval_legendre(l, mu[imu]) *
                                   orig_data[..., l]
                                   for l in range(self.num_orders))
                elif target_format == 'histogram':
                    mu = np.linspace(-1, 1, xsdata.num_orders + 1)
                    # This code will be written to utilize the vectorized
                    # integration capabilities instead of having an isotropic
                    # and angle representation path.
                    for h_bin in range(xsdata.num_orders):
                        mu_fine = np.linspace(mu[h_bin], mu[h_bin + 1], _NMU)
                        table_shape = new_data.shape[:-1] + (_NMU,)
                        table_fine = np.zeros(table_shape)
                        for imu in range(len(mu_fine)):
                            table_fine[..., imu] = \
                                np.sum((l + 0.5) *
                                       eval_legendre(l, mu_fine[imu]) *
                                       orig_data[..., l]
                                       for l in range(self.num_orders))
                        new_data[..., h_bin] = simps(table_fine, mu_fine)

                    # Remove the very small results from numerical precision
                    # issues (allowing conversions to be reproduced exactly)
                    new_data[..., np.abs(new_data) < 1.E-10] = 0.

            elif self.scatter_format == 'tabular':
                mu_self = np.linspace(-1, 1, self.num_orders)
                if target_format == 'legendre':
                    # find the legendre coefficients via integration. To best
                    # use the vectorized integration capabilities of scipy,
                    # this will be done with fixed sample integration routines.
                    mu_fine = np.linspace(-1, 1, _NMU)
                    y = [interp1d(mu_self, orig_data)(mu_fine) *
                         eval_legendre(l, mu_fine)
                         for l in range(xsdata.num_orders)]
                    for l in range(xsdata.num_orders):
                        new_data[..., l] = simps(y[l], mu_fine)

                    # Remove the very small results from numerical precision
                    # issues (allowing conversions to be reproduced exactly)
                    new_data[..., np.abs(new_data) < 1.E-10] = 0.
                elif target_format == 'tabular':
                    # Simply use an interpolating function to get the new data
                    mu = np.linspace(-1, 1, xsdata.num_orders)
                    new_data[..., :] = interp1d(mu_self, orig_data)(mu)
                elif target_format == 'histogram':
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

                    # Remove the very small results from numerical precision
                    # issues (allowing conversions to be reproduced exactly)
                    new_data[..., np.abs(new_data) < 1.E-10] = 0.
            elif self.scatter_format == 'histogram':
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
                if target_format == 'legendre':
                    # find the legendre coefficients via integration. To best
                    # use the vectorized integration capabilities of scipy,
                    # this will be done with fixed sample integration routines.
                    mu_fine = np.linspace(-1, 1, _NMU)
                    y = [interp(mu_fine) * norm * eval_legendre(l, mu_fine)
                         for l in range(xsdata.num_orders)]
                    for l in range(xsdata.num_orders):
                        new_data[..., l] = simps(y[l], mu_fine)

                    # Remove the very small results from numerical precision
                    # issues (allowing conversions to be reproduced exactly)
                    new_data[..., np.abs(new_data) < 1.E-10] = 0.
                elif target_format == 'tabular':
                    # Simply use an interpolating function to get the new data
                    mu = np.linspace(-1, 1, xsdata.num_orders)
                    new_data[..., :] = interp(mu) * norm
                elif target_format == 'histogram':
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

                    # Remove the very small results from numerical precision
                    # issues (allowing conversions to be reproduced exactly)
                    new_data[..., np.abs(new_data) < 1.E-10] = 0.
            xsdata.set_scatter_matrix(new_data, temp)

        return xsdata

    def convert_to_continuous_energy(self, temperature=294.):
        """Converts the XSdata object to an equivalent
        openmc.data.IncidentNeutron object.

        Parameters
        ----------
        temperature : float
            Temperature of dataset to print; defaults to 294K

        Returns
        -------
        openmc.data.IncidentNeutron
            The continuous-energy IncidentNeutron data library equivalent
            to the MGXS data in self
        """

        # Check if this can be performed successfully
        if self.representation == 'angle':
            raise ValueError("Cannot convert angle-dependent MGXS; convert to "
                             "an isotropic representation first")
        required_types = ['absorption', 'scatter_matrix']
        for type_check in required_types:
            if getattr(self, '_' + type_check)[0] is None:
                raise ValueError(type_check + ' data is required')
        if temperature not in self.temperatures:
            raise ValueError("Invalid temperature")

        def convert_xs(group_edges, values):
            cexs = openmc.data.Tabulated1D(group_edges,
                                           np.append(values[::-1],
                                                     [values[0]]),
                                           breakpoints=[len(group_edges)],
                                           interpolation=[1])
            return cexs

        # Build required metadata
        kTs = self.temperatures[:] * openmc.data.K_BOLTZMANN
        if self.atomic_weight_ratio:
            awr = self.atomic_weight_ratio
        else:
            awr = 1.

        # Get the temperature index
        iT = np.where(self.temperatures == temperature)[0][0]
        data = openmc.data.IncidentNeutron(self.name, 1, 1, 0, awr, kTs)
        strT = "{}K".format(int(round(self.temperatures[iT])))
        data.energy = {strT: self.energy_groups.group_edges}
        energy_midpoints = (self.energy_groups.group_edges[1:] +
                            self.energy_groups.group_edges[:-1]) / 2.

        # Elastic MGXS: must explicitly be 0 barns to avoid incorrect
        # elastic_scatter calculation with data available to us in MG library.
        el_rxn = openmc.data.Reaction(2)
        el_rxn.xs[strT] = \
            convert_xs(self.energy_groups.group_edges,
                       np.zeros(self.energy_groups.num_groups))
        data.reactions[2] = el_rxn

        # Absorption MGXS
        abs_rxn = openmc.data.Reaction(102)
        abs_rxn.xs[strT] = convert_xs(self.energy_groups.group_edges,
                                      np.subtract(self._absorption[iT],
                                                  self._fission[iT]))
        data.reactions[102] = abs_rxn

        if self._nu_fission[iT] is not None and self._fission[iT] is not None:
            fiss_rxn = openmc.data.Reaction(18)
            fiss_rxn.xs[strT] = convert_xs(self.energy_groups.group_edges,
                                           self._fission[iT])

            # Get nu_fission and chi from the presence of both, or the
            # nu_fission matrix
            if self._nu_fission[iT].shape == self.xs_shapes["[G][G']"]:
                nu_fiss = np.sum(self._nu_fission[iT], axis=1)[::-1]
                chi = self._nu_fission[iT][::-1] / nu_fiss
            else:
                nu_fiss = self._nu_fission[iT][::-1]
                chi = np.reshape(np.tile(self._chi[iT][::-1],
                                         self.energy_groups.num_groups),
                                 (self.energy_groups.num_groups,
                                  self.energy_groups.num_groups))
            nu = np.divide(nu_fiss, self._fission[iT][::-1])

            # Build a histogram distribution to represent the yield for the
            # incoming groups
            prod = openmc.data.Product()
            prod.yield_ = \
                openmc.data.Tabulated1D(self.energy_groups.group_edges[:-1],
                                        nu, breakpoints=[len(nu)],
                                        interpolation=[1])

            # Now build the outgoing energy distribution using discrete energy
            # lines within each of the outgoing groups
            chi_eouts = []
            for g in range(self.energy_groups.num_groups):
                chi_eouts.append(openmc.stats.Discrete(energy_midpoints,
                                                       chi[g, :]))
                # Ensure the distribution CDF starts with 0
                chi_eouts[-1].c = np.cumsum(chi[g, :]) - chi[g, 0]

            # Create the continuous distribution for the fission
            # The angular distribution will remain None (thus isotropic)
            chi_distrib = openmc.data.ContinuousTabular(
                [len(chi)], [1], self.energy_groups.group_edges[:-1],
                chi_eouts)
            prod.distribution = \
                [openmc.data.UncorrelatedAngleEnergy(energy=chi_distrib)]
            fiss_rxn.products = [prod]
            fiss_rxn.center_of_mass = False
            data.reactions[18] = fiss_rxn

        # Scattering Data
        # First convert the data to a tabular representation
        if self.scatter_format == 'tabular':
            tabular = self
        else:
            tabular = self.convert_scatter_format('tabular', 33)

        # Calculate the isotropic scattering matrix and use that to find the
        # total scattering x/s and outgoing energy distributions
        isotropic_matrix = np.mean(tabular._scatter_matrix[iT], axis=-1)[::-1,
                                                                         ::-1]
        scatt_xs = np.sum(isotropic_matrix, axis=1)
        energy = np.zeros((self.energy_groups.num_groups,
                           self.energy_groups.num_groups))
        for gin in range(self.energy_groups.num_groups):
            energy[gin, :] = isotropic_matrix[gin, :] / scatt_xs[gin]

        # Get the anisotropic but normalized angular distribution
        distrib = np.zeros((self.energy_groups.num_groups,
                            self.energy_groups.num_groups,
                            tabular.num_orders))
        for gin in range(self.energy_groups.num_groups):
            for gout in range(self.energy_groups.num_groups):
                distrib[gin, gout, :] = \
                    np.divide(tabular._scatter_matrix[iT][gin, gout, :],
                              isotropic_matrix[gin, gout])
        distrib = np.nan_to_num(distrib)

        # Incorporate the scattering multiplication, if required
        scatt_prod = openmc.data.Product()
        if self._multiplicity_matrix[iT]:
            yield_ = self._multiplicity_matrix[iT][::-1, ::-1]
            scatt_prod.yield_ = \
                openmc.data.Tabulated1D(self.energy_groups.group_edges[:-1],
                                        yield_, breakpoints=[len(yield_)],
                                        interpolation=[1])

        # Now build the outgoing energy distribution using discrete energy
        # lines within each of the outgoing groups
        scatt_eouts = []
        for g in range(self.energy_groups.num_groups):
            scatt_eouts.append(
                openmc.stats.Tabular(self.energy_groups.group_edges[:-1],
                                     energy[g, :],
                                     interpolation='histogram'))
            # Ensure the distribution CDF starts with 0
            scatt_eouts[-1].c = np.cumsum(energy[g, :]) - energy[g, 0]
        scatt_eouts.append(scatt_eouts[-1])

        # Build the angular distributions associated with each outgoing energy
        # group
        mu = np.linspace(-1., 1., tabular.num_orders)
        scatt_angles = []
        for gin in range(self.energy_groups.num_groups):
            scatt_angles.append([])
            for gout in range(self.energy_groups.num_groups):
                if energy[gin, gout] > 0.:
                    scatt_angles[gin].append(
                        openmc.stats.Tabular(mu, distrib[gin, gout, :]))
                    # Ensure the distribution CDF starts with 0
                    scatt_angles[gin][-1].c = \
                        np.cumsum(distrib[gin, gout, :]) - \
                        distrib[gin, gout, 0]
        scatt_angles.append(scatt_angles[-1])

        # Combine the energy and angle distributions in to a correlated
        # angle/energy object
        scatt_prod.distribution = \
            [openmc.data.CorrelatedAngleEnergy(
                breakpoints=[len(scatt_eouts)], interpolation=[1],
                energy=self.energy_groups.group_edges[:-1],
                energy_out=scatt_eouts, mu=scatt_angles)]

        # Finally build the reaction with the just-calculated information
        # This will be set to the (n,2n) reaction, though any scattering
        # reaction would suffice
        scatt_rxn = openmc.data.Reaction(16)
        scatt_rxn.xs[strT] = \
            convert_xs(self.energy_groups.group_edges, scatt_xs[::-1])
        scatt_rxn.products = [scatt_prod]
        scatt_rxn.center_of_mass = False
        data.reactions[16] = scatt_rxn

        return data

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
            if self.representation == 'isotropic':
                Np = 1
                Na = 1
            elif self.representation == 'angle':
                Np = self.num_polar
                Na = self.num_azimuthal

            g_out_bounds = np.zeros((Np, Na, G, 2), dtype=np.int)
            for p in range(Np):
                for a in range(Na):
                    for g_in in range(G):
                        if self.scatter_format == 'legendre':
                            if self.representation == 'isotropic':
                                matrix = \
                                    self._scatter_matrix[i][g_in, :, 0]
                            elif self.representation == 'angle':
                                matrix = \
                                    self._scatter_matrix[i][p, a, g_in, :, 0]
                        else:
                            if self.representation == 'isotropic':
                                matrix = \
                                    np.sum(self._scatter_matrix[i][g_in, :, :],
                                           axis=1)
                            elif self.representation == 'angle':
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
                    if self.representation == 'isotropic':
                        matrix = self._scatter_matrix[i][:, :, :]
                    elif self.representation == 'angle':
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
                        if self.representation == 'isotropic':
                            matrix = self._multiplicity_matrix[i][:, :]
                        elif self.representation == 'angle':
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
            if self.representation == 'isotropic':
                scatt_grp.create_dataset("g_min", data=g_out_bounds[0, 0, :, 0])
                scatt_grp.create_dataset("g_max", data=g_out_bounds[0, 0, :, 1])
            elif self.representation == 'angle':
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
            kT = kTs_group[temperature].value
            float_temperatures.append(kT / openmc.data.K_BOLTZMANN)

        attrs = group.attrs.keys()
        if 'representation' in attrs:
            representation = group.attrs['representation'].decode()
        else:
            representation = 'isotropic'

        data = cls(name, energy_groups, float_temperatures, representation,
                   num_delayed_groups)

        if 'scatter_format' in attrs:
            data.scatter_format = group.attrs['scatter_format'].decode()

        # Get the remaining optional attributes
        if 'atomic_weight_ratio' in attrs:
            data.atomic_weight_ratio = group.attrs['atomic_weight_ratio']
        if 'order' in attrs:
            data.order = group.attrs['order']
        if data.representation == 'angle':
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
                    getattr(data, set_func)(temperature_group[xs_type].value,
                                            float_temp)

            scatt_group = temperature_group['scatter_data']

            # Get scatter matrix and 'un-flatten' it
            g_max = scatt_group['g_max']
            g_min = scatt_group['g_min']
            flat_scatter = scatt_group['scatter_matrix'].value
            scatter_matrix = np.zeros(data.xs_shapes["[G][G'][Order]"])
            G = data.energy_groups.num_groups
            if data.representation == 'isotropic':
                Np = 1
                Na = 1
            elif data.representation == 'angle':
                Np = data.num_polar
                Na = data.num_azimuthal
            flat_index = 0
            for p in range(Np):
                for a in range(Na):
                    for g_in in range(G):
                        if data.representation == 'isotropic':
                            g_mins = g_min[g_in]
                            g_maxs = g_max[g_in]
                        elif data.representation == 'angle':
                            g_mins = g_min[p, a, g_in]
                            g_maxs = g_max[p, a, g_in]
                        for g_out in range(g_mins - 1, g_maxs):
                            for ang in range(data.num_orders):
                                if data.representation == 'isotropic':
                                    scatter_matrix[g_in, g_out, ang] = \
                                        flat_scatter[flat_index]
                                elif data.representation == 'angle':
                                    scatter_matrix[p, a, g_in, g_out, ang] = \
                                        flat_scatter[flat_index]
                                flat_index += 1
            data.set_scatter_matrix(scatter_matrix, float_temp)

            # Repeat for multiplicity
            if 'multiplicity_matrix' in scatt_group:
                flat_mult = scatt_group['multiplicity_matrix'].value
                mult_matrix = np.zeros(data.xs_shapes["[G][G']"])
                flat_index = 0
                for p in range(Np):
                    for a in range(Na):
                        for g_in in range(G):
                            if data.representation == 'isotropic':
                                g_mins = g_min[g_in]
                                g_maxs = g_max[g_in]
                            elif data.representation == 'angle':
                                g_mins = g_min[p, a, g_in]
                                g_maxs = g_max[p, a, g_in]
                            for g_out in range(g_mins - 1, g_maxs):
                                if data.representation == 'isotropic':
                                    mult_matrix[g_in, g_out] = \
                                        flat_mult[flat_index]
                                elif data.representation == 'angle':
                                    mult_matrix[p, a, g_in, g_out] = \
                                        flat_mult[flat_index]
                                flat_index += 1
                data.set_multiplicity_matrix(mult_matrix, float_temp)

        return data


class MGXSLibrary(object):
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
        """Produce a new MGXSLibrary object with the same data, but converted
        to the new representation

        Parameters
        ----------
        target_representation : {'isotropic', 'angle'}
            Representation of the MGXS (isotropic or angle-dependent flux
            weighting).
        num_polar : int, optional
            Number of equal width angular bins that the polar angular
            domain is subdivided into. This is required when
            :param:`target_representation` is "angle".
        num_azimuthal : int, optional
            Number of equal width angular bins that the azimuthal angular
            domain is subdivided into. This is required when
            :param:`target_representation` is "angle".

        Returns
        -------
        openmc.MGXSLibrary
            Multi-group Library with the same data as self, but represented as
            specified in :param:`target_representation`.

        """

        check_value('target_representation', target_representation,
                    _REPRESENTATIONS)
        if target_representation == 'angle':
            check_type('num_polar', num_polar, Integral)
            check_type('num_azimuthal', num_azimuthal, Integral)
            check_greater_than('num_polar', num_polar, 0)
            check_greater_than('num_azimuthal', num_azimuthal, 0)

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
            Multi-group Library with the same data as self, but with the
            scatter format represented as specified in :param:`target_format`
            and :param:`target_order`.

        """

        check_value('target_format', target_format, _SCATTER_TYPES)
        check_type('target_order', target_order, Integral)
        if target_format == 'legendre':
            check_greater_than('target_order', target_order, 0, equality=True)
        else:
            check_greater_than('target_order', target_order, 0)

        library = copy.deepcopy(self)
        for i, xsdata in enumerate(self.xsdatas):
            library.xsdatas[i] = \
                xsdata.convert_scatter_format(target_format, target_order)

        return library

    def convert_to_continuous_energy(self, h5_filename='ce_mgxs.h5',
                                     library_filename='ce_mgxs.xml'):
        """Converts the MGXSLibrary object to an equivalent
        library of openmc.data.IncidentNeutron objects

        Parameters
        ----------
        h5_filename : str
            HDF5 file to write with all the files
        library_filename : str
            cross_sections.xml file describing the HDF5 file

        """

        library = openmc.data.DataLibrary()
        data = []
        for i, xsdata in enumerate(self.xsdatas):
            data.append(xsdata.convert_to_continuous_energy())
            data[-1].export_to_hdf5(h5_filename)

        library.register_file(h5_filename)
        library.export_to_xml(library_filename)

    def export_to_hdf5(self, filename='mgxs.h5'):
        """Create an hdf5 file that can be used for a simulation.

        Parameters
        ----------
        filename : str
            Filename of file, default is mgxs.h5.

        """

        check_type('filename', filename, string_types)

        # Create and write to the HDF5 file
        file = h5py.File(filename, "w")
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

        group_structure = file.attrs['group structure']
        num_delayed_groups = file.attrs['delayed_groups']
        energy_groups = openmc.mgxs.EnergyGroups(group_structure)
        data = cls(energy_groups, num_delayed_groups)

        for group_name, group in file.items():
            data.add_xsdata(openmc.XSdata.from_hdf5(group, group_name,
                                                    energy_groups,
                                                    num_delayed_groups))

        return data
