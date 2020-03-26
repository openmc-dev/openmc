from collections.abc import Iterable
from collections import namedtuple
from difflib import get_close_matches
from numbers import Real
from io import StringIO
import itertools
import os
import re
import tempfile
from warnings import warn

import numpy as np
import h5py

import openmc.checkvalue as cv
from openmc.mixin import EqualityMixin
from openmc.stats import Discrete, Tabular
from . import HDF5_VERSION, HDF5_VERSION_MAJOR, endf
from .data import K_BOLTZMANN, ATOMIC_SYMBOL, EV_PER_MEV, NATURAL_ABUNDANCE
from .ace import Table, get_table, Library
from .angle_energy import AngleEnergy
from .function import Tabulated1D, Function1D
from .njoy import make_ace_thermal
from .thermal_angle_energy import (CoherentElasticAE, IncoherentElasticAE,
                                   IncoherentElasticAEDiscrete,
                                   IncoherentInelasticAEDiscrete,
                                   IncoherentInelasticAE)


_THERMAL_NAMES = {
    'c_Al27': ('al', 'al27', 'al-27'),
    'c_Al_in_Sapphire': ('asap00',),
    'c_Be': ('be', 'be-metal', 'be-met', 'be00'),
    'c_BeO': ('beo',),
    'c_Be_in_BeO': ('bebeo', 'be-beo', 'be-o', 'be/o', 'bbeo00'),
    'c_Be_in_Be2C': ('bebe2c',),
    'c_C6H6': ('benz', 'c6h6'),
    'c_C_in_SiC': ('csic', 'c-sic'),
    'c_Ca_in_CaH2': ('cah', 'cah00'),
    'c_D_in_D2O': ('dd2o', 'd-d2o', 'hwtr', 'hw', 'dhw00'),
    'c_D_in_D2O_ice': ('dice',),
    'c_Fe56': ('fe', 'fe56', 'fe-56'),
    'c_Graphite': ('graph', 'grph', 'gr', 'gr00'),
    'c_Graphite_10p': ('grph10',),
    'c_Graphite_30p': ('grph30',),
    'c_H_in_CaH2': ('hcah2', 'hca00'),
    'c_H_in_CH2': ('hch2', 'poly', 'pol', 'h-poly', 'pol00'),
    'c_H_in_CH4_liquid': ('lch4', 'lmeth'),
    'c_H_in_CH4_solid': ('sch4', 'smeth'),
    'c_H_in_CH4_solid_phase_II': ('sch4p2',),
    'c_H_in_H2O': ('hh2o', 'h-h2o', 'lwtr', 'lw', 'lw00'),
    'c_H_in_H2O_solid': ('hice', 'h-ice', 'ice00'),
    'c_H_in_C5O2H8': ('lucite', 'c5o2h8', 'h-luci'),
    'c_H_in_Mesitylene': ('mesi00',),
    'c_H_in_Toluene': ('tol00',),
    'c_H_in_YH2': ('hyh2', 'h-yh2'),
    'c_H_in_ZrH': ('hzrh', 'h-zrh', 'h-zr', 'h/zr', 'hzr', 'hzr00'),
    'c_Mg24': ('mg', 'mg24', 'mg00'),
    'c_O_in_Sapphire': ('osap00',),
    'c_O_in_BeO': ('obeo', 'o-beo', 'o-be', 'o/be', 'obeo00'),
    'c_O_in_D2O': ('od2o', 'o-d2o', 'ohw00'),
    'c_O_in_H2O_ice': ('oice', 'o-ice'),
    'c_O_in_UO2': ('ouo2', 'o-uo2', 'o2-u', 'o2/u', 'ouo200'),
    'c_N_in_UN': ('n-un',),
    'c_ortho_D': ('orthod', 'orthoD', 'dortho', 'od200'),
    'c_ortho_H': ('orthoh', 'orthoH', 'hortho', 'oh200'),
    'c_Si28': ('si00',),
    'c_Si_in_SiC': ('sisic', 'si-sic'),
    'c_SiO2_alpha': ('sio2', 'sio2a'),
    'c_SiO2_beta': ('sio2b',),
    'c_para_D': ('parad', 'paraD', 'dpara', 'pd200'),
    'c_para_H': ('parah', 'paraH', 'hpara', 'ph200'),
    'c_U_in_UN': ('u-un',),
    'c_U_in_UO2': ('uuo2', 'u-uo2', 'u-o2', 'u/o2', 'uuo200'),
    'c_Y_in_YH2': ('yyh2', 'y-yh2'),
    'c_Zr_in_ZrH': ('zrzrh', 'zr-zrh', 'zr-h', 'zr/h')
}


def _temperature_str(T):
    # round() normally returns an int when called with a single argument, but
    # numpy floats overload rounding to return another float
    return "{}K".format(int(round(T)))


def get_thermal_name(name):
    """Get proper S(a,b) table name, e.g. 'HH2O' -> 'c_H_in_H2O'

    Parameters
    ----------
    name : str
        Name of an ACE thermal scattering table

    Returns
    -------
    str
        GND-format thermal scattering name

    """
    if name in _THERMAL_NAMES:
        return name
    else:
        for proper_name, names in _THERMAL_NAMES.items():
            if name.lower() in names:
                return proper_name

        # Make an educated guess?? This actually works well for
        # JEFF-3.2 which stupidly uses names like lw00.32t,
        # lw01.32t, etc. for different temperatures

        # First, construct a list of all the values/keys in the names
        # dictionary
        all_names = itertools.chain(_THERMAL_NAMES.keys(),
                                    *_THERMAL_NAMES.values())

        matches = get_close_matches(name, all_names, cutoff=0.5)
        if matches:
            # Figure out the key for the corresponding match
            match = matches[0]
            if match not in _THERMAL_NAMES:
                for key, value_list in _THERMAL_NAMES.items():
                    if match in value_list:
                        match = key
                        break

            warn('Thermal scattering material "{}" is not recognized. '
                 'Assigning a name of {}.'.format(name, match))
            return match
        else:
            # OK, we give up. Just use the ACE name.
            warn('Thermal scattering material "{0}" is not recognized. '
                 'Assigning a name of c_{0}.'.format(name))
            return 'c_' + name


class CoherentElastic(Function1D):
    r"""Coherent elastic scattering data from a crystalline material

    The integrated cross section for coherent elastic scattering from a
    powdered crystalline material may be represented as:

    .. math::
        \sigma(E,T) = \frac{1}{E} \sum\limits_{i=1}^{E_i < E} s_i(T)

    where :math:`s_i(T)` is proportional the structure factor in [eV-b] at
    the moderator temperature :math:`T` in Kelvin.

    Parameters
    ----------
    bragg_edges : Iterable of float
        Bragg edge energies in eV
    factors : Iterable of float
        Partial sum of structure factors, :math:`\sum\limits_{i=1}^{E_i<E} s_i`

    Attributes
    ----------
    bragg_edges : Iterable of float
        Bragg edge energies in eV
    factors : Iterable of float
        Partial sum of structure factors, :math:`\sum\limits_{i=1}^{E_i<E} s_i`

    """

    def __init__(self, bragg_edges, factors):
        self.bragg_edges = bragg_edges
        self.factors = factors

    def __call__(self, E):
        idx = np.searchsorted(self.bragg_edges, E) - 1
        if isinstance(E, Iterable):
            E = np.asarray(E)
            nonzero = idx >= 0
            xs = np.zeros_like(E)
            xs[nonzero] = self.factors[idx[nonzero]] / E[nonzero]
            return xs
        else:
            return self.factors[idx] / E if idx >= 0 else 0.0

    def __len__(self):
        return len(self.bragg_edges)

    @property
    def bragg_edges(self):
        return self._bragg_edges

    @property
    def factors(self):
        return self._factors

    @bragg_edges.setter
    def bragg_edges(self, bragg_edges):
        cv.check_type('Bragg edges', bragg_edges, Iterable, Real)
        self._bragg_edges = np.asarray(bragg_edges)

    @factors.setter
    def factors(self, factors):
        cv.check_type('structure factor cumulative sums', factors,
                      Iterable, Real)
        self._factors = np.asarray(factors)

    def to_hdf5(self, group, name):
        """Write coherent elastic scattering to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to
        name : str
            Name of the dataset to create

        """
        dataset = group.create_dataset(name, data=np.vstack(
            [self.bragg_edges, self.factors]))
        dataset.attrs['type'] = np.string_(type(self).__name__)

    @classmethod
    def from_hdf5(cls, dataset):
        """Read coherent elastic scattering from an HDF5 dataset

        Parameters
        ----------
        dataset : h5py.Dataset
            HDF5 dataset to read from

        Returns
        -------
        openmc.data.CoherentElastic
            Coherent elastic scattering cross section

        """
        bragg_edges = dataset[0, :]
        factors = dataset[1, :]
        return cls(bragg_edges, factors)


class IncoherentElastic(Function1D):
    r"""Incoherent elastic scattering cross section

    Elastic scattering can be treated in the incoherent approximation for
    partially ordered systems such as ZrHx and polyethylene. The integrated
    cross section can be obtained as:

    .. math::
        \sigma(E,T) = \frac{\sigma_b}{2} \left ( \frac{1 - e^{-4EW'(T)}}
        {2EW'(T)} \right )

    where :math:`\sigma_b` is the characteristic bound cross section, and
    :math:`W'(T)` is the Debye-Waller integral divided by the atomic mass
    in [eV\ :math:`^{-1}`].

    Parameters
    ----------
    bound_xs : float
        Characteristic bound cross section in [b]
    debye_waller : float
        Debye-Waller integral in [eV\ :math:`^{-1}`]

    Attributes
    ----------
    bound_xs : float
        Characteristic bound cross section in [b]
    debye_waller : float
        Debye-Waller integral in [eV\ :math:`^{-1}`]

    """
    def __init__(self, bound_xs, debye_waller):
        self.bound_xs = bound_xs
        self.debye_waller = debye_waller

    def __call__(self, E):
        W = self.debye_waller
        return self.bound_xs / 2.0 * (1 - np.exp(-4*E*W)) / (2*E*W)

    def to_hdf5(self, group, name):
        """Write incoherent elastic scattering to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to
        name : str
            Name of the dataset to create

        """
        data = np.array([self.bound_xs, self.debye_waller])
        dataset = group.create_dataset(name, data=data)
        dataset.attrs['type'] = np.string_(type(self).__name__)

    @classmethod
    def from_hdf5(cls, dataset):
        """Read incoherent elastic scattering from an HDF5 dataset

        Parameters
        ----------
        dataset : h5py.Dataset
            HDF5 dataset to read from

        Returns
        -------
        openmc.data.IncoherentElastic
            Incoherent elastic scattering cross section

        """
        bound_xs, debye_waller = dataset[()]
        return cls(bound_xs, debye_waller)


class ThermalScatteringReaction(EqualityMixin):
    r"""Thermal scattering reaction

    This class is used to hold the integral and differential cross sections
    for either elastic or inelastic thermal scattering.

    Parameters
    ----------
    xs : dict of str to Function1D
        Integral cross section at each temperature
    distribution : dict of str to AngleEnergy
        Secondary angle-energy distribution at each temperature

    Attributes
    ----------
    xs : dict of str to Function1D
        Integral cross section at each temperature
    distribution : dict of str to AngleEnergy
        Secondary angle-energy distribution at each temperature

    """
    def __init__(self, xs, distribution):
        self.xs = xs
        self.distribution = distribution

    def to_hdf5(self, group, name):
        """Write thermal scattering reaction to HDF5

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to
        name : {'elastic', 'inelastic'}
            Name of reaction to write

        """
        for T, xs in self.xs.items():
            Tgroup = group.require_group(T)
            rx_group = Tgroup.create_group(name)
            xs.to_hdf5(rx_group, 'xs')
            dgroup = rx_group.create_group('distribution')
            self.distribution[T].to_hdf5(dgroup)

    @classmethod
    def from_hdf5(cls, group, name, temperatures):
        """Generate thermal scattering reaction data from HDF5

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from
        name : {'elastic', 'inelastic'}
            Name of the reaction to read
        temperatures : Iterable of str
            Temperatures to read

        Returns
        -------
        openmc.data.ThermalScatteringReaction
            Thermal scattering reaction data

        """
        xs = {}
        distribution = {}
        for T in temperatures:
            rx_group = group[T][name]
            xs[T] = Function1D.from_hdf5(rx_group['xs'])
            if isinstance(xs[T], CoherentElastic):
                distribution[T] = CoherentElasticAE(xs[T])
            else:
                distribution[T] = AngleEnergy.from_hdf5(rx_group['distribution'])
        return cls(xs, distribution)


class ThermalScattering(EqualityMixin):
    """A ThermalScattering object contains thermal scattering data as represented by
    an S(alpha, beta) table.

    Parameters
    ----------
    name : str
        Name of the material using GND convention, e.g. c_H_in_H2O
    atomic_weight_ratio : float
        Atomic mass ratio of the target nuclide.
    kTs : Iterable of float
        List of temperatures of the target nuclide in the data set.
        The temperatures have units of eV.

    Attributes
    ----------
    atomic_weight_ratio : float
        Atomic mass ratio of the target nuclide.
    energy_max : float
        Maximum energy for thermal scattering data in [eV]
    elastic : openmc.data.ThermalScatteringReaction or None
        Elastic scattering derived in the coherent or incoherent approximation
    inelastic : openmc.data.ThermalScatteringReaction
        Inelastic scattering cross section derived in the incoherent
        approximation
    name : str
        Name of the material using GND convention, e.g. c_H_in_H2O
    temperatures : Iterable of str
        List of string representations the temperatures of the target nuclide
        in the data set.  The temperatures are strings of the temperature,
        rounded to the nearest integer; e.g., '294K'
    kTs : Iterable of float
        List of temperatures of the target nuclide in the data set.
        The temperatures have units of eV.
    nuclides : Iterable of str
        Nuclide names that the thermal scattering data applies to

    """

    def __init__(self, name, atomic_weight_ratio, energy_max, kTs):
        self.name = name
        self.atomic_weight_ratio = atomic_weight_ratio
        self.energy_max = energy_max
        self.kTs = kTs
        self.elastic = None
        self.inelastic = None
        self.nuclides = []

    def __repr__(self):
        if hasattr(self, 'name'):
            return "<Thermal Scattering Data: {}>".format(self.name)
        else:
            return "<Thermal Scattering Data>"

    @property
    def temperatures(self):
        return [_temperature_str(kT / K_BOLTZMANN) for kT in self.kTs]

    def export_to_hdf5(self, path, mode='a', libver='earliest'):
        """Export table to an HDF5 file.

        Parameters
        ----------
        path : str
            Path to write HDF5 file to
        mode : {'r', r+', 'w', 'x', 'a'}
            Mode that is used to open the HDF5 file. This is the second argument
            to the :class:`h5py.File` constructor.
        libver : {'earliest', 'latest'}
            Compatibility mode for the HDF5 file. 'latest' will produce files
            that are less backwards compatible but have performance benefits.

        """
        # Open file and write version
        with h5py.File(str(path), mode, libver=libver) as f:
            f.attrs['filetype'] = np.string_('data_thermal')
            f.attrs['version'] = np.array(HDF5_VERSION)

            # Write basic data
            g = f.create_group(self.name)
            g.attrs['atomic_weight_ratio'] = self.atomic_weight_ratio
            g.attrs['energy_max'] = self.energy_max
            g.attrs['nuclides'] = np.array(self.nuclides, dtype='S')
            ktg = g.create_group('kTs')
            for i, temperature in enumerate(self.temperatures):
                ktg.create_dataset(temperature, data=self.kTs[i])

            # Write elastic/inelastic reaction data
            if self.elastic is not None:
                self.elastic.to_hdf5(g, 'elastic')
            self.inelastic.to_hdf5(g, 'inelastic')

    def add_temperature_from_ace(self, ace_or_filename, name=None):
        """Add data to the ThermalScattering object from an ACE file at a
        different temperature.

        Parameters
        ----------
        ace_or_filename : openmc.data.ace.Table or str
            ACE table to read from. If given as a string, it is assumed to be
            the filename for the ACE file.
        name : str
            GND-conforming name of the material, e.g. c_H_in_H2O. If none is
            passed, the appropriate name is guessed based on the name of the ACE
            table.

        Returns
        -------
        openmc.data.ThermalScattering
            Thermal scattering data

        """
        data = ThermalScattering.from_ace(ace_or_filename, name)

        # Check if temprature already exists
        strT = data.temperatures[0]
        if strT in self.temperatures:
            warn('S(a,b) data at T={} already exists.'.format(strT))
            return

        # Check that name matches
        if data.name != self.name:
            raise ValueError('Data provided for an incorrect material.')

        # Add temperature
        self.kTs += data.kTs

        # Add inelastic cross section and distributions
        if data.inelastic is not None:
            self.inelastic.xs.update(data.inelastic.xs)
            self.inelastic.distribution.update(data.inelastic.distribution)

        # Add elastic cross sectoin and angular distribution
        if data.elastic is not None:
            self.elastic.xs.update(data.elastic.xs)
            self.elastic.distribution.update(data.elastic.distribution)

    @classmethod
    def from_hdf5(cls, group_or_filename):
        """Generate thermal scattering data from HDF5 group

        Parameters
        ----------
        group_or_filename : h5py.Group or str
            HDF5 group containing interaction data. If given as a string, it is
            assumed to be the filename for the HDF5 file, and the first group
            is used to read from.

        Returns
        -------
        openmc.data.ThermalScattering
            Neutron thermal scattering data

        """
        if isinstance(group_or_filename, h5py.Group):
            group = group_or_filename
        else:
            h5file = h5py.File(str(group_or_filename), 'r')

            # Make sure version matches
            if 'version' in h5file.attrs:
                major, minor = h5file.attrs['version']
                if major != HDF5_VERSION_MAJOR:
                    raise IOError(
                        'HDF5 data format uses version {}.{} whereas your '
                        'installation of the OpenMC Python API expects version '
                        '{}.x.'.format(major, minor, HDF5_VERSION_MAJOR))
            else:
                raise IOError(
                    'HDF5 data does not indicate a version. Your installation of '
                    'the OpenMC Python API expects version {}.x data.'
                    .format(HDF5_VERSION_MAJOR))

            group = list(h5file.values())[0]

        name = group.name[1:]
        atomic_weight_ratio = group.attrs['atomic_weight_ratio']
        energy_max = group.attrs['energy_max']
        kTg = group['kTs']
        kTs = [dataset[()] for dataset in kTg.values()]

        table = cls(name, atomic_weight_ratio, energy_max, kTs)
        table.nuclides = [nuc.decode() for nuc in group.attrs['nuclides']]

        # Read thermal elastic scattering
        if 'elastic' in group[table.temperatures[0]]:
            table.elastic = ThermalScatteringReaction.from_hdf5(
                group, 'elastic', table.temperatures
            )

        # Read thermal inelastic scattering
        table.inelastic = ThermalScatteringReaction.from_hdf5(
            group, 'inelastic', table.temperatures
        )

        return table

    @classmethod
    def from_ace(cls, ace_or_filename, name=None):
        """Generate thermal scattering data from an ACE table

        Parameters
        ----------
        ace_or_filename : openmc.data.ace.Table or str
            ACE table to read from. If given as a string, it is assumed to be
            the filename for the ACE file.
        name : str
            GND-conforming name of the material, e.g. c_H_in_H2O. If none is
            passed, the appropriate name is guessed based on the name of the ACE
            table.

        Returns
        -------
        openmc.data.ThermalScattering
            Thermal scattering data

        """
        if isinstance(ace_or_filename, Table):
            ace = ace_or_filename
        else:
            ace = get_table(ace_or_filename)

        # Get new name that is GND-consistent
        ace_name, xs = ace.name.split('.')
        if not xs.endswith('t'):
            raise TypeError("{} is not a thermal scattering ACE table.".format(ace))
        name = get_thermal_name(ace_name)

        # Assign temperature to the running list
        kTs = [ace.temperature*EV_PER_MEV]

        # Incoherent inelastic scattering cross section
        idx = ace.jxs[1]
        n_energy = int(ace.xss[idx])
        energy = ace.xss[idx+1 : idx+1+n_energy]*EV_PER_MEV
        xs = ace.xss[idx+1+n_energy : idx+1+2*n_energy]
        inelastic_xs = Tabulated1D(energy, xs)
        energy_max = energy[-1]

        # Incoherent inelastic angle-energy distribution
        continuous = (ace.nxs[7] == 2)
        n_energy_out = ace.nxs[4]
        if not continuous:
            n_mu = ace.nxs[3]
            idx = ace.jxs[3]
            energy_out = ace.xss[idx:idx + n_energy * n_energy_out *
                (n_mu + 2): n_mu + 2]*EV_PER_MEV
            energy_out.shape = (n_energy, n_energy_out)

            mu_out = ace.xss[idx:idx + n_energy * n_energy_out * (n_mu + 2)]
            mu_out.shape = (n_energy, n_energy_out, n_mu+2)
            mu_out = mu_out[:, :, 1:]
            skewed = (ace.nxs[7] == 1)
            distribution = IncoherentInelasticAEDiscrete(energy_out, mu_out, skewed)
        else:
            n_mu = ace.nxs[3] - 1
            idx = ace.jxs[3]
            locc = ace.xss[idx:idx + n_energy].astype(int)
            n_energy_out = \
                ace.xss[idx + n_energy:idx + 2 * n_energy].astype(int)
            energy_out = []
            mu_out = []
            for i in range(n_energy):
                idx = locc[i]

                # Outgoing energy distribution for incoming energy i
                e = ace.xss[idx + 1:idx + 1 + n_energy_out[i]*(n_mu + 3):
                            n_mu + 3]*EV_PER_MEV
                p = ace.xss[idx + 2:idx + 2 + n_energy_out[i]*(n_mu + 3):
                            n_mu + 3]/EV_PER_MEV
                c = ace.xss[idx + 3:idx + 3 + n_energy_out[i]*(n_mu + 3):
                            n_mu + 3]
                eout_i = Tabular(e, p, 'linear-linear', ignore_negative=True)
                eout_i.c = c

                # Outgoing angle distribution for each
                # (incoming, outgoing) energy pair
                mu_i = []
                for j in range(n_energy_out[i]):
                    mu = ace.xss[idx + 4:idx + 4 + n_mu]
                    p_mu = 1. / n_mu * np.ones(n_mu)
                    mu_ij = Discrete(mu, p_mu)
                    mu_ij.c = np.cumsum(p_mu)
                    mu_i.append(mu_ij)
                    idx += 3 + n_mu

                energy_out.append(eout_i)
                mu_out.append(mu_i)

            # Create correlated angle-energy distribution
            breakpoints = [n_energy]
            interpolation = [2]
            energy = inelastic_xs.x
            distribution = IncoherentInelasticAE(
                breakpoints, interpolation, energy, energy_out, mu_out)

        table = cls(name, ace.atomic_weight_ratio, energy_max, kTs)
        T = table.temperatures[0]
        table.inelastic = ThermalScatteringReaction(
            {T: inelastic_xs}, {T: distribution}
        )

        # Incoherent/coherent elastic scattering cross section
        idx = ace.jxs[4]
        n_mu = ace.nxs[6] + 1
        if idx != 0:
            n_energy = int(ace.xss[idx])
            energy = ace.xss[idx + 1: idx + 1 + n_energy]*EV_PER_MEV
            P = ace.xss[idx + 1 + n_energy: idx + 1 + 2 * n_energy]

            if ace.nxs[5] == 4:
                # Coherent elastic
                xs = CoherentElastic(energy, P*EV_PER_MEV)
                distribution = CoherentElasticAE(xs)

                # Coherent elastic shouldn't have angular distributions listed
                assert n_mu == 0
            else:
                # Incoherent elastic
                xs = Tabulated1D(energy, P)

                # Angular distribution
                assert n_mu > 0
                idx = ace.jxs[6]
                mu_out = ace.xss[idx:idx + n_energy * n_mu]
                mu_out.shape = (n_energy, n_mu)
                distribution = IncoherentElasticAEDiscrete(mu_out)

            table.elastic = ThermalScatteringReaction({T: xs}, {T: distribution})

        # Get relevant nuclides -- NJOY only allows one to specify three
        # nuclides that the S(a,b) table applies to. Thus, for all elements
        # other than H and Fe, we automatically add all the naturally-occurring
        # isotopes.
        for zaid, awr in ace.pairs:
            if zaid > 0:
                Z, A = divmod(zaid, 1000)
                element = ATOMIC_SYMBOL[Z]
                if element in ['H', 'Fe']:
                    table.nuclides.append(element + str(A))
                else:
                    if element + '0' not in table.nuclides:
                        table.nuclides.append(element + '0')
                    for isotope in sorted(NATURAL_ABUNDANCE):
                        if re.match(r'{}\d+'.format(element), isotope):
                            if isotope not in table.nuclides:
                                table.nuclides.append(isotope)

        return table

    @classmethod
    def from_njoy(cls, filename, filename_thermal, temperatures=None,
                  evaluation=None, evaluation_thermal=None,
                  use_endf_data=True, **kwargs):
        """Generate thermal scattering data by running NJOY.

        Parameters
        ----------
        filename : str
            Path to ENDF neutron sublibrary file
        filename_thermal : str
            Path to ENDF thermal scattering sublibrary file
        temperatures : iterable of float
            Temperatures in Kelvin to produce data at. If omitted, data is
            produced at all temperatures in the ENDF thermal scattering
            sublibrary.
        evaluation : openmc.data.endf.Evaluation, optional
            If the ENDF neutron sublibrary file contains multiple material
            evaluations, this argument indicates which evaluation to use.
        evaluation_thermal : openmc.data.endf.Evaluation, optional
            If the ENDF thermal scattering sublibrary file contains multiple
            material evaluations, this argument indicates which evaluation to
            use.
        use_endf_data : bool
            If the material has incoherent elastic scattering, the ENDF data
            will be used rather than the ACE data.
        **kwargs
            Keyword arguments passed to :func:`openmc.data.njoy.make_ace_thermal`

        Returns
        -------
        data : openmc.data.ThermalScattering
            Thermal scattering data

        """
        with tempfile.TemporaryDirectory() as tmpdir:
            # Run NJOY to create an ACE library
            kwargs.setdefault('output_dir', tmpdir)
            kwargs.setdefault('ace', os.path.join(kwargs['output_dir'], 'ace'))
            kwargs['evaluation'] = evaluation
            kwargs['evaluation_thermal'] = evaluation_thermal
            make_ace_thermal(filename, filename_thermal, temperatures, **kwargs)

            # Create instance from ACE tables within library
            lib = Library(kwargs['ace'])
            data = cls.from_ace(lib.tables[0])
            for table in lib.tables[1:]:
                data.add_temperature_from_ace(table)

            # Load ENDF data to replace incoherent elastic
            if use_endf_data:
                data_endf = cls.from_endf(filename_thermal)
                if data_endf.elastic is not None:
                    # Get appropriate temperatures
                    if temperatures is None:
                        temperatures = data_endf.temperatures
                    else:
                        temperatures = [_temperature_str(t) for t in temperatures]

                    # Replace ACE data with ENDF data
                    rx, rx_endf = data.elastic, data_endf.elastic
                    for t in temperatures:
                        if isinstance(rx_endf.xs[t], IncoherentElastic):
                            rx.xs[t] = rx_endf.xs[t]
                            rx.distribution[t] = rx_endf.distribution[t]

        return data

    @classmethod
    def from_endf(cls, ev_or_filename):
        """Generate thermal scattering data from an ENDF file

        Parameters
        ----------
        ev_or_filename : openmc.data.endf.Evaluation or str
            ENDF evaluation to read from. If given as a string, it is assumed to
            be the filename for the ENDF file.

        Returns
        -------
        openmc.data.ThermalScattering
            Thermal scattering data

        """
        if isinstance(ev_or_filename, endf.Evaluation):
            ev = ev_or_filename
        else:
            ev = endf.Evaluation(ev_or_filename)

        # Read coherent/incoherent elastic data
        elastic = None
        if (7, 2) in ev.section:
            xs = {}
            distribution = {}

            file_obj = StringIO(ev.section[7, 2])
            lhtr = endf.get_head_record(file_obj)[2]
            if lhtr == 1:
                # coherent elastic

                # Get structure factor at first temperature
                params, S = endf.get_tab1_record(file_obj)
                strT = _temperature_str(params[0])
                n_temps = params[2]
                bragg_edges = S.x
                xs[strT] = CoherentElastic(bragg_edges, S.y)
                distribution = {strT: CoherentElasticAE(xs[strT])}

                # Get structure factor for subsequent temperatures
                for _ in range(n_temps):
                    params, S = endf.get_list_record(file_obj)
                    strT = _temperature_str(params[0])
                    xs[strT] = CoherentElastic(bragg_edges, S)
                    distribution[strT] = CoherentElasticAE(xs[strT])

            elif lhtr == 2:
                # incoherent elastic
                params, W = endf.get_tab1_record(file_obj)
                bound_xs = params[0]
                for T, debye_waller in zip(W.x, W.y):
                    strT = _temperature_str(T)
                    xs[strT] = IncoherentElastic(bound_xs, debye_waller)
                    distribution[strT] = IncoherentElasticAE(debye_waller)

            elastic = ThermalScatteringReaction(xs, distribution)

        # Read incoherent inelastic data
        assert (7, 4) in ev.section, 'No MF=7, MT=4 found in thermal scattering'
        file_obj = StringIO(ev.section[7, 4])
        params = endf.get_head_record(file_obj)
        data = {'symmetric': params[4] == 0}

        # Get information about principal atom
        params, B = endf.get_list_record(file_obj)
        data['log'] = bool(params[2])
        data['free_atom_xs'] = B[0]
        data['epsilon'] = B[1]
        data['A0'] = awr = B[2]
        data['e_max'] = energy_max = B[3]
        data['M0'] = B[5]

        # Get information about non-principal atoms
        n_non_principal = params[5]
        data['non_principal'] = []
        NonPrincipal = namedtuple('NonPrincipal', ['func', 'xs', 'A', 'M'])
        for i in range(1, n_non_principal + 1):
            func = {0.0: 'SCT', 1.0: 'free gas', 2.0: 'diffusive'}[B[6*i]]
            xs = B[6*i + 1]
            A = B[6*i + 2]
            M = B[6*i + 5]
            data['non_principal'].append(NonPrincipal(func, xs, A, M))

        # Get S(alpha,beta,T)
        kTs = []
        if data['free_atom_xs'] > 0.0:
            params, _ = endf.get_tab2_record(file_obj)
            n_beta = params[5]
            sab = {'beta': np.empty(n_beta)}
            for i in range(n_beta):
                params, S = endf.get_tab1_record(file_obj)
                T0, beta, lt = params[:3]
                if i == 0:
                    sab['alpha'] = alpha = S.x
                    sab[T0] = np.empty((alpha.size, n_beta))
                    kTs.append(K_BOLTZMANN * T0)
                sab['beta'][i] = beta
                sab[T0][:, i] = S.y
                for _ in range(lt):
                    params, S = endf.get_list_record(file_obj)
                    T = params[0]
                    if i == 0:
                        sab[T] = np.empty((alpha.size, n_beta))
                        kTs.append(K_BOLTZMANN * T)
                    sab[T][:, i] = S
            data['sab'] = sab

        # Get effective temperature for each atom
        _, Teff = endf.get_tab1_record(file_obj)
        data['effective_temperature'] = [Teff]
        for atom in data['non_principal']:
            if atom.func == 'SCT':
                _, Teff = endf.get_tab1_record(file_obj)
                data['effective_temperature'].append(Teff)

        name = ev.target['zsymam'].strip()
        instance = cls(name, awr, energy_max, kTs)
        if elastic is not None:
            instance.elastic = elastic

        # Currently we don't have a proper cross section or distribution for
        # incoherent inelastic, so we just create an empty object and attach
        # all the data as a dictionary
        instance.inelastic = ThermalScatteringReaction(None, None)
        instance.inelastic.data = data

        return instance
