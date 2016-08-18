from collections import Iterable
from difflib import get_close_matches
from numbers import Real
from warnings import warn

import numpy as np
import h5py

import openmc.checkvalue as cv
from .ace import Table, get_table
from .angle_energy import AngleEnergy
from .function import Tabulated1D
from .correlated import CorrelatedAngleEnergy
from openmc.stats import Discrete, Tabular


_THERMAL_NAMES = {'al': 'c_Al27', 'al27': 'c_Al27',
                  'be': 'c_Be',
                  'beo': 'c_BeO',
                  'bebeo': 'c_Be_in_BeO', 'be-o': 'c_Be_in_BeO', 'be/o': 'c_Be_in_BeO',
                  'benz': 'c_Benzine',
                  'cah': 'c_Ca_in_CaH2',
                  'dd2o': 'c_D_in_D2O', 'hwtr': 'c_D_in_D2O',
                  'fe': 'c_Fe56', 'fe56': 'c_Fe56',
                  'graph': 'c_Graphite', 'grph': 'c_Graphite',
                  'hca': 'c_H_in_CaH2',
                  'hch2': 'c_H_in_CH2', 'poly': 'c_H_in_CH2',
                  'hh2o': 'c_H_in_H2O', 'lwtr': 'c_H_in_H2O',
                  'hzrh': 'c_H_in_ZrH', 'h-zr': 'c_H_in_ZrH', 'h/zr': 'c_H_in_ZrH',
                  'lch4': 'c_liquid_CH4', 'lmeth': 'c_liquid_CH4',
                  'mg': 'c_Mg24',
                  'obeo': 'c_O_in_BeO', 'o-be': 'c_O_in_BeO', 'o/be': 'c_O_in_BeO',
                  'orthod': 'c_ortho_D', 'dortho': 'c_ortho_D',
                  'orthoh': 'c_ortho_H', 'hortho': 'c_ortho_H',
                  'ouo2': 'c_O_in_UO2', 'o2-u': 'c_O_in_UO2', 'o2/u': 'c_O_in_UO2',
                  'sio2': 'c_SiO2',
                  'parad': 'c_para_D', 'dpara': 'c_para_D',
                  'parah': 'c_para_H', 'hpara': 'c_para_H',
                  'sch4': 'c_solid_CH4', 'smeth': 'c_solid_CH4',
                  'uuo2': 'c_U_in_UO2', 'u-o2': 'c_U_in_UO2', 'u/o2': 'c_U_in_UO2',
                  'zrzrh': 'c_Zr_in_ZrH', 'zr-h': 'c_Zr_in_ZrH', 'zr/h': 'c_Zr_in_ZrH'}


def get_thermal_name(name):
    """Get proper S(a,b) table name, e.g. 'HH2O' -> 'c_H_in_H2O'"""

    if name in _THERMAL_NAMES.values():
        return name
    elif name.lower() in _THERMAL_NAMES:
        return _THERMAL_NAMES[name.lower()]
    else:
        # Make an educated guess?? This actually works well for
        # JEFF-3.2 which stupidly uses names like lw00.32t,
        # lw01.32t, etc. for different temperatures
        matches = get_close_matches(
            name.lower(), _THERMAL_NAMES.keys(), cutoff=0.5)
        if len(matches) > 0:
            return _THERMAL_NAMES[matches[0]]
        else:
            # OK, we give up. Just use the ACE name.
            return 'c_' + name


class CoherentElastic(object):
    r"""Coherent elastic scattering data from a crystalline material

    Parameters
    ----------
    bragg_edges : Iterable of float
        Bragg edge energies in MeV
    factors : Iterable of float
        Partial sum of structure factors, :math:`\sum\limits_{i=1}^{E_i<E} S_i`

    Attributes
    ----------
    bragg_edges : Iterable of float
        Bragg edge energies in MeV
    factors : Iterable of float
        Partial sum of structure factors, :math:`\sum\limits_{i=1}^{E_i<E} S_i`

    """

    def __init__(self, bragg_edges, factors):
        self.bragg_edges = bragg_edges
        self.factors = factors

    def __call__(self, E):
        if isinstance(E, Iterable):
            E = np.asarray(E)
        idx = np.searchsorted(self.bragg_edges, E)
        return self.factors[idx]/E

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
        dataset.attrs['type'] = np.string_('bragg')

    @classmethod
    def from_hdf5(cls, dataset):
        """Read coherent elastic scattering from an HDF5 dataset

        Parameters
        ----------
        group : h5py.Dataset
            HDF5 group to write to

        Returns
        -------
        openmc.data.CoherentElastic
            Coherent elastic scattering cross section

        """
        bragg_edges = dataset.value[0, :]
        factors = dataset.value[1, :]
        return cls(bragg_edges, factors)


class ThermalScattering(object):
    """A ThermalScattering object contains thermal scattering data as represented by
    an S(alpha, beta) table.

    Parameters
    ----------
    name : str
        ZAID identifier of the table, e.g. lwtr.10t.
    atomic_weight_ratio : float
        Atomic mass ratio of the target nuclide.
    temperature : float
        Temperature of the target nuclide in eV.

    Attributes
    ----------
    atomic_weight_ratio : float
        Atomic mass ratio of the target nuclide.
    elastic_xs : openmc.data.Tabulated1D or openmc.data.CoherentElastic
        Elastic scattering cross section derived in the coherent or incoherent
        approximation
    inelastic_xs : openmc.data.Tabulated1D
        Inelastic scattering cross section derived in the incoherent
        approximation
    name : str
        Name of the table, e.g. lwtr.20t.
    temperature : float
        Temperature of the target nuclide in eV.
    zaids : Iterable of int
        ZAID identifiers that the thermal scattering data applies to

    """

    def __init__(self, name, atomic_weight_ratio, temperature):
        self.name = name
        self.atomic_weight_ratio = atomic_weight_ratio
        self.temperature = temperature
        self.elastic_xs = None
        self.elastic_mu_out = None
        self.inelastic_xs = None
        self.inelastic_e_out = None
        self.inelastic_mu_out = None
        self.secondary_mode = None
        self.zaids = []

    def __repr__(self):
        if hasattr(self, 'name'):
            return "<Thermal Scattering Data: {0}>".format(self.name)
        else:
            return "<Thermal Scattering Data>"

    def export_to_hdf5(self, path, mode='a'):
        """Export table to an HDF5 file.

        Parameters
        ----------
        path : str
            Path to write HDF5 file to
        mode : {'r', r+', 'w', 'x', 'a'}
            Mode that is used to open the HDF5 file. This is the second argument
            to the :class:`h5py.File` constructor.

        """

        f = h5py.File(path, mode, libver='latest')

        # Write basic data
        g = f.create_group(self.name)
        g.attrs['atomic_weight_ratio'] = self.atomic_weight_ratio
        g.attrs['temperature'] = self.temperature
        g.attrs['zaids'] = self.zaids

        # Write thermal elastic scattering
        if self.elastic_xs is not None:
            elastic_group = g.create_group('elastic')
            self.elastic_xs.to_hdf5(elastic_group, 'xs')
            if self.elastic_mu_out is not None:
                elastic_group.create_dataset('mu_out', data=self.elastic_mu_out)

        # Write thermal inelastic scattering
        if self.inelastic_xs is not None:
            inelastic_group = g.create_group('inelastic')
            self.inelastic_xs.to_hdf5(inelastic_group, 'xs')
            inelastic_group.attrs['secondary_mode'] = np.string_(self.secondary_mode)
            if self.secondary_mode in ('equal', 'skewed'):
                inelastic_group.create_dataset('energy_out', data=self.inelastic_e_out)
                inelastic_group.create_dataset('mu_out', data=self.inelastic_mu_out)
            elif self.secondary_mode == 'continuous':
                self.inelastic_dist.to_hdf5(inelastic_group)

    @classmethod
    def from_hdf5(cls, group):
        """Generate thermal scattering data from HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from

        Returns
        -------
        openmc.data.ThermalScattering
            Neutron thermal scattering data

        """
        name = group.name[1:]
        atomic_weight_ratio = group.attrs['atomic_weight_ratio']
        temperature = group.attrs['temperature']
        table = cls(name, atomic_weight_ratio, temperature)
        table.zaids = group.attrs['zaids']

        # Read thermal elastic scattering
        if 'elastic' in group:
            elastic_group = group['elastic']

            # Cross section
            elastic_xs_type = elastic_group['xs'].attrs['type'].decode()
            if elastic_xs_type == 'tab1':
                table.elastic_xs = Tabulated1D.from_hdf5(elastic_group['xs'])
            elif elastic_xs_type == 'bragg':
                table.elastic_xs = CoherentElastic.from_hdf5(elastic_group['xs'])

            # Angular distribution
            if 'mu_out' in elastic_group:
                table.elastic_mu_out = elastic_group['mu_out'].value

        # Read thermal inelastic scattering
        if 'inelastic' in group:
            inelastic_group = group['inelastic']
            table.secondary_mode = inelastic_group.attrs['secondary_mode'].decode()
            table.inelastic_xs = Tabulated1D.from_hdf5(inelastic_group['xs'])
            if table.secondary_mode in ('equal', 'skewed'):
                table.inelastic_e_out = inelastic_group['energy_out']
                table.inelastic_mu_out = inelastic_group['mu_out']
            elif table.secondary_mode == 'continuous':
                table.inelastic_dist = AngleEnergy.from_hdf5(inelastic_group)

        return table

    @classmethod
    def from_ace(cls, ace_or_filename, name=None):
        """Generate thermal scattering data from an ACE table

        Parameters
        ----------
        ace : openmc.data.ace.Table or str
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
        if name is None:
            if ace_name.lower() in _THERMAL_NAMES:
                name = _THERMAL_NAMES[ace_name.lower()] + '.' + xs
            else:
                # Make an educated guess?? This actually works well for JEFF-3.2
                # which stupidly uses names like lw00.32t, lw01.32t, etc. for
                # different temperatures
                matches = get_close_matches(
                    ace_name.lower(), _THERMAL_NAMES.keys(), cutoff=0.5)
                if len(matches) > 0:
                    name = _THERMAL_NAMES[matches[0]] + '.' + xs
                else:
                    # OK, we give up. Just use the ACE name.
                    name = 'c_' + ace.name
                warn('Thermal scattering material "{}" is not recognized. '
                     'Assigning a name of {}.'.format(ace.name, name))

        table = cls(name, ace.atomic_weight_ratio, ace.temperature)

        # Incoherent inelastic scattering cross section
        idx = ace.jxs[1]
        n_energy = int(ace.xss[idx])
        energy = ace.xss[idx+1 : idx+1+n_energy]
        xs = ace.xss[idx+1+n_energy : idx+1+2*n_energy]
        table.inelastic_xs = Tabulated1D(energy, xs)

        if ace.nxs[7] == 0:
            table.secondary_mode = 'equal'
        elif ace.nxs[7] == 1:
            table.secondary_mode = 'skewed'
        elif ace.nxs[7] == 2:
            table.secondary_mode = 'continuous'

        n_energy_out = ace.nxs[4]
        if table.secondary_mode in ('equal', 'skewed'):
            n_mu = ace.nxs[3]
            idx = ace.jxs[3]
            table.inelastic_e_out = ace.xss[idx:idx+n_energy*n_energy_out*(n_mu+2):n_mu+2]
            table.inelastic_e_out.shape = (n_energy, n_energy_out)

            table.inelastic_mu_out = ace.xss[idx:idx+n_energy*n_energy_out*(n_mu+2)]
            table.inelastic_mu_out.shape = (n_energy, n_energy_out, n_mu+2)
            table.inelastic_mu_out = table.inelastic_mu_out[:, :, 1:]
        else:
            n_mu = ace.nxs[3] - 1
            idx = ace.jxs[3]
            locc = ace.xss[idx:idx + n_energy].astype(int)
            n_energy_out = ace.xss[idx + n_energy:idx + 2*n_energy].astype(int)
            energy_out = []
            mu_out = []
            for i in range(n_energy):
                idx = locc[i]

                # Outgoing energy distribution for incoming energy i
                e = ace.xss[idx + 1:idx + 1 + n_energy_out[i]*(n_mu + 3):n_mu + 3]
                p = ace.xss[idx + 2:idx + 2 + n_energy_out[i]*(n_mu + 3):n_mu + 3]
                c = ace.xss[idx + 3:idx + 3 + n_energy_out[i]*(n_mu + 3):n_mu + 3]
                eout_i = Tabular(e, p, 'linear-linear', ignore_negative=True)
                eout_i.c = c

                # Outgoing angle distribution for each (incoming, outgoing) energy pair
                mu_i = []
                for j in range(n_energy_out[i]):
                    mu = ace.xss[idx + 4:idx + 4 + n_mu]
                    p_mu = 1./n_mu*np.ones(n_mu)
                    mu_ij = Discrete(mu, p_mu)
                    mu_ij.c = np.cumsum(p_mu)
                    mu_i.append(mu_ij)
                    idx += 3 + n_mu

                energy_out.append(eout_i)
                mu_out.append(mu_i)

            # Create correlated angle-energy distribution
            breakpoints = [n_energy]
            interpolation = [2]
            energy = table.inelastic_xs.x
            table.inelastic_dist = CorrelatedAngleEnergy(
                breakpoints, interpolation, energy, energy_out, mu_out)

        # Incoherent/coherent elastic scattering cross section
        idx = ace.jxs[4]
        if idx != 0:
            n_energy = int(ace.xss[idx])
            energy = ace.xss[idx+1 : idx+1+n_energy]
            P = ace.xss[idx+1+n_energy : idx+1+2*n_energy]

            if ace.nxs[5] == 4:
                table.elastic_xs = CoherentElastic(energy, P)
            else:
                table.elastic_xs = Tabulated1D(energy, P)

            # Angular distribution
            n_mu = ace.nxs[6]
            if n_mu != -1:
                idx = ace.jxs[6]
                table.elastic_mu_out = ace.xss[idx:idx + n_energy*n_mu]
                table.elastic_mu_out.shape = (n_energy, n_mu)

        # Get relevant ZAIDs
        pairs = np.fromiter(map(lambda p: p[0], ace.pairs), int)
        table.zaids = pairs[np.nonzero(pairs)]

        return table
