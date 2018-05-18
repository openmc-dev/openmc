from collections import OrderedDict, Mapping, Callable
from copy import deepcopy
from io import StringIO
from numbers import Integral, Real
import os

import h5py
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline

from openmc.mixin import EqualityMixin
import openmc.checkvalue as cv
from . import HDF5_VERSION
from .ace import Table, get_metadata, get_table
from .data import ATOMIC_SYMBOL, EV_PER_MEV
from .endf import Evaluation, get_head_record, get_tab1_record, get_list_record
from .function import Tabulated1D


_SUBSHELLS = ['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5',
              'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2',
              'O3', 'O4', 'O5', 'O6', 'O7', 'O8', 'O9', 'P1', 'P2',
              'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'P11',
              'Q1', 'Q2', 'Q3']


# Helper function to map designator to subshell string or None
def _subshell(i):
    if i == 0:
        return None
    else:
        return _SUBSHELLS[i - 1]


_REACTION_NAME = {
    501: 'Total photon interaction',
    502: 'Photon coherent scattering',
    504: 'Photon incoherent scattering',
    515: 'Pair production, electron field',
    516: 'Total pair production',
    517: 'Pair production, nuclear field',
    522: 'Photoelectric absorption',
    526: 'Electro-atomic scattering',
    527: 'Electro-atomic bremsstrahlung',
    528: 'Electro-atomic excitation',
    534: 'K (1s1/2) subshell photoelectric',
    535: 'L1 (2s1/2) subshell photoelectric',
    536: 'L2 (2p1/2) subshell photoelectric',
    537: 'L3 (2p3/2) subshell photoelectric',
    538: 'M1 (3s1/2) subshell photoelectric',
    539: 'M2 (3p1/2) subshell photoelectric',
    540: 'M3 (3p3/2) subshell photoelectric',
    541: 'M4 (3d3/2) subshell photoelectric',
    542: 'M5 (3d5/2) subshell photoelectric',
    543: 'N1 (4s1/2) subshell photoelectric',
    544: 'N2 (4p1/2) subshell photoelectric',
    545: 'N3 (4p3/2) subshell photoelectric',
    546: 'N4 (4d3/2) subshell photoelectric',
    547: 'N5 (4d5/2) subshell photoelectric',
    548: 'N6 (4f5/2) subshell photoelectric',
    549: 'N7 (4f7/2) subshell photoelectric',
    550: 'O1 (5s1/2) subshell photoelectric',
    551: 'O2 (5p1/2) subshell photoelectric',
    552: 'O3 (5p3/2) subshell photoelectric',
    553: 'O4 (5d3/2) subshell photoelectric',
    554: 'O5 (5d5/2) subshell photoelectric',
    555: 'O6 (5f5/2) subshell photoelectric',
    556: 'O7 (5f7/2) subshell photoelectric',
    557: 'O8 (5g7/2) subshell photoelectric',
    558: 'O9 (5g9/2) subshell photoelectric',
    559: 'P1 (6s1/2) subshell photoelectric',
    560: 'P2 (6p1/2) subshell photoelectric',
    561: 'P3 (6p3/2) subshell photoelectric',
    562: 'P4 (6d3/2) subshell photoelectric',
    563: 'P5 (6d5/2) subshell photoelectric',
    564: 'P6 (6f5/2) subshell photoelectric',
    565: 'P7 (6f7/2) subshell photoelectric',
    566: 'P8 (6g7/2) subshell photoelectric',
    567: 'P9 (6g9/2) subshell photoelectric',
    568: 'P10 (6h9/2) subshell photoelectric',
    569: 'P11 (6h11/2) subshell photoelectric',
    570: 'Q1 (7s1/2) subshell photoelectric',
    571: 'Q2 (7p1/2) subshell photoelectric',
    572: 'Q3 (7p3/2) subshell photoelectric'
}

# Compton profiles are read from a pre-generated HDF5 file when they are first
# needed. The dictionary stores an array of electron momentum values (at which
# the profiles are tabulated) with the key 'pz' and the profile for each element
# is a 2D array with shape (n_shells, n_momentum_values) stored on the key Z
_COMPTON_PROFILES = {}

# Stopping powers are read from a pre-generated HDF5 file when they are first
# needed. The dictionary stores an array of energy values at which the other
# quantities are tabulated with the key 'energy' and for each element has the
# mass density, the mean excitation energy, and arrays containing the collision
# stopping powers, radiative stopping powers, and the density effect parameter
# stored on the key 'Z'.
_STOPPING_POWERS = {}

# Scaled bremsstrahlung DCSs are read from a data file provided by Selzter and
# Berger when they are first needed. The dictionary stores an array of n
# incident electron kinetic energies with key 'electron_energies', an array of
# k reduced photon energies with key 'photon_energies', and the cross sections
# for each element are in a 2D array with shape (n, k) stored on the key 'Z'.
_BREMSSTRAHLUNG = {}

# Reduced screening radii for Z = 1-99 from F. Salvat, J. M. Fernández-Varea,
# and J. Sempau, "PENELOPE-2011: A Code System for Monte Carlo Simulation of
# Electron and Photon Transport," OECD-NEA, Issy-les-Moulineaux, France (2011).
_REDUCED_SCREENING_RADIUS = [
    122.81, 73.167, 69.228, 67.301, 64.696, 61.228, 57.524, 54.033, 50.787,
    47.851, 46.373, 45.401, 44.503, 43.815, 43.074, 42.321, 41.586, 40.953,
    40.524, 40.256, 39.756, 39.144, 38.462, 37.778, 37.174, 36.663, 35.986,
    35.317, 34.688, 34.197, 33.786, 33.422, 33.068, 32.740, 32.438, 32.143,
    31.884, 31.622, 31.438, 31.142, 30.950, 30.758, 30.561, 30.285, 30.097,
    29.832, 29.581, 29.411, 29.247, 29.085, 28.930, 28.721, 28.580, 28.442,
    28.312, 28.139, 27.973, 27.819, 27.675, 27.496, 27.285, 27.093, 26.911,
    26.705, 26.516, 26.304, 26.108, 25.929, 25.730, 25.577, 25.403, 25.245,
    25.100, 24.941, 24.790, 24.655, 24.506, 24.391, 24.262, 24.145, 24.039,
    23.922, 23.813, 23.712, 23.621, 23.523, 23.430, 23.331, 23.238, 23.139,
    23.048, 22.967, 22.833, 22.694, 22.624, 22.545, 22.446, 22.358, 22.264
]

class AtomicRelaxation(EqualityMixin):
    """Atomic relaxation data.

    This class stores the binding energy, number of electrons, and electron
    transitions possible from ioniziation for each subshell with an atom. All of
    the data originates from an ENDF-6 atomic relaxation sub-library
    (NSUB=6). Instances of this class are not normally instantiated directly but
    rather created using the factory method :math:`AtomicRelaxation.from_endf`.

    Parameters
    ----------
    binding_energy : dict
        Dictionary indicating the binding energy in eV (values) for given
        subshells (keys). The subshells should be given as strings, e.g., 'K',
        'L1', 'L2', etc.
    num_electrons : dict
        Dictionary indicating the number of electrons in a subshell when neutral
        (values) for given subshells (keys). The subshells should be given as
        strings, e.g., 'K', 'L1', 'L2', etc.
    transitions : pandas.DataFrame
        Dictionary indicating allowed transitions and their probabilities
        (values) for given subshells (keys). The subshells should be given as
        strings, e.g., 'K', 'L1', 'L2', etc. The transitions are represented as
        a DataFrame with columns indicating the secondary and tertiary subshell,
        the energy of the transition in eV, and the fractional probability of
        the transition.

    Attributes
    ----------
    binding_energy : dict
        Dictionary indicating the binding energy in eV (values) for given
        subshells (keys). The subshells should be given as strings, e.g., 'K',
        'L1', 'L2', etc.
    num_electrons : dict
        Dictionary indicating the number of electrons in a subshell when neutral
        (values) for given subshells (keys). The subshells should be given as
        strings, e.g., 'K', 'L1', 'L2', etc.
    transitions : pandas.DataFrame
        Dictionary indicating allowed transitions and their probabilities
        (values) for given subshells (keys). The subshells should be given as
        strings, e.g., 'K', 'L1', 'L2', etc. The transitions are represented as
        a DataFrame with columns indicating the secondary and tertiary subshell,
        the energy of the transition in eV, and the fractional probability of
        the transition.

    See Also
    --------
    IncidentPhoton

    """
    def __init__(self, binding_energy, num_electrons, transitions):
        self.binding_energy = binding_energy
        self.num_electrons = num_electrons
        self.transitions = transitions
        self.compton_profile = OrderedDict()

    @property
    def binding_energy(self):
        return self._binding_energy

    @property
    def num_electrons(self):
        return self._num_electrons

    @property
    def subshells(self):
        return list(sorted(self.binding_energy.keys()))

    @property
    def transitions(self):
        return self._transitions

    @binding_energy.setter
    def binding_energy(self, binding_energy):
        cv.check_type('binding energies', binding_energy, Mapping)
        for subshell, energy in binding_energy.items():
            cv.check_value('subshell', subshell, _SUBSHELLS)
            cv.check_type('binding energy', energy, Real)
            cv.check_greater_than('binding energy', energy, 0.0, True)
        self._binding_energy = binding_energy

    @num_electrons.setter
    def num_electrons(self, num_electrons):
        cv.check_type('number of electrons', num_electrons, Mapping)
        for subshell, num in num_electrons.items():
            cv.check_value('subshell', subshell, _SUBSHELLS)
            cv.check_type('number of electrons', num, Real)
            cv.check_greater_than('number of electrons', num, 0.0, True)
        self._num_electrons = num_electrons

    @transitions.setter
    def transitions(self, transitions):
        cv.check_type('transitions', transitions, Mapping)
        for subshell, df in transitions.items():
            cv.check_value('subshell', subshell, _SUBSHELLS)
            cv.check_type('transitions', df, pd.DataFrame)
        self._transitions = transitions

    @classmethod
    def from_ace(cls, ace):
        """Generate atomic relaxation data from an ACE file

        Parameters
        ----------
        ace : openmc.data.ace.Table
            ACE table to read from

        Returns
        -------
        openmc.data.AtomicRelaxation
            Atomic relaxation data

        """
        # Create data dictionaries
        binding_energy = {}
        num_electrons = {}
        transitions = {}

        # Get shell designators
        n = ace.nxs[7]
        idx = ace.jxs[11]
        shells = [_subshell(int(i)) for i in ace.xss[idx : idx+n]]

        # Get number of electrons for each shell
        idx = ace.jxs[12]
        for shell, num in zip(shells, ace.xss[idx : idx+n]):
            num_electrons[shell] = num

        # Get binding energy for each shell
        idx = ace.jxs[13]
        for shell, e in zip(shells, ace.xss[idx : idx+n]):
            binding_energy[shell] = e*EV_PER_MEV

        # Get transition table
        columns = ['secondary', 'tertiary', 'energy (eV)', 'probability']
        idx = ace.jxs[18]
        for i, subi in enumerate(shells):
            n_transitions = int(ace.xss[ace.jxs[15] + i])
            if n_transitions > 0:
                records = []
                for j in range(n_transitions):
                    subj = _subshell(int(ace.xss[idx]))
                    subk = _subshell(int(ace.xss[idx + 1]))
                    etr = ace.xss[idx + 2]*EV_PER_MEV
                    if j == 0:
                        ftr = ace.xss[idx + 3]
                    else:
                        ftr = ace.xss[idx + 3] - ace.xss[idx - 1]
                    records.append((subj, subk, etr, ftr))
                    idx += 4

                # Create dataframe for transitions
                transitions[subi] = pd.DataFrame.from_records(
                    records, columns=columns)

        return cls(binding_energy, num_electrons, transitions)

    @classmethod
    def from_endf(cls, ev_or_filename):
        """Generate atomic relaxation data from an ENDF evaluation

        Parameters
        ----------
        ev_or_filename : str or openmc.data.endf.Evaluation
            ENDF atomic relaxation evaluation to read from. If given as a
            string, it is assumed to be the filename for the ENDF file.

        Returns
        -------
        openmc.data.AtomicRelaxation
            Atomic relaxation data

        """
        if isinstance(ev_or_filename, Evaluation):
            ev = ev_or_filename
        else:
            ev = Evaluation(ev_or_filename)

        # Atomic relaxation data is always MF=28, MT=533
        if (28, 533) not in ev.section:
            raise IOError('{} does not appear to be an atomic relaxation '
                          'sublibrary.'.format(ev))

        # Determine number of subshells
        file_obj = StringIO(ev.section[28, 533])
        params = get_head_record(file_obj)
        n_subshells = params[4]

        # Create data dictionaries
        binding_energy = {}
        num_electrons = {}
        transitions = {}
        columns = ['secondary', 'tertiary', 'energy (eV)', 'probability']

        # Read data for each subshell
        for i in range(n_subshells):
            params, list_items = get_list_record(file_obj)
            subi = _subshell(int(params[0]))
            n_transitions = int(params[5])
            binding_energy[subi] = list_items[0]
            num_electrons[subi] = list_items[1]

            if n_transitions > 0:
                # Read transition data
                records = []
                for j in range(n_transitions):
                    subj = _subshell(int(list_items[6*(j+1)]))
                    subk = _subshell(int(list_items[6*(j+1) + 1]))
                    etr = list_items[6*(j+1) + 2]
                    ftr = list_items[6*(j+1) + 3]
                    records.append((subj, subk, etr, ftr))

                # Create dataframe for transitions
                transitions[subi] = pd.DataFrame.from_records(
                    records, columns=columns)

        # Return instance of class
        return cls(binding_energy, num_electrons, transitions)

    def to_hdf5(self, group):
        raise NotImplementedError


class IncidentPhoton(EqualityMixin):
    """Photon interaction data.

    This class stores photo-atomic, photo-nuclear, atomic relaxation,
    Compton profile, stopping power, and bremsstrahlung data assembled from
    different sources. To create an instance, the factory method
    :meth:`IncidentPhoton.from_endf` can be used. To add atomic relaxation or
    Compton profile data, set the :attr:`IncidentPhoton.atomic_relaxation` and
    :attr:`IncidentPhoton.compton_profiles` attributes directly.

    Parameters
    ----------
    atomic_number : int
        Number of protons in the target nucleus

    Attributes
    ----------
    atomic_number : int
        Number of protons in the target nucleus
    atomic_relaxation : openmc.data.AtomicRelaxation or None
        Atomic relaxation data
    bremsstrahlung : dict
        Dictionary of bremsstrahlung DCS data with keys 'electron_energy'
        (incident electron kinetic energy values in eV), 'photon_energy'
        (ratio of the energy of the emitted photon to the incident electron
        kinetic energy), and 'dcs' (cross sectin values in mb). The cross
        sections are in scaled form: :math:`(\beta^2/Z^2) E_k (d\sigma/dE_k)`,
        where :math:`E_k` is the energy of the emitted photon.
    compton_profiles : dict
        Dictionary of Compton profile data with keys 'num_electrons' (number of
        electrons in each subshell), 'binding_energy' (ionization potential of
        each subshell), and 'J' (Hartree-Fock Compton profile as a function of
        the projection of the electron momentum on the scattering vector,
        :math:`p_z` for each subshell). Note that subshell occupancies may not
        match the atomic relaxation data.
    reactions : collections.OrderedDict
        Contains the cross sections for each photon reaction. The keys are MT
        values and the values are instances of :class:`PhotonReaction`.
    reduced_screening_radius : float
        Reduced screening radius :math:`R m_e c/\hbar`, where R is the screening
        radius for an atom of atomic number Z under the assumption that the
        Coulomb field of the nucleus is exponentially screened by atomic electrons.
        :math:`\hbar/m_e c` is the Compton wavelength of the electron.
    stopping_powers : dict
        Dictionary of stopping power data with keys 'energy' (in eV), 'density'
        (mass density in g/cm:sup:`3`), 'I' (mean excitation energy),
        's_collision' (collision stopping power in eV cm:sup:`2`/g),
        's_radiative' (radiative stopping power in eV cm:sup:`2`/g), and
        'density_effect' (density effect parameter).
    summed_reactions : collections.OrderedDict
        Contains summed cross sections. The keys are MT values and the values
        are instances of :class:`PhotonReaction`.

    """

    def __init__(self, atomic_number):
        self.atomic_number = atomic_number
        self._atomic_relaxation = None
        self.reactions = OrderedDict()
        self.summed_reactions = OrderedDict()
        self.compton_profiles = {}
        self.stopping_powers = {}
        self.bremsstrahlung = {}

    def __contains__(self, mt):
        return mt in self.reactions or mt in self.summed_reactions

    def __getitem__(self, mt):
        if mt in self.reactions:
            return self.reactions[mt]
        elif mt in self.summed_reactions:
            return self.summed_reactions[mt]
        else:
            raise KeyError('No reaction with MT={}.'.format(mt))

    def __repr__(self):
        return "<IncidentPhoton: {}>".format(self.name)

    def __iter__(self):
        return iter(self.reactions.values())

    @property
    def atomic_number(self):
        return self._atomic_number

    @property
    def atomic_relaxation(self):
        return self._atomic_relaxation

    @property
    def name(self):
        return ATOMIC_SYMBOL[self.atomic_number]

    @property
    def reduced_screening_radius(self):
        if self.atomic_number < 100:
            return _REDUCED_SCREENING_RADIUS[self.atomic_number - 1]
        else:
            raise IndexError('No reduced screening radius for '
                             'Z={}.'.format(self.atomic_number))

    @atomic_number.setter
    def atomic_number(self, atomic_number):
        cv.check_type('atomic number', atomic_number, Integral)
        cv.check_greater_than('atomic number', atomic_number, 0, True)
        self._atomic_number = atomic_number

    @atomic_relaxation.setter
    def atomic_relaxation(self, atomic_relaxation):
        cv.check_type('atomic relaxation data', atomic_relaxation,
                      AtomicRelaxation)
        self._atomic_relaxation = atomic_relaxation

    @classmethod
    def from_ace(cls, ace_or_filename):
        """Generate incident photon data from an ACE table

        Parameters
        ----------
        ace_or_filename : str or openmc.data.ace.Table
            ACE table to read from. If given as a string, it is assumed to be
            the filename for the ACE file.

        Returns
        -------
        openmc.data.IncidentPhoton
            Photon interaction data

        """
        # First obtain the data for the first provided ACE table/file
        if isinstance(ace_or_filename, Table):
            ace = ace_or_filename
        else:
            ace = get_table(ace_or_filename)

        # Get atomic number based on name of ACE table
        zaid = ace.name.split('.')[0]
        Z = get_metadata(int(zaid))[2]

        # Read each reaction
        data = cls(Z)
        for mt in (502, 504, 515, 522):
            data.reactions[mt] = PhotonReaction.from_ace(ace, mt)

        # Compton profiles
        n_shell = ace.nxs[5]
        if n_shell != 0:
            # Get number of electrons in each shell
            idx = ace.jxs[6]
            data.compton_profiles['num_electrons'] = ace.xss[idx : idx+n_shell]

            # Get binding energy for each shell
            idx = ace.jxs[7]
            data.compton_profiles['binding_energy'] = ace.xss[idx : idx+n_shell]

            # Create Compton profile for each electron shell
            profiles = []
            for k in range(n_shell):
                # Get number of momentum values and interpolation scheme
                loca = int(ace.xss[ace.jxs[9] + k])
                jj = int(ace.xss[ace.jxs[10] + loca - 1])
                m = int(ace.xss[ace.jxs[10] + loca])

                # Read momentum and PDF
                idx = ace.jxs[10] + loca + 1
                pz = ace.xss[idx : idx+m]
                pdf = ace.xss[idx+m : idx+2*m]

                # Create proflie function
                J_k = Tabulated1D(pz, pdf, [m], [jj])
                profiles.append(J_k)
            data.compton_profiles['J'] = profiles

        # Subshell photoelectric xs and atomic relaxation data
        if ace.nxs[7] > 0:
            data.atomic_relaxation = AtomicRelaxation.from_ace(ace)

            # Get subshell designators
            n_subshells = ace.nxs[7]
            idx = ace.jxs[11]
            designators = [int(i) for i in ace.xss[idx : idx+n_subshells]]

            # Get energy grid for subshell photoionization
            n_energy = ace.nxs[3]
            idx = ace.jxs[1]
            energy = np.exp(ace.xss[idx : idx+n_energy])*EV_PER_MEV

            # Get cross section for each subshell
            idx = ace.jxs[16]
            for d in designators:
                # Create photon reaction
                mt = 533 + d
                rx = PhotonReaction(mt)
                data.reactions[mt] = rx

                # Store cross section
                xs = ace.xss[idx : idx+n_energy].copy()
                nonzero = (xs != 0.0)
                xs[nonzero] = np.exp(xs[nonzero])
                rx.xs = Tabulated1D(energy, xs, [n_energy], [5])
                idx += n_energy

                # Copy binding energy
                shell = _subshell(d)
                e = data.atomic_relaxation.binding_energy[shell]
                rx.subshell_binding_energy = e

        return data

    @classmethod
    def from_endf(cls, photoatomic, relaxation=None):
        """Generate incident photon data from an ENDF evaluation

        Parameters
        ----------
        photoatomic : str or openmc.data.endf.Evaluation
            ENDF photoatomic data evaluation to read from. If given as a string,
            it is assumed to be the filename for the ENDF file.
        relaxation : str or openmc.data.endf.Evaluation, optional
            ENDF atomic relaxation data evaluation to read from. If given as a
            string, it is assumed to be the filename for the ENDF file.

        Returns
        -------
        openmc.data.IncidentPhoton
            Photon interaction data

        """
        if isinstance(photoatomic, Evaluation):
            ev = photoatomic
        else:
            ev = Evaluation(photoatomic)

        Z = ev.target['atomic_number']
        data = cls(Z)

        # Read each reaction
        for mf, mt, nc, mod in ev.reaction_list:
            if mf == 23:
                data.reactions[mt] = PhotonReaction.from_endf(ev, mt)

        # Add atomic relaxation data if it hasn't been added already
        if relaxation is not None:
            data.atomic_relaxation = AtomicRelaxation.from_endf(relaxation)

        # If Compton profile data hasn't been loaded, do so
        if not _COMPTON_PROFILES:
            filename = os.path.join(os.path.dirname(__file__), 'compton_profiles.h5')
            with h5py.File(filename, 'r') as f:
                _COMPTON_PROFILES['pz'] = f['pz'].value
                for i in range(1, 101):
                    group = f['{:03}'.format(i)]
                    num_electrons = group['num_electrons'].value
                    binding_energy = group['binding_energy'].value*EV_PER_MEV
                    J = group['J'].value
                    _COMPTON_PROFILES[i] = {'num_electrons': num_electrons,
                                            'binding_energy': binding_energy,
                                            'J': J}

        # Add Compton profile data
        pz = _COMPTON_PROFILES['pz']
        profile = _COMPTON_PROFILES[Z]
        data.compton_profiles['num_electrons'] = profile['num_electrons']
        data.compton_profiles['binding_energy'] = profile['binding_energy']
        data.compton_profiles['J'] = [Tabulated1D(pz, J_k) for J_k in profile['J']]

        # Load stopping power data if it has not yet been loaded
        if not _STOPPING_POWERS:
            filename = os.path.join(os.path.dirname(__file__), 'stopping_powers.h5')
            with h5py.File(filename, 'r') as f:
                # Units are in MeV; convert to eV
                _STOPPING_POWERS['energy'] = f['energy'].value*EV_PER_MEV
                for i in range(1, 99):
                    group = f['{:03}'.format(i)]
                    _STOPPING_POWERS[i] = {'density': group.attrs['density'],
                                           'I': group.attrs['I'],
                                           's_collision': group['s_collision'].value,
                                           's_radiative': group['s_radiative'].value,
                                           'density_effect': group['density_effect'].value}

                    # Units are in MeV cm^2/g; convert to eV cm^2/g
                    _STOPPING_POWERS[i]['s_collision'] *= EV_PER_MEV
                    _STOPPING_POWERS[i]['s_radiative'] *= EV_PER_MEV

        # Add stopping power data
        if Z < 99:
            data.stopping_powers['energy'] = _STOPPING_POWERS['energy']
            data.stopping_powers.update(_STOPPING_POWERS[Z])

        # Load bremsstrahlung data if it has not yet been loaded
        if not _BREMSSTRAHLUNG:
            filename = os.path.join(os.path.dirname(__file__), 'BREMX.DAT')
            brem = open(filename, 'r').read().split()

            # Incident electron kinetic energy grid in eV
            _BREMSSTRAHLUNG['electron_energy'] = np.logspace(3, 9, 200)
            log_energy = np.log(_BREMSSTRAHLUNG['electron_energy'])

            # Get number of tabulated electron and photon energy values
            n = int(brem[37])
            k = int(brem[38])

            # Index in data
            p = 39

            # Get log of incident electron kinetic energy values, used for cubic
            # spline interpolation in log energy. Units are in MeV, so convert to eV.
            logx = np.log(np.fromiter(brem[p:p+n], float, n)*EV_PER_MEV)
            p += n

            # Get reduced photon energy values
            _BREMSSTRAHLUNG['photon_energy'] = np.fromiter(brem[p:p+k], float, k)
            p += k

            for i in range(1, 101):
                dcs = np.empty([len(log_energy), k])

                # Get the scaled cross section values for each electron energy and
                # reduced photon energy for this Z
                y = np.reshape(np.fromiter(brem[p:p+n*k], float, n*k), (n, k))
                p += k*n

                for j in range(k):
                    # Cubic spline interpolation in log energy and linear DCS
                    cs = CubicSpline(logx, y[:,j])

                    # Get scaled DCS values (millibarns) on new energy grid
                    dcs[:,j] = cs(log_energy)

                _BREMSSTRAHLUNG[i] = {'dcs': dcs}

        # Add bremsstrahlung DCS data
        data.bremsstrahlung['electron_energy'] = _BREMSSTRAHLUNG['electron_energy']
        data.bremsstrahlung['photon_energy'] = _BREMSSTRAHLUNG['photon_energy']
        data.bremsstrahlung['dcs'] = _BREMSSTRAHLUNG[Z]['dcs']

        return data

    def export_to_hdf5(self, path, mode='a'):
        """Export incident photon data to an HDF5 file.

        Parameters
        ----------
        path : str
            Path to write HDF5 file to
        mode : {'r', r+', 'w', 'x', 'a'}
            Mode that is used to open the HDF5 file. This is the second argument
            to the :class:`h5py.File` constructor.

        """
        # Open file and write version
        f = h5py.File(path, mode, libver='latest')
        f.attrs['filetype'] = np.string_('data_photon')
        if 'version' not in f.attrs:
            f.attrs['version'] = np.array(HDF5_VERSION)

        group = f.create_group(self.name)
        group.attrs['Z'] = Z = self.atomic_number

        # Determine union energy grid
        union_grid = np.array([])
        for rx in self:
            union_grid = np.union1d(union_grid, rx.xs.x)
        group.create_dataset('energy', data=union_grid)

        # Write coherent scattering cross section
        rx = self.reactions[502]
        coh_group = group.create_group('coherent')
        coh_group.create_dataset('xs', data=rx.xs(union_grid))
        if rx.scattering_factor is not None:
            # Create integrated form factor
            ff = deepcopy(rx.scattering_factor)
            ff.x *= ff.x
            ff.y *= ff.y/Z**2
            int_ff = Tabulated1D(ff.x, ff.integral())
            int_ff.to_hdf5(coh_group, 'integrated_scattering_factor')
        if rx.anomalous_real is not None:
            rx.anomalous_real.to_hdf5(coh_group, 'anomalous_real')
        if rx.anomalous_imag is not None:
            rx.anomalous_imag.to_hdf5(coh_group, 'anomalous_imag')

        # Write incoherent scattering cross section
        rx = self[504]
        incoh_group = group.create_group('incoherent')
        incoh_group.create_dataset('xs', data=rx.xs(union_grid))
        if rx.scattering_factor is not None:
            rx.scattering_factor.to_hdf5(incoh_group, 'scattering_factor')

        # Write electron-field pair production cross section
        if 515 in self:
            pair_group = group.create_group('pair_production_electron')
            pair_group.create_dataset('xs', data=self[515].xs(union_grid))

        # Write nuclear-field pair production cross section
        if 517 in self:
            pair_group = group.create_group('pair_production_nuclear')
            pair_group.create_dataset('xs', data=self[517].xs(union_grid))

        # Write photoelectric cross section
        photoelec_group = group.create_group('photoelectric')
        photoelec_group.create_dataset('xs', data=self[522].xs(union_grid))

        # Write photoionization cross sections
        shell_group = group.create_group('subshells')
        designators = []
        for mt, rx in self.reactions.items():
            if mt >= 534 and mt <= 572:
                # Get name of subshell
                shell = _SUBSHELLS[mt - 534]
                designators.append(shell)
                sub_group = shell_group.create_group(shell)

                if self.atomic_relaxation is not None:
                    relax = self.atomic_relaxation
                    # Write subshell binding energy and number of electrons
                    sub_group.attrs['binding_energy'] = relax.binding_energy[shell]
                    sub_group.attrs['num_electrons'] = relax.num_electrons[shell]

                    # Write transition data with replacements
                    if shell in relax.transitions:
                        shell_values = _SUBSHELLS.copy()
                        shell_values.insert(0, None)
                        df = relax.transitions[shell].replace(
                            shell_values, range(len(shell_values)))
                        sub_group.create_dataset('transitions', data=df.as_matrix())

                # Determine threshold
                threshold = rx.xs.x[0]
                idx = np.searchsorted(union_grid, threshold, side='right') - 1

                # Interpolate cross section onto union grid and write
                photoionization = rx.xs(union_grid[idx:])
                sub_group.create_dataset('xs', data=photoionization)
                assert len(union_grid) == len(photoionization) + idx
                sub_group['xs'].attrs['threshold_idx'] = idx

        shell_group.attrs['designators'] = np.array(designators, dtype='S')

        # Write reduced screening radius
        if Z < 100:
            group.attrs['reduced_screening_radius'] = self.reduced_screening_radius

        # Write Compton profiles
        if self.compton_profiles:
            compton_group = group.create_group('compton_profiles')

            profile = self.compton_profiles
            compton_group.create_dataset('num_electrons',
                                         data=profile['num_electrons'])
            compton_group.create_dataset('binding_energy',
                                         data=profile['binding_energy'])

            # Get electron momentum values
            compton_group.create_dataset('pz', data=profile['J'][0].x)

            # Create/write 2D array of profiles
            J = np.array([Jk.y for Jk in profile['J']])
            compton_group.create_dataset('J', data=J)

        # Write stopping powers
        if self.stopping_powers:
            s_group = group.create_group('stopping_powers')

            for key, value in self.stopping_powers.items():
                if key in ('density', 'I'):
                    s_group.attrs[key] = value
                else:
                    s_group.create_dataset(key, data=value)

        # Write bremsstrahlung
        if self.bremsstrahlung:
            brem_group = group.create_group('bremsstrahlung')

            brem = self.bremsstrahlung
            brem_group.create_dataset('electron_energy',
                                      data=brem['electron_energy'])
            brem_group.create_dataset('photon_energy',
                                      data=brem['photon_energy'])
            brem_group.create_dataset('dcs', data=brem['dcs'])


class PhotonReaction(EqualityMixin):
    """Photon-induced reaction

    Parameters
    ----------
    mt : int
        The ENDF MT number for this reaction.

    Attributes
    ----------
    anomalous_real : openmc.data.Tabulated1D
        Real part of the anomalous scattering factor
    anomlaous_imag : openmc.data.Tabulated1D
        Imaginary part of the anomalous scatttering factor
    mt : int
        The ENDF MT number for this reaction.
    scattering_factor : openmc.data.Tabulated1D
        Coherent or incoherent form factor.
    xs : Callable
        Cross section as a function of incident photon energy

    """

    def __init__(self, mt):
        self.mt = mt
        self._xs = None
        self._scattering_factor = None
        self._anomalous_real = None
        self._anomalous_imag = None

    def __repr__(self):
        if self.mt in _REACTION_NAME:
            return "<Photon Reaction: MT={} {}>".format(
                self.mt, _REACTION_NAME[self.mt])
        else:
            return "<Photon Reaction: MT={}>".format(self.mt)

    @property
    def anomalous_real(self):
        return self._anomalous_real

    @property
    def anomalous_imag(self):
        return self._anomalous_imag

    @property
    def scattering_factor(self):
        return self._scattering_factor

    @property
    def xs(self):
        return self._xs

    @anomalous_real.setter
    def anomalous_real(self, anomalous_real):
        cv.check_type('real part of anomalous scattering factor',
                      anomalous_real, Callable)
        self._anomalous_real = anomalous_real

    @anomalous_imag.setter
    def anomalous_imag(self, anomalous_imag):
        cv.check_type('imaginary part of anomalous scattering factor',
                      anomalous_imag, Callable)
        self._anomalous_imag = anomalous_imag

    @scattering_factor.setter
    def scattering_factor(self, scattering_factor):
        cv.check_type('scattering factor', scattering_factor, Callable)
        self._scattering_factor = scattering_factor

    @xs.setter
    def xs(self, xs):
        cv.check_type('reaction cross section', xs, Callable)
        self._xs = xs

    @classmethod
    def from_ace(cls, ace, mt):
        """Generate photon reaction from an ACE table

        Parameters
        ----------
        ace : openmc.data.ace.Table
            ACE table to read from
        mt : int
            The MT value of the reaction to get data for

        Returns
        -------
        openmc.data.PhotonReaction
            Photon reaction data

        """
        # Create instance
        rx = cls(mt)

        # Get energy grid (stored as logarithms)
        n = ace.nxs[3]
        idx = ace.jxs[1]
        energy = np.exp(ace.xss[idx : idx+n])*EV_PER_MEV

        # Get index for appropriate reaction
        if mt == 502:
            # Coherent scattering
            idx = ace.jxs[1] + 2*n
        elif mt == 504:
            # Incoherent scattering
            idx = ace.jxs[1] + n
        elif mt == 515:
            # Pair production
            idx = ace.jxs[1] + 4*n
        elif mt == 522:
            # Photoelectric
            idx = ace.jxs[1] + 3*n
        else:
            raise ValueError('ACE photoatomic cross sections do not have '
                             'data for MT={}.'.format(mt))

        # Store cross section
        xs = ace.xss[idx : idx+n].copy()
        nonzero = (xs != 0.0)
        xs[nonzero] = np.exp(xs[nonzero])
        rx.xs = Tabulated1D(energy, xs, [n], [5])

        # Get form factors for incoherent/coherent scattering
        new_format = (ace.nxs[6] > 0)
        if mt == 502:
            idx = ace.jxs[3]
            if new_format:
                n = (ace.jxs[4] - ace.jxs[3]) // 3
                x = ace.xss[idx : idx+n]
                idx += n
            else:
                x = np.array([
                    0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.12,
                    0.15, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
                    0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
                    1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4,
                    3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6,
                    5.8, 6.0])
                n = x.size
            ff = ace.xss[idx+n : idx+2*n]
            rx.scattering_factor = Tabulated1D(x, ff)

        elif mt == 504:
            idx = ace.jxs[2]
            if new_format:
                n = (ace.jxs[3] - ace.jxs[2]) // 2
                x = ace.xss[idx : idx+n]
                idx += n
            else:
                x = np.array([
                    0.0, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6,
                    0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 8.0
                ])
                n = x.size
            ff = ace.xss[idx : idx+n]
            rx.scattering_factor = Tabulated1D(x, ff)

        return rx

    @classmethod
    def from_endf(cls, ev, mt):
        """Generate photon reaction from an ENDF evaluation

        Parameters
        ----------
        ev : openmc.data.endf.Evaluation
            ENDF photo-atomic interaction data evaluation
        mt : int
            The MT value of the reaction to get data for

        Returns
        -------
        openmc.data.PhotonReaction
            Photon reaction data

        """
        rx = cls(mt)

        # Read photon cross section
        if (23, mt) in ev.section:
            file_obj = StringIO(ev.section[23, mt])
            get_head_record(file_obj)
            params, rx.xs = get_tab1_record(file_obj)

            # Set subshell binding energy and/or fluorescence yield
            if mt >= 534 and mt <= 599:
                rx.subshell_binding_energy = params[0]
            if mt >= 534 and mt <= 572:
                rx.fluorescence_yield = params[1]

        # Read form factors / scattering functions
        if (27, mt) in ev.section:
            file_obj = StringIO(ev.section[27, mt])
            get_head_record(file_obj)
            params, rx.scattering_factor = get_tab1_record(file_obj)

        # Check for anomalous scattering factor
        if mt == 502:
            if (27, 506) in ev.section:
                file_obj = StringIO(ev.section[27, 506])
                get_head_record(file_obj)
                params, rx.anomalous_real = get_tab1_record(file_obj)

            if (27, 505) in ev.section:
                file_obj = StringIO(ev.section[27, 505])
                get_head_record(file_obj)
                params, rx.anomalous_imag = get_tab1_record(file_obj)

        return rx
