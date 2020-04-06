from collections import defaultdict
from collections.abc import MutableSequence, Iterable
import io

import numpy as np
from numpy.polynomial import Polynomial
import pandas as pd

from .data import NEUTRON_MASS
from .endf import get_head_record, get_cont_record, get_tab1_record, get_list_record
try:
    from .reconstruct import wave_number, penetration_shift, reconstruct_mlbw, \
        reconstruct_slbw, reconstruct_rm
    _reconstruct = True
except ImportError:
    _reconstruct = False
import openmc.checkvalue as cv


class Resonances:
    """Resolved and unresolved resonance data

    Parameters
    ----------
    ranges : list of openmc.data.ResonanceRange
        Distinct energy ranges for resonance data

    Attributes
    ----------
    ranges : list of openmc.data.ResonanceRange
        Distinct energy ranges for resonance data
    resolved : openmc.data.ResonanceRange or None
        Resolved resonance range
    unresolved : openmc.data.Unresolved or None
        Unresolved resonance range

    """

    def __init__(self, ranges):
        self.ranges = ranges

    def __iter__(self):
        for r in self.ranges:
            yield r

    @property
    def ranges(self):
        return self._ranges

    @property
    def resolved(self):
        resolved_ranges = [r for r in self.ranges
                           if not isinstance(r, Unresolved)]
        if len(resolved_ranges) > 1:
            raise ValueError('More than one resolved range present')
        elif len(resolved_ranges) == 0:
            return None
        else:
            return resolved_ranges[0]

    @property
    def unresolved(self):
        for r in self.ranges:
            if isinstance(r, Unresolved):
                return r
        else:
            return None

    @ranges.setter
    def ranges(self, ranges):
        cv.check_type('resonance ranges', ranges, MutableSequence)
        self._ranges = cv.CheckedList(ResonanceRange, 'resonance ranges',
                                      ranges)

    @classmethod
    def from_endf(cls, ev):
        """Generate resonance data from an ENDF evaluation.

        Parameters
        ----------
        ev : openmc.data.endf.Evaluation
            ENDF evaluation

        Returns
        -------
        openmc.data.Resonances
            Resonance data

        """
        file_obj = io.StringIO(ev.section[2, 151])

        # Determine whether discrete or continuous representation
        items = get_head_record(file_obj)
        n_isotope = items[4]  # Number of isotopes

        ranges = []
        for iso in range(n_isotope):
            items = get_cont_record(file_obj)
            abundance = items[1]
            fission_widths = (items[3] == 1)  # fission widths are given?
            n_ranges = items[4]  # number of resonance energy ranges

            for j in range(n_ranges):
                items = get_cont_record(file_obj)
                resonance_flag = items[2]  # flag for resolved (1)/unresolved (2)
                formalism = items[3]  # resonance formalism

                if resonance_flag in (0, 1):
                    # resolved resonance region
                    erange = _FORMALISMS[formalism].from_endf(ev, file_obj, items)

                elif resonance_flag == 2:
                    # unresolved resonance region
                    erange = Unresolved.from_endf(file_obj, items, fission_widths)

                # erange.material = self
                ranges.append(erange)

        return cls(ranges)


class ResonanceRange:
    """Resolved resonance range

    Parameters
    ----------
    target_spin : float
        Intrinsic spin, :math:`I`, of the target nuclide
    energy_min : float
        Minimum energy of the resolved resonance range in eV
    energy_max : float
        Maximum energy of the resolved resonance range in eV
    channel : dict
        Dictionary whose keys are l-values and values are channel radii as a
        function of energy
    scattering : dict
        Dictionary whose keys are l-values and values are scattering radii as a
        function of energy

    Attributes
    ----------
    channel_radius : dict
        Dictionary whose keys are l-values and values are channel radii as a
        function of energy
    energy_max : float
        Maximum energy of the resolved resonance range in eV
    energy_min : float
        Minimum energy of the resolved resonance range in eV
    scattering_radius : dict
        Dictionary whose keys are l-values and values are scattering radii as a
        function of energ
    target_spin : float
        Intrinsic spin, :math:`I`, of the target nuclide

    """
    def __init__(self, target_spin, energy_min, energy_max, channel, scattering):
        self.target_spin = target_spin
        self.energy_min = energy_min
        self.energy_max = energy_max
        self.channel_radius = channel
        self.scattering_radius = scattering

        self._prepared = False
        self._parameter_matrix = {}

    def __copy__(self):
        cls = type(self)
        new_copy = cls.__new__(cls)
        new_copy.__dict__.update(self.__dict__)
        new_copy._prepared = False
        return new_copy

    @classmethod
    def from_endf(cls, ev, file_obj, items):
        """Create resonance range from an ENDF evaluation.

        This factory method is only used when LRU=0, indicating that only a
        scattering radius appears in MF=2, MT=151. All subclasses of
        ResonanceRange override this method with their own.

        Parameters
        ----------
        ev : openmc.data.endf.Evaluation
            ENDF evaluation
        file_obj : file-like object
            ENDF file positioned at the second record of a resonance range
            subsection in MF=2, MT=151
        items : list
            Items from the CONT record at the start of the resonance range
            subsection

        Returns
        -------
        openmc.data.ResonanceRange
            Resonance range data

        """
        energy_min, energy_max = items[0:2]

        # For scattering radius-only, NRO must be zero
        assert items[4] == 0

        # Get energy-independent scattering radius
        items = get_cont_record(file_obj)
        target_spin = items[0]
        ap = Polynomial((items[1],))

        # Calculate channel radius from ENDF-102 equation D.14
        a = Polynomial((0.123 * (NEUTRON_MASS*ev.target['mass'])**(1./3.) + 0.08,))

        return cls(target_spin, energy_min, energy_max, {0: a}, {0: ap})

    def reconstruct(self, energies):
        """Evaluate cross section at specified energies.

        Parameters
        ----------
        energies : float or Iterable of float
            Energies at which the cross section should be evaluated

        Returns
        -------
        3-tuple of float or numpy.ndarray
            Elastic, capture, and fission cross sections at the specified
            energies

        """
        if not _reconstruct:
            raise RuntimeError("Resonance reconstruction not available.")

        # Pre-calculate penetrations and shifts for resonances
        if not self._prepared:
            self._prepare_resonances()

        if isinstance(energies, Iterable):
            elastic = np.zeros_like(energies)
            capture = np.zeros_like(energies)
            fission = np.zeros_like(energies)

            for i, E in enumerate(energies):
                xse, xsg, xsf = self._reconstruct(self, E)
                elastic[i] = xse
                capture[i] = xsg
                fission[i] = xsf
        else:
            elastic, capture, fission = self._reconstruct(self, energies)

        return {2: elastic, 102: capture, 18: fission}


class MultiLevelBreitWigner(ResonanceRange):
    """Multi-level Breit-Wigner resolved resonance formalism data.

    Multi-level Breit-Wigner resolved resonance data is identified by LRF=2 in
    the ENDF-6 format.

    Parameters
    ----------
    target_spin : float
        Intrinsic spin, :math:`I`, of the target nuclide
    energy_min : float
        Minimum energy of the resolved resonance range in eV
    energy_max : float
        Maximum energy of the resolved resonance range in eV
    channel : dict
        Dictionary whose keys are l-values and values are channel radii as a
        function of energy
    scattering : dict
        Dictionary whose keys are l-values and values are scattering radii as a
        function of energy

    Attributes
    ----------
    atomic_weight_ratio : float
        Atomic weight ratio of the target nuclide given as a function of
        l-value. Note that this may be different than the value for the
        evaluation as a whole.
    channel_radius : dict
        Dictionary whose keys are l-values and values are channel radii as a
        function of energy
    energy_max : float
        Maximum energy of the resolved resonance range in eV
    energy_min : float
        Minimum energy of the resolved resonance range in eV
    parameters : pandas.DataFrame
        Energies, spins, and resonances widths for each resonance
    q_value : dict
        Q-value to be added to incident particle's center-of-mass energy to
        determine the channel energy for use in the penetrability factor. The
        keys of the dictionary are l-values.
    scattering_radius : dict
        Dictionary whose keys are l-values and values are scattering radii as a
        function of energy
    target_spin : float
        Intrinsic spin, :math:`I`, of the target nuclide

    """

    def __init__(self, target_spin, energy_min, energy_max, channel, scattering):
        super().__init__(target_spin, energy_min, energy_max, channel,
                         scattering)
        self.parameters = None
        self.q_value = {}
        self.atomic_weight_ratio = None

        # Set resonance reconstruction function
        if _reconstruct:
            self._reconstruct = reconstruct_mlbw
        else:
            self._reconstruct = None

    @classmethod
    def from_endf(cls, ev, file_obj, items):
        """Create MLBW data from an ENDF evaluation.

        Parameters
        ----------
        ev : openmc.data.endf.Evaluation
            ENDF evaluation
        file_obj : file-like object
            ENDF file positioned at the second record of a resonance range
            subsection in MF=2, MT=151
        items : list
            Items from the CONT record at the start of the resonance range
            subsection

        Returns
        -------
        openmc.data.MultiLevelBreitWigner
            Multi-level Breit-Wigner resonance parameters

        """

        # Read energy-dependent scattering radius if present
        energy_min, energy_max = items[0:2]
        nro, naps = items[4:6]
        if nro != 0:
            params, ape = get_tab1_record(file_obj)

        # Other scatter radius parameters
        items = get_cont_record(file_obj)
        target_spin = items[0]
        ap = Polynomial((items[1],))  # energy-independent scattering-radius
        NLS = items[4]  # number of l-values

        # Read resonance widths, J values, etc
        channel_radius = {}
        scattering_radius = {}
        q_value = {}
        records = []
        for l in range(NLS):
            items, values = get_list_record(file_obj)
            l_value = items[2]
            awri = items[0]
            q_value[l_value] = items[1]
            competitive = items[3]

            # Calculate channel radius from ENDF-102 equation D.14
            a = Polynomial((0.123 * (NEUTRON_MASS*awri)**(1./3.) + 0.08,))

            # Construct scattering and channel radius
            if nro == 0:
                scattering_radius[l_value] = ap
                if naps == 0:
                    channel_radius[l_value] = a
                elif naps == 1:
                    channel_radius[l_value] = ap
            elif nro == 1:
                scattering_radius[l_value] = ape
                if naps == 0:
                    channel_radius[l_value] = a
                elif naps == 1:
                    channel_radius[l_value] = ape
                elif naps == 2:
                    channel_radius[l_value] = ap

            energy = values[0::6]
            spin = values[1::6]
            gt = np.asarray(values[2::6])
            gn = np.asarray(values[3::6])
            gg = np.asarray(values[4::6])
            gf = np.asarray(values[5::6])
            if competitive > 0:
                gx = gt - (gn + gg + gf)
            else:
                gx = np.zeros_like(gt)

            for i, E in enumerate(energy):
                records.append([energy[i], l_value, spin[i], gt[i], gn[i],
                                gg[i], gf[i], gx[i]])

        columns = ['energy', 'L', 'J', 'totalWidth', 'neutronWidth',
                   'captureWidth', 'fissionWidth', 'competitiveWidth']
        parameters = pd.DataFrame.from_records(records, columns=columns)

        # Create instance of class
        mlbw = cls(target_spin, energy_min, energy_max,
                   channel_radius, scattering_radius)
        mlbw.q_value = q_value
        mlbw.atomic_weight_ratio = awri
        mlbw.parameters = parameters

        return mlbw

    def _prepare_resonances(self):
        df = self.parameters.copy()

        # Penetration and shift factors
        p = np.zeros(len(df))
        s = np.zeros(len(df))

        # Penetration and shift factors for competitive reaction
        px = np.zeros(len(df))
        sx = np.zeros(len(df))

        l_values = []
        competitive = []

        A = self.atomic_weight_ratio
        for i, E, l, J, gt, gn, gg, gf, gx in df.itertuples():
            if l not in l_values:
                l_values.append(l)
                competitive.append(gx > 0)

            # Determine penetration and shift corresponding to resonance energy
            k = wave_number(A, E)
            rho = k*self.channel_radius[l](E)
            rhohat = k*self.scattering_radius[l](E)
            p[i], s[i] = penetration_shift(l, rho)

            # Determine penetration at modified energy for competitive reaction
            if gx > 0:
                Ex = E + self.q_value[l]*(A + 1)/A
                rho = k*self.channel_radius[l](Ex)
                rhohat = k*self.scattering_radius[l](Ex)
                px[i], sx[i] = penetration_shift(l, rho)
            else:
                px[i] = sx[i] = 0.0

        df['p'] = p
        df['s'] = s
        df['px'] = px
        df['sx'] = sx

        self._l_values = np.array(l_values)
        self._competitive = np.array(competitive)
        for l in l_values:
            self._parameter_matrix[l] = df[df.L == l].values

        self._prepared = True


class SingleLevelBreitWigner(MultiLevelBreitWigner):
    """Single-level Breit-Wigner resolved resonance formalism data.

    Single-level Breit-Wigner resolved resonance data is is identified by LRF=1
    in the ENDF-6 format.

    Parameters
    ----------
    target_spin : float
        Intrinsic spin, :math:`I`, of the target nuclide
    energy_min : float
        Minimum energy of the resolved resonance range in eV
    energy_max : float
        Maximum energy of the resolved resonance range in eV
    channel : dict
        Dictionary whose keys are l-values and values are channel radii as a
        function of energy
    scattering : dict
        Dictionary whose keys are l-values and values are scattering radii as a
        function of energy

    Attributes
    ----------
    atomic_weight_ratio : float
        Atomic weight ratio of the target nuclide given as a function of
        l-value. Note that this may be different than the value for the
        evaluation as a whole.
    channel_radius : dict
        Dictionary whose keys are l-values and values are channel radii as a
        function of energy
    energy_max : float
        Maximum energy of the resolved resonance range in eV
    energy_min : float
        Minimum energy of the resolved resonance range in eV
    parameters : pandas.DataFrame
        Energies, spins, and resonances widths for each resonance
    q_value : dict
        Q-value to be added to incident particle's center-of-mass energy to
        determine the channel energy for use in the penetrability factor. The
        keys of the dictionary are l-values.
    scattering_radius : dict
        Dictionary whose keys are l-values and values are scattering radii as a
        function of energy
    target_spin : float
        Intrinsic spin, :math:`I`, of the target nuclide

    """

    def __init__(self, target_spin, energy_min, energy_max, channel, scattering):
        super().__init__(target_spin, energy_min, energy_max, channel,
                         scattering)

        # Set resonance reconstruction function
        if _reconstruct:
            self._reconstruct = reconstruct_slbw
        else:
            self._reconstruct = None


class ReichMoore(ResonanceRange):
    """Reich-Moore resolved resonance formalism data.

    Reich-Moore resolved resonance data is identified by LRF=3 in the ENDF-6
    format.

    Parameters
    ----------
    target_spin : float
        Intrinsic spin, :math:`I`, of the target nuclide
    energy_min : float
        Minimum energy of the resolved resonance range in eV
    energy_max : float
        Maximum energy of the resolved resonance range in eV
    channel : dict
        Dictionary whose keys are l-values and values are channel radii as a
        function of energy
    scattering : dict
        Dictionary whose keys are l-values and values are scattering radii as a
        function of energy

    Attributes
    ----------
    angle_distribution : bool
        Indicate whether parameters can be used to compute angular distributions
    atomic_weight_ratio : float
        Atomic weight ratio of the target nuclide given as a function of
        l-value. Note that this may be different than the value for the
        evaluation as a whole.
    channel_radius : dict
        Dictionary whose keys are l-values and values are channel radii as a
        function of energy
    energy_max : float
        Maximum energy of the resolved resonance range in eV
    energy_min : float
        Minimum energy of the resolved resonance range in eV
    num_l_convergence : int
        Number of l-values which must be used to converge the calculation
    scattering_radius : dict
        Dictionary whose keys are l-values and values are scattering radii as a
        function of energy
    parameters : pandas.DataFrame
        Energies, spins, and resonances widths for each resonance
    target_spin : float
        Intrinsic spin, :math:`I`, of the target nuclide

    """

    def __init__(self, target_spin, energy_min, energy_max, channel, scattering):
        super().__init__(target_spin, energy_min, energy_max, channel,
                         scattering)
        self.parameters = None
        self.angle_distribution = False
        self.num_l_convergence = 0

        # Set resonance reconstruction function
        if _reconstruct:
            self._reconstruct = reconstruct_rm
        else:
            self._reconstruct = None

    @classmethod
    def from_endf(cls, ev, file_obj, items):
        """Create Reich-Moore resonance data from an ENDF evaluation.

        Parameters
        ----------
        ev : openmc.data.endf.Evaluation
            ENDF evaluation
        file_obj : file-like object
            ENDF file positioned at the second record of a resonance range
            subsection in MF=2, MT=151
        items : list
            Items from the CONT record at the start of the resonance range
            subsection

        Returns
        -------
        openmc.data.ReichMoore
            Reich-Moore resonance parameters

        """
        # Read energy-dependent scattering radius if present
        energy_min, energy_max = items[0:2]
        nro, naps = items[4:6]
        if nro != 0:
            params, ape = get_tab1_record(file_obj)

        # Other scatter radius parameters
        items = get_cont_record(file_obj)
        target_spin = items[0]
        ap = Polynomial((items[1],))
        angle_distribution = (items[3] == 1)  # Flag for angular distribution
        NLS = items[4]  # Number of l-values
        num_l_convergence = items[5]  # Number of l-values for convergence

        # Read resonance widths, J values, etc
        channel_radius = {}
        scattering_radius = {}
        records = []
        for i in range(NLS):
            items, values = get_list_record(file_obj)
            apl = Polynomial((items[1],)) if items[1] != 0.0 else ap
            l_value = items[2]
            awri = items[0]

            # Calculate channel radius from ENDF-102 equation D.14
            a = Polynomial((0.123 * (NEUTRON_MASS*awri)**(1./3.) + 0.08,))

            # Construct scattering and channel radius
            if nro == 0:
                scattering_radius[l_value] = apl
                if naps == 0:
                    channel_radius[l_value] = a
                elif naps == 1:
                    channel_radius[l_value] = apl
            elif nro == 1:
                if naps == 0:
                    channel_radius[l_value] = a
                    scattering_radius[l_value] = ape
                elif naps == 1:
                    channel_radius[l_value] = scattering_radius[l_value] = ape
                elif naps == 2:
                    channel_radius[l_value] = apl
                    scattering_radius[l_value] = ape

            energy = values[0::6]
            spin = values[1::6]
            gn = values[2::6]
            gg = values[3::6]
            gfa = values[4::6]
            gfb = values[5::6]

            for i, E in enumerate(energy):
                records.append([energy[i], l_value, spin[i], gn[i], gg[i],
                                gfa[i], gfb[i]])

        # Create pandas DataFrame with resonance data
        columns = ['energy', 'L', 'J', 'neutronWidth', 'captureWidth',
                   'fissionWidthA', 'fissionWidthB']
        parameters = pd.DataFrame.from_records(records, columns=columns)

        # Create instance of ReichMoore
        rm = cls(target_spin, energy_min, energy_max,
                 channel_radius, scattering_radius)
        rm.parameters = parameters
        rm.angle_distribution = angle_distribution
        rm.num_l_convergence = num_l_convergence
        rm.atomic_weight_ratio = awri

        return rm

    def _prepare_resonances(self):
        df = self.parameters.copy()

        # Penetration and shift factors
        p = np.zeros(len(df))
        s = np.zeros(len(df))

        l_values = []
        lj_values = []

        A = self.atomic_weight_ratio
        for i, E, l, J, gn, gg, gfa, gfb in df.itertuples():
            if l not in l_values:
                l_values.append(l)
            if (l, abs(J)) not in lj_values:
                lj_values.append((l, abs(J)))

            # Determine penetration and shift corresponding to resonance energy
            k = wave_number(A, E)
            rho = k*self.channel_radius[l](E)
            rhohat = k*self.scattering_radius[l](E)
            p[i], s[i] = penetration_shift(l, rho)

        df['p'] = p
        df['s'] = s

        self._l_values = np.array(l_values)
        for (l, J) in lj_values:
            self._parameter_matrix[l, J] = df[(df.L == l) &
                                              (abs(df.J) == J)].values

        self._prepared = True


class RMatrixLimited(ResonanceRange):
    """R-matrix limited resolved resonance formalism data.

    R-matrix limited resolved resonance data is identified by LRF=7 in the
    ENDF-6 format.

    Parameters
    ----------
    energy_min : float
        Minimum energy of the resolved resonance range in eV
    energy_max : float
        Maximum energy of the resolved resonance range in eV
    particle_pairs : list of dict
        List of particle pairs. Each particle pair is represented by a
        dictionary that contains the mass, atomic number, spin, and parity of
        each particle as well as other characteristics.
    spin_groups : list of dict
        List of spin groups. Each spin group is characterized by channels,
        resonance energies, and resonance widths.

    Attributes
    ----------
    reduced_width : bool
        Flag indicating whether channel widths in eV or reduced-width amplitudes
        in eV^1/2 are given
    formalism : int
        Flag to specify which formulae for the R-matrix are to be used
    particle_pairs : list of dict
        List of particle pairs. Each particle pair is represented by a
        dictionary that contains the mass, atomic number, spin, and parity of
        each particle as well as other characteristics.
    spin_groups : list of dict
        List of spin groups. Each spin group is characterized by channels,
        resonance energies, and resonance widths.

    """

    def __init__(self, energy_min, energy_max, particle_pairs, spin_groups):
        super().__init__(0.0, energy_min, energy_max, None, None)
        self.reduced_width = False
        self.formalism = 3
        self.particle_pairs = particle_pairs
        self.spin_groups = spin_groups

    @classmethod
    def from_endf(cls, ev, file_obj, items):
        """Read R-Matrix limited resonance data from an ENDF evaluation.

        Parameters
        ----------
        ev : openmc.data.endf.Evaluation
            ENDF evaluation
        file_obj : file-like object
            ENDF file positioned at the second record of a resonance range
            subsection in MF=2, MT=151
        items : list
            Items from the CONT record at the start of the resonance range
            subsection

        Returns
        -------
        openmc.data.RMatrixLimited
            R-matrix limited resonance parameters

        """
        energy_min, energy_max = items[0:2]

        items = get_cont_record(file_obj)
        reduced_width = (items[2] == 1)  # reduced width amplitude?
        formalism = items[3]  # Specify which formulae are used
        n_spin_groups = items[4]  # Number of Jpi values (NJS)

        particle_pairs = []
        spin_groups = []

        items, values = get_list_record(file_obj)
        n_pairs = items[5]//2  # Number of particle pairs (NPP)
        for i in range(n_pairs):
            first = {'mass': values[12*i],
                     'z': int(values[12*i + 2]),
                     'spin': values[12*i + 4],
                     'parity': values[12*i + 10]}
            second = {'mass': values[12*i + 1],
                      'z': int(values[12*i + 3]),
                      'spin': values[12*i + 5],
                      'parity': values[12*i + 11]}

            q_value = values[12*i + 6]
            penetrability = values[12*i + 7]
            shift = values[12*i + 8]
            mt = int(values[12*i + 9])

            particle_pairs.append(ParticlePair(
                first, second, q_value, penetrability, shift, mt))

        # loop over spin groups
        for i in range(n_spin_groups):
            items, values = get_list_record(file_obj)
            J = items[0]
            if J == 0.0:
                parity = '+' if items[1] == 1.0 else '-'
            else:
                parity = '+' if J > 0. else '-'
                J = abs(J)
            kbk = items[2]
            kps = items[3]
            n_channels = items[5]
            channels = []
            for j in range(n_channels):
                channel = {}
                channel['particle_pair'] = particle_pairs[
                    int(values[6*j]) - 1]
                channel['l'] = values[6*j + 1]
                channel['spin'] = values[6*j + 2]
                channel['boundary'] = values[6*j + 3]
                channel['effective_radius'] = values[6*j + 4]
                channel['true_radius'] = values[6*j + 5]
                channels.append(channel)

            # Read resonance energies and widths
            items, values = get_list_record(file_obj)
            n_resonances = items[3]
            records = []
            m = n_channels//6 + 1
            for j in range(n_resonances):
                energy = values[6*m*j]
                records.append([energy] + [values[6*m*j + k + 1]
                                           for k in range(n_channels)])

            # Determine column names
            columns = ['energy']
            for channel in channels:
                mt = channel['particle_pair'].mt
                if mt == 2:
                    columns.append('neutronWidth')
                elif mt == 18:
                    columns.append('fissionWidth')
                elif mt == 102:
                    columns.append('captureWidth')
                else:
                    columns.append('width (MT={})'.format(mt))

            # Create Pandas dataframe with resonance parameters
            parameters = pd.DataFrame.from_records(records, columns=columns)

            # Construct SpinGroup instance and add to list
            sg = SpinGroup(J, parity, channels, parameters)
            spin_groups.append(sg)

            # Optional extension (Background R-Matrix)
            if kbk > 0:
                items, values = get_list_record(file_obj)
                lbk = items[4]
                if lbk == 1:
                    params, rbr = get_tab1_record(file_obj)
                    params, rbi = get_tab1_record(file_obj)

            # Optional extension (Tabulated phase shifts)
            if kps > 0:
                items, values = get_list_record(file_obj)
                lps = items[4]
                if lps == 1:
                    params, psr = get_tab1_record(file_obj)
                    params, psi = get_tab1_record(file_obj)

        rml = cls(energy_min, energy_max, particle_pairs, spin_groups)
        rml.reduced_width = reduced_width
        rml.formalism = formalism

        return rml


class ParticlePair:
    def __init__(self, first, second, q_value, penetrability,
                 shift, mt):
        self.first = first
        self.second = second
        self.q_value = q_value
        self.penetrability = penetrability
        self.shift = shift
        self.mt = mt


class SpinGroup:
    """Resonance spin group

    Attributes
    ----------
    spin : float
        Total angular momentum (nuclear spin)
    parity : {'+', '-'}
        Even (+) or odd(-) parity
    channels : list of openmc.data.Channel
        Available channels
    parameters : pandas.DataFrame
        Energies/widths for each resonance/channel

    """

    def __init__(self, spin, parity, channels, parameters):
        self.spin = spin
        self.parity = parity
        self.channels = channels
        self.parameters = parameters

    def __repr__(self):
        return '<SpinGroup: Jpi={}{}>'.format(self.spin, self.parity)


class Unresolved(ResonanceRange):
    """Unresolved resonance parameters as identified by LRU=2 in MF=2.

    Parameters
    ----------
    target_spin : float
        Intrinsic spin, :math:`I`, of the target nuclide
    energy_min : float
        Minimum energy of the unresolved resonance range in eV
    energy_max : float
        Maximum energy of the unresolved resonance range in eV
    channel : openmc.data.Function1D
        Channel radii as a function of energy
    scattering : openmc.data.Function1D
        Scattering radii as a function of energy

    Attributes
    ----------
    add_to_background : bool
        If True, file 3 contains partial cross sections to be added to the
        average unresolved cross sections calculated from parameters.
    atomic_weight_ratio : float
        Atomic weight ratio of the target nuclide
    channel_radius : openmc.data.Function1D
        Channel radii as a function of energy
    energies : Iterable of float
        Energies at which parameters are tabulated
    energy_max : float
        Maximum energy of the unresolved resonance range in eV
    energy_min : float
        Minimum energy of the unresolved resonance range in eV
    parameters : list of pandas.DataFrame
        Average resonance parameters at each energy
    scattering_radius : openmc.data.Function1D
        Scattering radii as a function of energy
    target_spin : float
        Intrinsic spin, :math:`I`, of the target nuclide

    """

    def __init__(self, target_spin, energy_min, energy_max, channel, scattering):
        super().__init__(target_spin, energy_min, energy_max, channel,
                         scattering)
        self.energies = None
        self.parameters = None
        self.add_to_background = False
        self.atomic_weight_ratio = None

    @classmethod
    def from_endf(cls, file_obj, items, fission_widths):
        """Read unresolved resonance data from an ENDF evaluation.

        Parameters
        ----------
        file_obj : file-like object
            ENDF file positioned at the second record of a resonance range
            subsection in MF=2, MT=151
        items : list
            Items from the CONT record at the start of the resonance range
            subsection
        fission_widths : bool
            Whether fission widths are given

        Returns
        -------
        openmc.data.Unresolved
            Unresolved resonance region parameters

        """
        # Read energy-dependent scattering radius if present
        energy_min, energy_max = items[0:2]
        nro, naps = items[4:6]
        if nro != 0:
            params, ape = get_tab1_record(file_obj)

        # Get SPI, AP, and LSSF
        formalism = items[3]
        if not (fission_widths and formalism == 1):
            items = get_cont_record(file_obj)
            target_spin = items[0]
            if nro == 0:
                ap = Polynomial((items[1],))
            add_to_background = (items[2] == 0)

        if not fission_widths and formalism == 1:
            # Case A -- fission widths not given, all parameters are
            # energy-independent
            NLS = items[4]
            columns = ['L', 'J', 'd', 'amun', 'gn0', 'gg']
            records = []
            for ls in range(NLS):
                items, values = get_list_record(file_obj)
                awri = items[0]
                l = items[2]
                NJS = items[5]
                for j in range(NJS):
                    d, j, amun, gn0, gg = values[6*j:6*j + 5]
                    records.append([l, j, d, amun, gn0, gg])
            parameters = pd.DataFrame.from_records(records, columns=columns)
            energies = None

        elif fission_widths and formalism == 1:
            # Case B -- fission widths given, only fission widths are
            # energy-dependent
            items, energies = get_list_record(file_obj)
            target_spin = items[0]
            if nro == 0:
                ap = Polynomial((items[1],))
            add_to_background = (items[2] == 0)
            NE, NLS = items[4:6]
            records = []
            columns = ['L', 'J', 'E', 'd', 'amun', 'amuf', 'gn0', 'gg', 'gf']
            for ls in range(NLS):
                items = get_cont_record(file_obj)
                awri = items[0]
                l = items[2]
                NJS = items[4]
                for j in range(NJS):
                    items, values = get_list_record(file_obj)
                    muf = items[3]
                    d = values[0]
                    j = values[1]
                    amun = values[2]
                    gn0 = values[3]
                    gg = values[4]
                    gfs = values[6:]
                    for E, gf in zip(energies, gfs):
                        records.append([l, j, E, d, amun, muf, gn0, gg, gf])
            parameters = pd.DataFrame.from_records(records, columns=columns)

        elif formalism == 2:
            # Case C -- all parameters are energy-dependent
            NLS = items[4]
            columns = ['L', 'J', 'E', 'd', 'amux', 'amun', 'amuf', 'gx', 'gn0',
                       'gg', 'gf']
            records = []
            for ls in range(NLS):
                items = get_cont_record(file_obj)
                awri = items[0]
                l = items[2]
                NJS = items[4]
                for j in range(NJS):
                    items, values = get_list_record(file_obj)
                    ne = items[5]
                    j = items[0]
                    amux = values[2]
                    amun = values[3]
                    amuf = values[5]
                    energies = []
                    for k in range(1, ne + 1):
                        E = values[6*k]
                        d = values[6*k + 1]
                        gx = values[6*k + 2]
                        gn0 = values[6*k + 3]
                        gg = values[6*k + 4]
                        gf = values[6*k + 5]
                        energies.append(E)
                        records.append([l, j, E, d, amux, amun, amuf, gx, gn0,
                                        gg, gf])
            parameters = pd.DataFrame.from_records(records, columns=columns)

        # Calculate channel radius from ENDF-102 equation D.14
        a = Polynomial((0.123 * (NEUTRON_MASS*awri)**(1./3.) + 0.08,))

        # Determine scattering and channel radius
        if nro == 0:
            scattering_radius = ap
            if naps == 0:
                channel_radius = a
            elif naps == 1:
                channel_radius = ap
        elif nro == 1:
            scattering_radius = ape
            if naps == 0:
                channel_radius = a
            elif naps == 1:
                channel_radius = ape
            elif naps == 2:
                channel_radius = ap

        urr = cls(target_spin, energy_min, energy_max, channel_radius,
                  scattering_radius)
        urr.parameters = parameters
        urr.add_to_background = add_to_background
        urr.atomic_weight_ratio = awri
        urr.energies = energies

        return urr


_FORMALISMS = {0: ResonanceRange,
               1: SingleLevelBreitWigner,
               2: MultiLevelBreitWigner,
               3: ReichMoore,
               7: RMatrixLimited}

_RESOLVED = (SingleLevelBreitWigner, MultiLevelBreitWigner,
             ReichMoore, RMatrixLimited)
