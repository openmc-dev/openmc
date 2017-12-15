from collections import defaultdict, MutableSequence, Iterable
import io

import numpy as np
from numpy.polynomial import Polynomial
import pandas as pd

from .data import NEUTRON_MASS
from .endf import get_head_record, get_cont_record, get_tab1_record, get_list_record
import openmc.checkvalue as cv
from .resonance import ResonanceRange

class ResonanceCovariance(object):
    """Resolved resonance covariance data

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

    @ranges.setter
    def ranges(self, ranges):
        cv.check_type('resonance ranges', ranges, MutableSequence)
        self._ranges = cv.CheckedList(ResonanceRange, 'resonance ranges',
                                      ranges)

    @classmethod
    def from_endf(cls, ev):
        """Generate resonance covariance data from an ENDF evaluation.

        Parameters
        ----------
        ev : openmc.data.endf.Evaluation
            ENDF evaluation

        Returns
        -------
        openmc.data.ResonanceCovariance
            Resonance covariance data

        """
        file_obj = io.StringIO(ev.section[32, 151])

        # Determine whether discrete or continuous representation
        items = get_head_record(file_obj)
        n_isotope = items[4] # Number of isotopes

        ranges = []
        for iso in range(n_isotope):
            items = get_cont_record(file_obj)
            abundance = items[1]
            fission_widths = (items[3] == 1) # fission widths are given?
            n_ranges = items[4] # number of resonance energy ranges

            for j in range(n_ranges):
                items = get_cont_record(file_obj)
                resonance_flag = items[2]  # flag for resolved (1)/unresolved (2)
                formalism = items[3]  # resonance formalism

                # Throw error for unsupported formalisms
                if formalism in [0,1,2,7]:
                    raise TypeError('LRF= ', formalism,
                                    ' covariance not supported for this formalism')

                if resonance_flag in (0, 1):
                    # resolved resonance region
                    erange = _FORMALISMS[formalism].from_endf(ev, file_obj, items)

                elif resonance_flag == 2:
                    raise TypeError('Unresolved resonance not supported')

                #erange.material = self
                ranges.append(erange)

        return cls(ranges)


class ReichMooreCovariance(ResonanceRange):
    """Reich-Moore resolved resonance formalism covariance data.

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
    cov_paramaters: list
        The parameters that are included in the covariance matrix
    covariance_matrix : array
        The covariance matrix contained within the ENDF evaluation
    
    
    """

    def __init__(self, energy_min, energy_max):
        self.parameters = None
        self.covariance = None

    @classmethod
    def from_endf(cls, ev, file_obj, items):
        """Create Reich-Moore resonance covariance data from an ENDF evaluation.
        Includes the resonance parameters contained separately in File 32.

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
        openmc.data.ReichMooreCovariance
            Reich-Moore resonance covariance parameters

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
        LCOMP = items[3]  # Flag for compatibility 0,1,2 - 2 is compact form
        NLS = items[4]  # Number of l-values

        
        # Build covariance matrix for General Resolved Resonance Formats
        if LCOMP == 1:
            items = get_cont_record(file_obj)
            awri = items[0]
            num_short_range = items[4] #Number of short range type resonance 
                                       #covariances
            num_long_range = items[5] #Number of long range type resonance
                                      #covariances
            # Read resonance widths, J values, etc
            channel_radius = {}
            scattering_radius = {}
            records = []
            for i in range(num_short_range):
                items, values = get_list_record(file_obj)
                num_parameters = items[2]
                num_res = items[5]
                num_par_vals = num_res*6
                res_values = values[:num_par_vals]
                cov_values = values[num_par_vals:]
                
                energy = res_values[0::6]
                spin = res_values[1::6]
                gn = res_values[2::6]
                gg = res_values[3::6]
                gfa = res_values[4::6]
                gfb = res_values[5::6]

                for i, E in enumerate(energy):
                    records.append([energy[i], spin[i], gn[i], gg[i],
                                    gfa[i], gfb[i]])

                #Build the upper-triangular covariance matrix
                cov_dim = num_parameters*num_res
                cov = np.zeros([cov_dim,cov_dim])
                indices = np.triu_indices(cov_dim)
                cov[indices] = cov_values

            # Create pandas DataFrame with resonance data
            columns = ['energy', 'J', 'neutronWidth', 'captureWidth',
                       'fissionWidthA', 'fissionWidthB']
            parameters = pd.DataFrame.from_records(records, columns=columns)

            # Create instance of ReichMooreCovariance
            rmc = cls(energy_min, energy_max)
            rmc.parameters = parameters
            rmc.covariance = cov

            return rmc

        elif LCOMP in [0,2]:
            TypeError('LCOMP = ',LCOMP,' not supported')


# _FORMALISMS = {0: ResonanceRange,
#                1: SingleLevelBreitWigner,
#                2: MultiLevelBreitWigner,
#                3: ReichMoore,
#                7: RMatrixLimited}
_FORMALISMS = {3: ReichMooreCovariance}
                       

