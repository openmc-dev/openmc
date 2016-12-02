from six import string_types

from numbers import Integral, Real

import h5py
import numpy as np
from scipy.special import wofz
from six import string_types

from . import WMP_VERSION
from .data import K_BOLTZMANN
import openmc.checkvalue as cv
from openmc.mixin import EqualityMixin


# Formalisms
FORM_MLBW = 2
FORM_RM   = 3
FORM_RML  = 7

# Constants that determine which value to access
MP_EA = 0       # Pole

# Reich-Moore indices
RM_RT = 1       # Residue total
RM_RA = 2       # Residue absorption
RM_RF = 3       # Residue fission

# Multi-level Breit Wigner indices
MLBW_RT = 1     # Residue total
MLBW_RX = 2     # Residue compettitive
MLBW_RA = 3     # Residue absorption
MLBW_RF = 4     # Residue fission

# Polynomial fit indices
FIT_T = 0       # Total
FIT_A = 1       # Absorption
FIT_F = 2       # Fission


def _faddeeva(z):
    """Evaluate the complex Faddeeva function.

    Technically, the value we want is given by the equation:
    w(z) = I/Pi * Integrate[Exp[-t^2]/(z-t), {t, -Infinity, Infinity}]
    as shown in Equation 63 from Hwang, R. N. "A rigorous pole
    representation of multilevel cross sections and its practical
    applications." Nuclear Science and Engineering 96.3 (1987): 192-209.

    The scipy.special.wofz function evaluates w(z) = exp(-z^2)erfc(-iz). These
    two forms of the Faddeeva function are related by a transformation.

    If we call the integral form w_int, and the function form w_fun:
    For imag(z) > 0, w_int(z) = w_fun(z)
    For imag(z) < 0, w_int(z) = -conjg(w_fun(conjg(z)))

    Parameters
    ----------
    z : Complex
        Argument to the Faddeeva function.

    Returns
    -------
    Complex
        I/Pi * Integrate[Exp[-t^2]/(z-t), {t, -Infinity, Infinity}]

    """
    if np.angle(z) > 0:
        return wofz(z)
    else:
        return -np.conj(wofz(z))


def _broaden_wmp_polynomials(En, dopp, n):
    """Evaluate Doppler-broadened windowed multipole curvefit.

    The curvefit is a polynomial of the form
    a/En + b/sqrt(En) + c + d sqrt(En) ...

    Parameters
    ----------
    En : Real
        Energy to evaluate at.
    dopp : Real
        sqrt(atomic weight ratio / kT) in units of eV.
    n : Integral
        Number of components to the polynomial.

    Returns
    -------
    np.ndarray
        The value of each Doppler-broadened curvefit polynomial term.

    """
    sqrtE = np.sqrt(En)
    beta = sqrtE * dopp
    half_inv_dopp2 = 0.5 / dopp**2
    quarter_inv_dopp4 = half_inv_dopp2**2

    if beta > 6.0:
        # Save time, ERF(6) is 1 to machine precision.
        # beta/sqrtpi*exp(-beta**2) is also approximately 1 machine epsilon.
        erfBeta = 1.0
        exp_m_beta2 = 0.0
    else:
        erfBeta = np.erf(beta)
        exp_m_beta2 = np.exp(-beta**2)

    # Assume that, for sure, we'll use a second order (1/E, 1/V, const)
    # fit, and no less.

    factors = np.zeros(n)

    factors[0] = erfBeta / En
    factors[1] = 1.0 / sqrtE
    factors[2] = (factors[0] * (half_inv_dopp2 + En)
                  + exp_m_beta2 / (beta * np.sqrt(np.pi)))

    # Perform recursive broadening of high order components.
    for i in range(1, n-2):
        if i != 1:
            factors[i+2] = (-factors[i-2] * (i - 1.0) * i * quarter_inv_dopp4
                + factors[i] * (En + (1.0 + 2.0 * i) * half_inv_dopp2))
        else:
            # Although it's mathematically identical, factors[0] will contain
            # nothing, and we don't want to have to worry about memory.
            factors[i+2] = factors[i]*(En + (1.0 + 2.0 * i) * half_inv_dopp2)

    return factors


class WindowedMultipole(EqualityMixin):
    """Resonant cross sections represented in the windowed multipole format.

    Attributes
    ----------
    num_l : Integral
        Number of possible l quantum states for this nuclide.
    fit_order : Integral
        Order of the windowed curvefit.
    fissionable : bool
        Whether or not the target nuclide has fission data.
    formalism : str
        The R-matrix formalism used to reconstruct resonances.  Either 'MLBW'
        for multi-level Breit Wigner or 'RM' for Reich-Moore.
    spacing : Real
        The width of each window in sqrt(E)-space.  For example, the frst window
        will end at (sqrt(start_E) + spacing)**2 and the second window at
        (sqrt(start_E) + 2*spacing)**2.
    sqrtAWR : Real
        Square root of the atomic weight ratio of the target nuclide.
    start_E : Real
        Lowest energy in eV the library is valid for.
    end_E : Real
        Highest energy in eV the library is valid for.
    data : np.ndarray
        A 2D array of complex poles and residues.  data[i, 0] gives the energy
        at which pole i is located.  data[i, 1:] gives the residues associated
        with the i-th pole.  There are 3 residues for Reich-Moore data, one each
        for the total, absorption, and fission channels.  Multi-level
        Breit Wigner data has an additional residue for the competitive channel.
    pseudo_k0RS : np.ndarray
        A 1D array of Real values.  There is one value for each valid l
        quantum number.  The values are equal to
        sqrt(2 m / hbar) * AWR / (AWR + 1) * r
        where m is the neutron mass, AWR is the atomic weight ratio, and r
        is the l-dependent scattering radius.
    l_value : np.ndarray
        A 1D array of Integral values equal to the l quantum number for each
        pole + 1.
    w_start : np.ndarray
        A 1D array of Integral values.  w_start[i] - 1 is the index of the first
        pole in window i.
    w_end : np.ndarray
        A 1D array of Integral values.  w_end[i] - 1 is the index of the last
        pole in window i.
    broaden_poly : np.ndarray
        A 1D array of boolean values indicating whether or not the polynomial
        curvefit in that window should be Doppler broadened.
    curvefit : np.ndarray
        A 3D array of Real curvefit polynomial coefficients.  curvefit[i, 0, :]
        gives coefficients for the total cross section in window i.
        curvefit[i, 1, :] gives absorption coefficients and curvefit[i, 2, :]
        gives fission coefficients.  The polynomial terms are increasing powers
        of sqrt(E) starting with 1/E e.g:
        a/E + b/sqrt(E) + c + d sqrt(E) + ...

    """
    def __init__(self):
        self.num_l = None
        self.fit_order = None
        self.fissionable = None
        self.formalism = None
        self.spacing = None
        self.sqrtAWR = None
        self.start_E = None
        self.end_E = None
        self.data = None
        self.pseudo_k0RS = None
        self.l_value = None
        self.w_start = None
        self.w_end = None
        self.broaden_poly = None
        self.curvefit = None

    @property
    def num_l(self):
        return self._num_l

    @property
    def fit_order(self):
        return self._fit_order

    @property
    def fissionable(self):
        return self._fissionable

    @property
    def formalism(self):
        return self._formalism

    @property
    def spacing(self):
        return self._spacing

    @property
    def sqrtAWR(self):
        return self._sqrtAWR

    @property
    def start_E(self):
        return self._start_E

    @property
    def end_E(self):
        return self._end_E

    @property
    def data(self):
        return self._data

    @property
    def pseudo_k0RS(self):
        return self._pseudo_k0RS

    @property
    def l_value(self):
        return self._l_value

    @property
    def w_start(self):
        return self._w_start

    @property
    def w_end(self):
        return self._w_end

    @property
    def broaden_poly(self):
        return self._broaden_poly

    @property
    def curvefit(self):
        return self._curvefit

    @num_l.setter
    def num_l(self, num_l):
        if num_l is not None:
            cv.check_type('num_l', num_l, Integral)
            cv.check_greater_than('num_l', num_l, 1, equality=True)
        self._num_l = num_l

    @fit_order.setter
    def fit_order(self, fit_order):
        if fit_order is not None:
            cv.check_type('fit_order', fit_order, Integral)
            cv.check_greater_than('fit_order', fit_order, 1, equality=True)
        self._fit_order = fit_order

    @fissionable.setter
    def fissionable(self, fissionable):
        if fissionable is not None:
            cv.check_type('fissionable', fissionable, bool)
        self._fissionable = fissionable

    @formalism.setter
    def formalism(self, formalism):
        if formalism is not None:
            cv.check_type('formalism', formalism, string_types)
            cv.check_value('formalism', formalism, ('MLBW', 'RM'))
        self._formalism = formalism

    @spacing.setter
    def spacing(self, spacing):
        if spacing is not None:
            cv.check_type('spacing', spacing, Real)
            cv.check_greater_than('spacing', spacing, 0.0, equality=False)
        self._spacing = spacing

    @sqrtAWR.setter
    def sqrtAWR(self, sqrtAWR):
        if sqrtAWR is not None:
            cv.check_type('sqrtAWR', sqrtAWR, Real)
            cv.check_greater_than('sqrtAWR', sqrtAWR, 0.0, equality=False)
        self._sqrtAWR = sqrtAWR

    @start_E.setter
    def start_E(self, start_E):
        if start_E is not None:
            cv.check_type('start_E', start_E, Real)
            cv.check_greater_than('start_E', start_E, 0.0, equality=True)
        self._start_E = start_E

    @end_E.setter
    def end_E(self, end_E):
        if end_E is not None:
            cv.check_type('end_E', end_E, Real)
            cv.check_greater_than('end_E', end_E, 0.0, equality=False)
        self._end_E = end_E

    @data.setter
    def data(self, data):
        if data is not None:
            cv.check_type('data', data, np.ndarray)
            if len(data.shape) != 2:
                raise ValueError('Multipole data arrays must be 2D')
            if data.shape[1] not in (4, 5):  # 4 for RM, 5 for MLBW
                raise ValueError('The second dimension of multipole data arrays'
                                 ' must have a length of either 4 or 5')
            if not np.issubdtype(data.dtype, np.complex):
                raise TypeError('Multipole data arrays must be complex dtype')
        self._data = data

    @pseudo_k0RS.setter
    def pseudo_k0RS(self, pseudo_k0RS):
        if pseudo_k0RS is not None:
            cv.check_type('pseudo_k0RS', pseudo_k0RS, np.ndarray)
            if len(pseudo_k0RS.shape) != 1:
                raise ValueError('Multipole pseudo_k0RS arrays must be 1D')
            if not np.issubdtype(pseudo_k0RS.dtype, np.float):
                raise TypeError('Multipole data arrays must be float dtype')
        self._pseudo_k0RS = pseudo_k0RS

    @l_value.setter
    def l_value(self, l_value):
        if l_value is not None:
            cv.check_type('l_value', l_value, np.ndarray)
            if len(l_value.shape) != 1:
                raise ValueError('Multipole l_value arrays must be 1D')
            if not np.issubdtype(l_value.dtype, np.integer):
                raise TypeError('Multipole l_value arrays must be integer'
                                ' dtype')
        self._l_value = l_value

    @w_start.setter
    def w_start(self, w_start):
        if w_start is not None:
            cv.check_type('w_start', w_start, np.ndarray)
            if len(w_start.shape) != 1:
                raise ValueError('Multipole w_start arrays must be 1D')
            if not np.issubdtype(w_start.dtype, np.integer):
                raise TypeError('Multipole w_start arrays must be integer'
                                ' dtype')
        self._w_start = w_start

    @w_end.setter
    def w_end(self, w_end):
        if w_end is not None:
            cv.check_type('w_end', w_end, np.ndarray)
            if len(w_end.shape) != 1:
                raise ValueError('Multipole w_end arrays must be 1D')
            if not np.issubdtype(w_end.dtype, np.integer):
                raise TypeError('Multipole w_end arrays must be integer dtype')
        self._w_end = w_end

    @broaden_poly.setter
    def broaden_poly(self, broaden_poly):
        if broaden_poly is not None:
            cv.check_type('broaden_poly', broaden_poly, np.ndarray)
            if len(broaden_poly.shape) != 1:
                raise ValueError('Multipole broaden_poly arrays must be 1D')
            if not np.issubdtype(broaden_poly.dtype, np.bool):
                raise TypeError('Multipole broaden_poly arrays must be boolean'
                                ' dtype')
        self._broaden_poly = broaden_poly

    @curvefit.setter
    def curvefit(self, curvefit):
        if curvefit is not None:
            cv.check_type('curvefit', curvefit, np.ndarray)
            if len(curvefit.shape) != 3:
                raise ValueError('Multipole curvefit arrays must be 3D')
            if curvefit.shape[2] != 3:  # One each for sigT, sigA, sigF
                raise ValueError('The third dimension of multipole curvefit'
                                 ' arrays must have a length of 3')
            if not np.issubdtype(curvefit.dtype, np.float):
                raise TypeError('Multipole curvefit arrays must be float dtype')
        self._curvefit = curvefit

    @classmethod
    def from_hdf5(cls, group_or_filename):
        if isinstance(group_or_filename, h5py.Group):
            group = group_or_filename
        else:
            h5file = h5py.File(group_or_filename, 'r')
            version = h5file['version'].value
            if version != WMP_VERSION:
                raise ValueError('The given WMP data uses version '
                    + str(version) + ' whereas your installation of the OpenMC '
                    'Python API expects version ' + str(WMP_VERSION))
            group = h5file['nuclide']

        out = cls()

        # Read scalar values.  Note that group['max_w'] is ignored.

        length = group['length'].value
        windows = group['windows'].value
        out.num_l = group['num_l'].value
        out.fit_order = group['fit_order'].value
        out.fissionable = bool(group['fissionable'].value)

        if group['formalism'].value == FORM_MLBW:
            out.formalism = 'MLBW'
        elif group['formalism'].value == FORM_RM:
            out.formalism = 'RM'
        else:
            raise ValueError('Unrecognized/Unsupported R-matrix formalism')

        out.spacing = group['spacing'].value
        out.sqrtAWR = group['sqrtAWR'].value
        out.start_E = group['start_E'].value
        out.end_E = group['end_E'].value

        # Read arrays.

        err = "WMP '{:}' array shape is not consistent with the '{:}' value"

        out.data = group['data'].value
        if out.data.shape[0] != length:
            raise ValueError(err.format('data', 'length'))

        out.pseudo_k0RS = group['pseudo_K0RS'].value
        if out.pseudo_k0RS.shape[0] != out.num_l:
            raise ValueError(err.format('pseudo_k0RS', 'num_l'))

        out.l_value = group['l_value'].value
        if out.l_value.shape[0] != length:
            raise ValueError(err.format('l_value', 'length'))

        out.w_start = group['w_start'].value
        if out.w_start.shape[0] != windows:
            raise ValueError(err.format('w_start', 'windows'))

        out.w_end = group['w_end'].value
        if out.w_end.shape[0] != windows:
            raise ValueError(err.format('w_end', 'windows'))

        out.broaden_poly = group['broaden_poly'].value.astype(np.bool)
        if out.broaden_poly.shape[0] != windows:
            raise ValueError(err.format('broaden_poly', 'windows'))

        out.curvefit = group['curvefit'].value
        if out.curvefit.shape[0] != windows:
            raise ValueError(err.format('curvefit', 'windows'))
        if out.curvefit.shape[1] != out.fit_order + 1:
            raise ValueError(err.format('curvefit', 'fit_order'))

        return out

    def _evaluate(self, E, T):
        """Return total, absorption, and fission XS at a single energy and temp.

        Parameters
        ----------
        E : Real
            Energy of the incident neutron in eV.
        T : Real
            Temperature of the target in K.

        Returns
        -------
        3-tuple of Real
            Total, absorption, and fission microscopic cross sections at the
            given energy and temperature.

        """

        # ======================================================================
        # Bookkeeping

        # Define some frequently used variables.
        sqrtkT = np.sqrt(K_BOLTZMANN * T)
        sqrtE = np.sqrt(E)
        invE = 1.0 / E
        dopp = self.sqrtAWR / sqrtkT

        # Locate us.
        i_window = int(np.floor((sqrtE - np.sqrt(self.start_E)) / self.spacing))
        startw = self.w_start[i_window] - 1
        endw = self.w_end[i_window]

        # Fill in factors.
        if startw <= endw:
          twophi = np.zeros(self.num_l, dtype=np.float)
          sigT_factor = np.zeros(self.num_l, dtype=np.cfloat)

          for iL in range(1, self.num_l+1):
              twophi[iL-1] = self.pseudo_k0RS[iL-1] * sqrtE
              if iL == 2:
                twophi[iL-1] = twophi[iL-1] - np.arctan(twophi[iL-1])
              elif iL == 3:
                arg = 3.0 * twophi[iL-1] / (3.0 - twophi[iL-1]**2)
                twophi[iL-1] = twophi[iL-1] - np.arctan(arg)
              elif iL == 4:
                arg = (twophi[iL-1] * (15.0 - twophi[iL-1]**2)
                       / (15.0 - 6.0 * twophi[iL-1]**2))
                twophi[iL-1] = twophi[iL-1] - np.arctan(arg)

          twophi = 2.0 * twophi
          sigT_factor = np.cos(twophi) - 1j*np.sin(twophi)

        # Initialize the ouptut cross sections.
        sigT = 0.0
        sigA = 0.0
        sigF = 0.0

        # ======================================================================
        # Add the contribution from the curvefit polynomial.

        if sqrtkT != 0 and self.broaden_poly[i_window]:
            # Broaden the curvefit.
            broadened_polynomials = _broaden_wmp_polynomials(E, dopp,
                                                             self.fit_order + 1)
            for i_poly in range(self.fit_order+1):
                sigT += (self.curvefit[i_window, i_poly, FIT_T]
                         * broadened_polynomials[i_poly])
                sigA += (self.curvefit[i_window, i_poly, FIT_A]
                         * broadened_polynomials[i_poly])
                sigF += (self.curvefit[i_window, i_poly, FIT_F]
                         * broadened_polynomials[i_poly])
        else:
            temp = invE
            for i_poly in range(self.fit_order+1):
                sigT += self.curvefit[i_window, i_poly, FIT_T] * temp
                sigA += self.curvefit[i_window, i_poly, FIT_A] * temp
                sigF += self.curvefit[i_window, i_poly, FIT_F] * temp
                temp *= sqrtE

        # ======================================================================
        # Add the contribution from the poles in this window.

        if sqrtkT == 0.0:
          # If at 0K, use asymptotic form.
          for i_pole in range(startw, endw):
              psi_chi = -1j / (self.data[i_pole, MP_EA] - sqrtE)
              c_temp = psi_chi / E
              if self.formalism == 'MLBW':
                  sigT += ((self.data[i_pole, MLBW_RT] * c_temp *
                            sigT_factor[self.l_value[i_pole]-1]).real
                           + (self.data[i_pole, MLBW_RX] * c_temp).real)
                  sigA += (self.data[i_pole, MLBW_RA] * c_temp).real
                  sigF += (self.data[i_pole, MLBW_RF] * c_temp).real
              elif self.formalism == 'RM':
                  sigT += (self.data[i_pole, RM_RT] * c_temp *
                           sigT_factor[self.l_value[i_pole]-1]).real
                  sigA += (self.data[i_pole, RM_RA] * c_temp).real
                  sigF += (self.data[i_pole, RM_RF] * c_temp).real
              else:
                  raise ValueError('Unrecognized/Unsupported R-matrix'
                                   ' formalism')

        else:
          # At temperature, use Faddeeva function-based form.
          for i_pole in range(startw, endw):
              Z = (sqrtE - self.data[i_pole, MP_EA]) * dopp
              w_val = _faddeeva(Z) * dopp * invE * np.sqrt(np.pi)
              if self.formalism == 'MLBW':
                  sigT += ((self.data[i_pole, MLBW_RT] *
                            sigT_factor[self.l_value[i_pole]-1] +
                            self.data[i_pole, MLBW_RX]) * w_val).real
                  sigA += (self.data[i_pole, MLBW_RA] * w_val).real
                  sigF += (self.data[i_pole, MLBW_RF] * w_val).real
              elif self.formalism == 'RM':
                  sigT += (self.data[i_pole, RM_RT] * w_val *
                           sigT_factor[self.l_value[i_pole]-1]).real
                  sigA += (self.data[i_pole, RM_RA] * w_val).real
                  sigF += (self.data[i_pole, RM_RF] * w_val).real
              else:
                  raise ValueError('Unrecognized/Unsupported R-matrix'
                                   ' formalism')

        return sigT, sigA, sigF

    def __call__(self, E, T):
        """Return total, absorption, and fission XS at given energy and temp.

        Parameters
        ----------
        E : Real or Iterable of Real
            Energy of the incident neutron in eV.
        T : Real
            Temperature of the target in K.

        Returns
        -------
        3-tuple of Real or 3-tuple of numpy.ndarray
            Total, absorption, and fission microscopic cross sections at the
            given energy and temperature.

        """

        fun = np.vectorize(lambda x: self._evaluate(x, T))
        return fun(E)
