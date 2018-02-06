from numbers import Integral, Real
from math import exp, erf, pi, sqrt

import h5py
import numpy as np

from . import WMP_VERSION
from .data import K_BOLTZMANN
import openmc.checkvalue as cv
from openmc.mixin import EqualityMixin


# Formalisms
_FORM_MLBW = 2
_FORM_RM   = 3

# Constants that determine which value to access
_MP_EA = 0       # Pole

# Reich-Moore indices
_RM_RT = 1       # Residue total
_RM_RA = 2       # Residue absorption
_RM_RF = 3       # Residue fission

# Multi-level Breit Wigner indices
_MLBW_RT = 1     # Residue total
_MLBW_RX = 2     # Residue compettitive
_MLBW_RA = 3     # Residue absorption
_MLBW_RF = 4     # Residue fission

# Polynomial fit indices
_FIT_T = 0       # Total
_FIT_A = 1       # Absorption
_FIT_F = 2       # Fission


def _faddeeva(z):
    r"""Evaluate the complex Faddeeva function.

    Technically, the value we want is given by the equation:

    .. math::
        w(z) = \frac{i}{\pi} \int_{-\infty}^{\infty} \frac{1}{z - t}
        \exp(-t^2) \text{d}t

    as shown in Equation 63 from Hwang, R. N. "A rigorous pole
    representation of multilevel cross sections and its practical
    applications." Nuclear Science and Engineering 96.3 (1987): 192-209.

    The :func:`scipy.special.wofz` function evaluates
    :math:`w(z) = \exp(-z^2) \text{erfc}(-iz)`. These two forms of the Faddeeva
    function are related by a transformation.

    If we call the integral form :math:`w_\text{int}`, and the function form
    :math:`w_\text{fun}`:

    .. math::
        w_\text{int}(z) =
        \begin{cases}
            w_\text{fun}(z) & \text{for } \text{Im}(z) > 0\\
            -w_\text{fun}(z^*)^* & \text{for } \text{Im}(z) < 0
        \end{cases}

    Parameters
    ----------
    z : complex
        Argument to the Faddeeva function.

    Returns
    -------
    complex
        :math:`\frac{i}{\pi} \int_{-\infty}^{\infty} \frac{1}{z - t} \exp(-t^2)
        \text{d}t`

    """
    from scipy.special import wofz
    if np.angle(z) > 0:
        return wofz(z)
    else:
        return -np.conj(wofz(z.conjugate()))


def _broaden_wmp_polynomials(E, dopp, n):
    r"""Evaluate Doppler-broadened windowed multipole curvefit.

    The curvefit is a polynomial of the form :math:`\frac{a}{E}
    + \frac{b}{\sqrt{E}} + c + d \sqrt{E} + \ldots`

    Parameters
    ----------
    E : Real
        Energy to evaluate at.
    dopp : Real
        sqrt(atomic weight ratio / kT) in units of eV.
    n : Integral
        Number of components to the polynomial.

    Returns
    -------
    numpy.ndarray
        The value of each Doppler-broadened curvefit polynomial term.

    """
    sqrtE = sqrt(E)
    beta = sqrtE * dopp
    half_inv_dopp2 = 0.5 / dopp**2
    quarter_inv_dopp4 = half_inv_dopp2**2

    if beta > 6.0:
        # Save time, ERF(6) is 1 to machine precision.
        # beta/sqrtpi*exp(-beta**2) is also approximately 1 machine epsilon.
        erf_beta = 1.0
        exp_m_beta2 = 0.0
    else:
        erf_beta = erf(beta)
        exp_m_beta2 = exp(-beta**2)

    # Assume that, for sure, we'll use a second order (1/E, 1/V, const)
    # fit, and no less.

    factors = np.zeros(n)

    factors[0] = erf_beta / E
    factors[1] = 1.0 / sqrtE
    factors[2] = (factors[0] * (half_inv_dopp2 + E)
                  + exp_m_beta2 / (beta * sqrt(pi)))

    # Perform recursive broadening of high order components. range(1, n-2)
    # replaces a do i = 1, n-3.  All indices are reduced by one due to the
    # 1-based vs. 0-based indexing.
    for i in range(1, n-2):
        if i != 1:
            factors[i+2] = (-factors[i-2] * (i - 1.0) * i * quarter_inv_dopp4
                + factors[i] * (E + (1.0 + 2.0 * i) * half_inv_dopp2))
        else:
            factors[i+2] = factors[i]*(E + (1.0 + 2.0 * i) * half_inv_dopp2)

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
    formalism : {'MLBW', 'RM'}
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
            cv.check_less_than('num_l', num_l, 4, equality=True)
            # There is an if block in _evaluate that assumes num_l <= 4.
        self._num_l = num_l

    @fit_order.setter
    def fit_order(self, fit_order):
        if fit_order is not None:
            cv.check_type('fit_order', fit_order, Integral)
            cv.check_greater_than('fit_order', fit_order, 2, equality=True)
            # _broaden_wmp_polynomials assumes the curve fit has at least 3
            # terms.
        self._fit_order = fit_order

    @fissionable.setter
    def fissionable(self, fissionable):
        if fissionable is not None:
            cv.check_type('fissionable', fissionable, bool)
        self._fissionable = fissionable

    @formalism.setter
    def formalism(self, formalism):
        if formalism is not None:
            cv.check_type('formalism', formalism, str)
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
            if data.shape[1] not in (3, 4, 5):  # 3 or 4 for RM, 4 or 5 for MLBW
                raise ValueError('The second dimension of multipole data arrays'
                                 ' must have a length of 3, 4 or 5')
            if not np.issubdtype(data.dtype, complex):
                raise TypeError('Multipole data arrays must be complex dtype')
        self._data = data

    @pseudo_k0RS.setter
    def pseudo_k0RS(self, pseudo_k0RS):
        if pseudo_k0RS is not None:
            cv.check_type('pseudo_k0RS', pseudo_k0RS, np.ndarray)
            if len(pseudo_k0RS.shape) != 1:
                raise ValueError('Multipole pseudo_k0RS arrays must be 1D')
            if not np.issubdtype(pseudo_k0RS.dtype, float):
                raise TypeError('Multipole data arrays must be float dtype')
        self._pseudo_k0RS = pseudo_k0RS

    @l_value.setter
    def l_value(self, l_value):
        if l_value is not None:
            cv.check_type('l_value', l_value, np.ndarray)
            if len(l_value.shape) != 1:
                raise ValueError('Multipole l_value arrays must be 1D')
            if not np.issubdtype(l_value.dtype, int):
                raise TypeError('Multipole l_value arrays must be integer'
                                ' dtype')
        self._l_value = l_value

    @w_start.setter
    def w_start(self, w_start):
        if w_start is not None:
            cv.check_type('w_start', w_start, np.ndarray)
            if len(w_start.shape) != 1:
                raise ValueError('Multipole w_start arrays must be 1D')
            if not np.issubdtype(w_start.dtype, int):
                raise TypeError('Multipole w_start arrays must be integer'
                                ' dtype')
        self._w_start = w_start

    @w_end.setter
    def w_end(self, w_end):
        if w_end is not None:
            cv.check_type('w_end', w_end, np.ndarray)
            if len(w_end.shape) != 1:
                raise ValueError('Multipole w_end arrays must be 1D')
            if not np.issubdtype(w_end.dtype, int):
                raise TypeError('Multipole w_end arrays must be integer dtype')
        self._w_end = w_end

    @broaden_poly.setter
    def broaden_poly(self, broaden_poly):
        if broaden_poly is not None:
            cv.check_type('broaden_poly', broaden_poly, np.ndarray)
            if len(broaden_poly.shape) != 1:
                raise ValueError('Multipole broaden_poly arrays must be 1D')
            if not np.issubdtype(broaden_poly.dtype, bool):
                raise TypeError('Multipole broaden_poly arrays must be boolean'
                                ' dtype')
        self._broaden_poly = broaden_poly

    @curvefit.setter
    def curvefit(self, curvefit):
        if curvefit is not None:
            cv.check_type('curvefit', curvefit, np.ndarray)
            if len(curvefit.shape) != 3:
                raise ValueError('Multipole curvefit arrays must be 3D')
            if curvefit.shape[2] not in (2, 3):  # sig_t, sig_a (maybe sig_f)
                raise ValueError('The third dimension of multipole curvefit'
                                 ' arrays must have a length of 2 or 3')
            if not np.issubdtype(curvefit.dtype, float):
                raise TypeError('Multipole curvefit arrays must be float dtype')
        self._curvefit = curvefit

    @classmethod
    def from_hdf5(cls, group_or_filename):
        """Construct a WindowedMultipole object from an HDF5 group or file.

        Parameters
        ----------
        group_or_filename : h5py.Group or str
            HDF5 group containing multipole data. If given as a string, it is
            assumed to be the filename for the HDF5 file, and the first group is
            used to read from.

        Returns
        -------
        openmc.data.WindowedMultipole
            Resonant cross sections represented in the windowed multipole
            format.

        """
        if isinstance(group_or_filename, h5py.Group):
            group = group_or_filename
        else:
            h5file = h5py.File(group_or_filename, 'r')
            try:
                version = h5file['version'].value.decode()
            except AttributeError:
                version = h5file['version'].value[0].decode()
            if version != WMP_VERSION:
                raise ValueError('The given WMP data uses version '
                    + version + ' whereas your installation of the OpenMC '
                    'Python API expects version ' + WMP_VERSION)
            group = h5file['nuclide']

        out = cls()

        # Read scalar values.  Note that group['max_w'] is ignored.

        length = group['length'].value
        windows = group['windows'].value
        out.num_l = group['num_l'].value
        out.fit_order = group['fit_order'].value
        out.fissionable = bool(group['fissionable'].value)

        if group['formalism'].value == _FORM_MLBW:
            out.formalism = 'MLBW'
        elif group['formalism'].value == _FORM_RM:
            out.formalism = 'RM'
        else:
            raise ValueError('Unrecognized/Unsupported R-matrix formalism')

        out.spacing = group['spacing'].value
        out.sqrtAWR = group['sqrtAWR'].value
        out.start_E = group['start_E'].value
        out.end_E = group['end_E'].value

        # Read arrays.

        err = "WMP '{}' array shape is not consistent with the '{}' value"

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

        # Note that all the file 3 data (group['reactions/MT...']) are ignored.

        return out

    def _evaluate(self, E, T):
        """Compute total, absorption, and fission cross sections.

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

        if E < self.start_E: return (0, 0, 0)
        if E > self.end_E: return (0, 0, 0)

        # ======================================================================
        # Bookkeeping

        # Define some frequently used variables.
        sqrtkT = sqrt(K_BOLTZMANN * T)
        sqrtE = sqrt(E)
        invE = 1.0 / E

        # Locate us.  The i_window calc omits a + 1 present in F90 because of
        # the 1-based vs. 0-based indexing.  Similarly startw needs to be
        # decreased by 1.  endw does not need to be decreased because
        # range(startw, endw) does not include endw.
        i_window = int(np.floor((sqrtE - sqrt(self.start_E)) / self.spacing))
        startw = self.w_start[i_window] - 1
        endw = self.w_end[i_window]

        # Fill in factors.  Because of the unique interference dips in scatering
        # resonances, the total cross section has a special "factor" that does
        # not appear in the absorption and fission equations.
        if startw <= endw:
            twophi = np.zeros(self.num_l, dtype=np.float)
            sig_t_factor = np.zeros(self.num_l, dtype=np.cfloat)

            for iL in range(self.num_l):
                twophi[iL] = self.pseudo_k0RS[iL] * sqrtE
                if iL == 1:
                    twophi[iL] = twophi[iL] - np.arctan(twophi[iL])
                elif iL == 2:
                    arg = 3.0 * twophi[iL] / (3.0 - twophi[iL]**2)
                    twophi[iL] = twophi[iL] - np.arctan(arg)
                elif iL == 3:
                    arg = (twophi[iL] * (15.0 - twophi[iL]**2)
                           / (15.0 - 6.0 * twophi[iL]**2))
                    twophi[iL] = twophi[iL] - np.arctan(arg)

            twophi = 2.0 * twophi
            sig_t_factor = np.cos(twophi) - 1j*np.sin(twophi)

        # Initialize the ouptut cross sections.
        sig_t = 0.0
        sig_a = 0.0
        sig_f = 0.0

        # ======================================================================
        # Add the contribution from the curvefit polynomial.

        if sqrtkT != 0 and self.broaden_poly[i_window]:
            # Broaden the curvefit.
            dopp = self.sqrtAWR / sqrtkT
            broadened_polynomials = _broaden_wmp_polynomials(E, dopp,
                                                             self.fit_order + 1)
            for i_poly in range(self.fit_order+1):
                sig_t += (self.curvefit[i_window, i_poly, _FIT_T]
                          * broadened_polynomials[i_poly])
                sig_a += (self.curvefit[i_window, i_poly, _FIT_A]
                          * broadened_polynomials[i_poly])
                if self.fissionable:
                    sig_f += (self.curvefit[i_window, i_poly, _FIT_F]
                              * broadened_polynomials[i_poly])
        else:
            temp = invE
            for i_poly in range(self.fit_order+1):
                sig_t += self.curvefit[i_window, i_poly, _FIT_T] * temp
                sig_a += self.curvefit[i_window, i_poly, _FIT_A] * temp
                if self.fissionable:
                    sig_f += self.curvefit[i_window, i_poly, _FIT_F] * temp
                temp *= sqrtE

        # ======================================================================
        # Add the contribution from the poles in this window.

        if sqrtkT == 0.0:
            # If at 0K, use asymptotic form.
            for i_pole in range(startw, endw):
                psi_chi = -1j / (self.data[i_pole, _MP_EA] - sqrtE)
                c_temp = psi_chi / E
                if self.formalism == 'MLBW':
                    sig_t += ((self.data[i_pole, _MLBW_RT] * c_temp *
                               sig_t_factor[self.l_value[i_pole]-1]).real
                              + (self.data[i_pole, _MLBW_RX] * c_temp).real)
                    sig_a += (self.data[i_pole, _MLBW_RA] * c_temp).real
                    if self.fissionable:
                        sig_f += (self.data[i_pole, _MLBW_RF] * c_temp).real
                elif self.formalism == 'RM':
                    sig_t += (self.data[i_pole, _RM_RT] * c_temp *
                              sig_t_factor[self.l_value[i_pole]-1]).real
                    sig_a += (self.data[i_pole, _RM_RA] * c_temp).real
                    if self.fissionable:
                        sig_f += (self.data[i_pole, _RM_RF] * c_temp).real
                else:
                    raise ValueError('Unrecognized/Unsupported R-matrix'
                                     ' formalism')

        else:
            # At temperature, use Faddeeva function-based form.
            dopp = self.sqrtAWR / sqrtkT
            for i_pole in range(startw, endw):
                Z = (sqrtE - self.data[i_pole, _MP_EA]) * dopp
                w_val = _faddeeva(Z) * dopp * invE * sqrt(pi)
                if self.formalism == 'MLBW':
                    sig_t += ((self.data[i_pole, _MLBW_RT] *
                               sig_t_factor[self.l_value[i_pole]-1] +
                               self.data[i_pole, _MLBW_RX]) * w_val).real
                    sig_a += (self.data[i_pole, _MLBW_RA] * w_val).real
                    if self.fissionable:
                        sig_f += (self.data[i_pole, _MLBW_RF] * w_val).real
                elif self.formalism == 'RM':
                    sig_t += (self.data[i_pole, _RM_RT] * w_val *
                              sig_t_factor[self.l_value[i_pole]-1]).real
                    sig_a += (self.data[i_pole, _RM_RA] * w_val).real
                    if self.fissionable:
                        sig_f += (self.data[i_pole, _RM_RF] * w_val).real
                else:
                    raise ValueError('Unrecognized/Unsupported R-matrix'
                                     ' formalism')

        return sig_t, sig_a, sig_f

    def __call__(self, E, T):
        """Compute total, absorption, and fission cross sections.

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
