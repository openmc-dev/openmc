from numbers import Integral, Real
from math import exp, erf, pi, sqrt

import h5py
import numpy as np

from . import WMP_VERSION, WMP_VERSION_MAJOR
from .data import K_BOLTZMANN
import openmc.checkvalue as cv
from openmc.mixin import EqualityMixin


# Constants that determine which value to access
_MP_EA = 0       # Pole

# Residue indices
_MP_RS = 1       # Residue scattering
_MP_RA = 2       # Residue absorption
_MP_RF = 3       # Residue fission

# Polynomial fit indices
_FIT_S = 0       # Scattering
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

    Parameters
    ----------
    name : str
        Name of the nuclide using the GND naming convention

    Attributes
    ----------
    fit_order : Integral
        Order of the windowed curvefit.
    fissionable : bool
        Whether or not the target nuclide has fission data.
    spacing : Real
        The width of each window in sqrt(E)-space.  For example, the frst window
        will end at (sqrt(E_min) + spacing)**2 and the second window at
        (sqrt(E_min) + 2*spacing)**2.
    sqrtAWR : Real
        Square root of the atomic weight ratio of the target nuclide.
    E_min : Real
        Lowest energy in eV the library is valid for.
    E_max : Real
        Highest energy in eV the library is valid for.
    data : np.ndarray
        A 2D array of complex poles and residues.  data[i, 0] gives the energy
        at which pole i is located.  data[i, 1:] gives the residues associated
        with the i-th pole.  There are 3 residues, one each for the scattering,
        absorption, and fission channels.
    windows : np.ndarray
        A 2D array of Integral values.  windows[i, 0] - 1 is the index of the
        first pole in window i. windows[i, 1] - 1 is the index of the last pole
        in window i.
    broaden_poly : np.ndarray
        A 1D array of boolean values indicating whether or not the polynomial
        curvefit in that window should be Doppler broadened.
    curvefit : np.ndarray
        A 3D array of Real curvefit polynomial coefficients.  curvefit[i, 0, :]
        gives coefficients for the scattering cross section in window i.
        curvefit[i, 1, :] gives absorption coefficients and curvefit[i, 2, :]
        gives fission coefficients.  The polynomial terms are increasing powers
        of sqrt(E) starting with 1/E e.g:
        a/E + b/sqrt(E) + c + d sqrt(E) + ...

    """
    def __init__(self, name):
        self.name = name
        self.spacing = None
        self.sqrtAWR = None
        self.E_min = None
        self.E_max = None
        self.data = None
        self.windows = None
        self.broaden_poly = None
        self.curvefit = None

    @property
    def name(self):
        return self._name

    @property
    def fit_order(self):
        return self.curvefit.shape[1] - 1

    @property
    def fissionable(self):
        return self.data.shape[1] == 4

    @property
    def spacing(self):
        return self._spacing

    @property
    def sqrtAWR(self):
        return self._sqrtAWR

    @property
    def E_min(self):
        return self._E_min

    @property
    def E_max(self):
        return self._E_max

    @property
    def data(self):
        return self._data

    @property
    def windows(self):
        return self._windows

    @property
    def broaden_poly(self):
        return self._broaden_poly

    @property
    def curvefit(self):
        return self._curvefit

    @name.setter
    def name(self, name):
        cv.check_type('name', name, str)
        self._name = name

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

    @E_min.setter
    def E_min(self, E_min):
        if E_min is not None:
            cv.check_type('E_min', E_min, Real)
            cv.check_greater_than('E_min', E_min, 0.0, equality=True)
        self._E_min = E_min

    @E_max.setter
    def E_max(self, E_max):
        if E_max is not None:
            cv.check_type('E_max', E_max, Real)
            cv.check_greater_than('E_max', E_max, 0.0, equality=False)
        self._E_max = E_max

    @data.setter
    def data(self, data):
        if data is not None:
            cv.check_type('data', data, np.ndarray)
            if len(data.shape) != 2:
                raise ValueError('Multipole data arrays must be 2D')
            if data.shape[1] not in (3, 4):
                raise ValueError(
                     'data.shape[1] must be 3 or 4. One value for the pole.'
                     ' One each for the scattering and absorption residues. '
                     'Possibly one more for a fission residue.')
            if not np.issubdtype(data.dtype, np.complexfloating):
                raise TypeError('Multipole data arrays must be complex dtype')
        self._data = data

    @windows.setter
    def windows(self, windows):
        if windows is not None:
            cv.check_type('windows', windows, np.ndarray)
            if len(windows.shape) != 2:
                raise ValueError('Multipole windows arrays must be 2D')
            if not np.issubdtype(windows.dtype, np.integer):
                raise TypeError('Multipole windows arrays must be integer'
                                ' dtype')
        self._windows = windows

    @broaden_poly.setter
    def broaden_poly(self, broaden_poly):
        if broaden_poly is not None:
            cv.check_type('broaden_poly', broaden_poly, np.ndarray)
            if len(broaden_poly.shape) != 1:
                raise ValueError('Multipole broaden_poly arrays must be 1D')
            if not np.issubdtype(broaden_poly.dtype, np.bool_):
                raise TypeError('Multipole broaden_poly arrays must be boolean'
                                ' dtype')
        self._broaden_poly = broaden_poly

    @curvefit.setter
    def curvefit(self, curvefit):
        if curvefit is not None:
            cv.check_type('curvefit', curvefit, np.ndarray)
            if len(curvefit.shape) != 3:
                raise ValueError('Multipole curvefit arrays must be 3D')
            if curvefit.shape[2] not in (2, 3):  # sig_s, sig_a (maybe sig_f)
                raise ValueError('The third dimension of multipole curvefit'
                                 ' arrays must have a length of 2 or 3')
            if not np.issubdtype(curvefit.dtype, np.floating):
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
            h5file = h5py.File(str(group_or_filename), 'r')

            # Make sure version matches
            if 'version' in h5file.attrs:
                major, minor = h5file.attrs['version']
                if major != WMP_VERSION_MAJOR:
                    raise IOError(
                        'WMP data format uses version {}. {} whereas your '
                        'installation of the OpenMC Python API expects version '
                        '{}.x.'.format(major, minor, WMP_VERSION_MAJOR))
            else:
                raise IOError(
                    'WMP data does not indicate a version. Your installation of '
                    'the OpenMC Python API expects version {}.x data.'
                    .format(WMP_VERSION_MAJOR))

            group = list(h5file.values())[0]

        name = group.name[1:]
        out = cls(name)

        # Read scalars.

        out.spacing = group['spacing'][()]
        out.sqrtAWR = group['sqrtAWR'][()]
        out.E_min = group['E_min'][()]
        out.E_max = group['E_max'][()]

        # Read arrays.

        err = "WMP '{}' array shape is not consistent with the '{}' array shape"

        out.data = group['data'][()]

        out.windows = group['windows'][()]

        out.broaden_poly = group['broaden_poly'][...].astype(np.bool)
        if out.broaden_poly.shape[0] != out.windows.shape[0]:
            raise ValueError(err.format('broaden_poly', 'windows'))

        out.curvefit = group['curvefit'][()]
        if out.curvefit.shape[0] != out.windows.shape[0]:
            raise ValueError(err.format('curvefit', 'windows'))

        # _broaden_wmp_polynomials assumes the curve fit has at least 3 terms.
        if out.fit_order < 2:
            raise ValueError("Windowed multipole is only supported for "
                             "curvefits with 3 or more terms.")

        return out

    def _evaluate(self, E, T):
        """Compute scattering, absorption, and fission cross sections.

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

        if E < self.E_min: return (0, 0, 0)
        if E > self.E_max: return (0, 0, 0)

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
        i_window = int(np.floor((sqrtE - sqrt(self.E_min)) / self.spacing))
        startw = self.windows[i_window, 0] - 1
        endw = self.windows[i_window, 1]

        # Initialize the ouptut cross sections.
        sig_s = 0.0
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
                sig_s += (self.curvefit[i_window, i_poly, _FIT_S]
                          * broadened_polynomials[i_poly])
                sig_a += (self.curvefit[i_window, i_poly, _FIT_A]
                          * broadened_polynomials[i_poly])
                if self.fissionable:
                    sig_f += (self.curvefit[i_window, i_poly, _FIT_F]
                              * broadened_polynomials[i_poly])
        else:
            temp = invE
            for i_poly in range(self.fit_order+1):
                sig_s += self.curvefit[i_window, i_poly, _FIT_S] * temp
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
                sig_s += (self.data[i_pole, _MP_RS] * c_temp).real
                sig_a += (self.data[i_pole, _MP_RA] * c_temp).real
                if self.fissionable:
                    sig_f += (self.data[i_pole, _MP_RF] * c_temp).real

        else:
            # At temperature, use Faddeeva function-based form.
            dopp = self.sqrtAWR / sqrtkT
            for i_pole in range(startw, endw):
                Z = (sqrtE - self.data[i_pole, _MP_EA]) * dopp
                w_val = _faddeeva(Z) * dopp * invE * sqrt(pi)
                sig_s += (self.data[i_pole, _MP_RS] * w_val).real
                sig_a += (self.data[i_pole, _MP_RA] * w_val).real
                if self.fissionable:
                    sig_f += (self.data[i_pole, _MP_RF] * w_val).real

        return sig_s, sig_a, sig_f

    def __call__(self, E, T):
        """Compute scattering, absorption, and fission cross sections.

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

    def export_to_hdf5(self, path, mode='a', libver='earliest'):
        """Export windowed multipole data to an HDF5 file.

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

        # Open file and write version.
        with h5py.File(str(path), mode, libver=libver) as f:
            f.attrs['filetype'] = np.string_('data_wmp')
            f.attrs['version'] = np.array(WMP_VERSION)

            g = f.create_group(self.name)

            # Write scalars.
            g.create_dataset('spacing', data=np.array(self.spacing))
            g.create_dataset('sqrtAWR', data=np.array(self.sqrtAWR))
            g.create_dataset('E_min', data=np.array(self.E_min))
            g.create_dataset('E_max', data=np.array(self.E_max))

            # Write arrays.
            g.create_dataset('data', data=self.data)
            g.create_dataset('windows', data=self.windows)
            g.create_dataset('broaden_poly',
                             data=self.broaden_poly.astype(np.int8))
            g.create_dataset('curvefit', data=self.curvefit)
