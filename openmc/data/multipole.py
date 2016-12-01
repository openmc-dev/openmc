from six import string_types

import h5py
import numpy as np
from scipy.special import wofz

from . import WMP_VERSION
from .data import K_BOLTZMANN
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


def faddeeva(z):
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

    """
    if np.angle(z) > 0:
        return wofz(z)
    else:
        return -np.conj(wofz(z))


def broaden_wmp_polynomials(En, dopp, n):
    sqrtE = sqrt(En)
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

    factors[0] = erfbeta / En
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

        # Scalars.
        out.length = group['length'].value
        out.windows = group['windows'].value
        out.num_l = group['num_l'].value
        out.fit_order = group['fit_order'].value
        out.max_w = group['max_w'].value
        out.fissionable = group['fissionable'].value
        out.formalism = group['formalism'].value
        out.spacing = group['spacing'].value
        out.sqrtAWR = group['sqrtAWR'].value
        out.start_E = group['start_E'].value
        out.end_E = group['end_E'].value

        # Arrays.
        out.data = group['data'].value
        out.pseudo_k0RS = group['pseudo_K0RS'].value
        out.l_value = group['l_value'].value
        out.w_start = group['w_start'].value
        out.w_end = group['w_end'].value
        out.broaden_poly = group['broaden_poly'].value
        out.curvefit = group['curvefit'].value

        return out

    def evaluate(self, E, T):
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
            broadened_polynomials = broaden_wmp_polynomials(E, dopp,
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
            if self.formalism == FORM_MLBW:
                sigT += ((self.data[i_pole, MLBW_RT] * c_temp *
                          sigT_factor[self.l_value[i_pole]-1]).real
                         + (self.data[i_pole, MLBW_RX] * c_temp).real)
                sigA += (self.data[i_pole, MLBW_RA] * c_temp).real
                sigF += (self.data[i_pole, MLBW_RF] * c_temp).real
            elif self.formalism == FORM_RM:
                sigT += (self.data[i_pole, RM_RT] * c_temp *
                         sigT_factor[self.l_value[i_pole]-1]).real
                sigA += (self.data[i_pole, RM_RA] * c_temp).real
                sigF += (self.data[i_pole, RM_RF] * c_temp).real
        else:
          # At temperature, use Faddeeva function-based form.
          for i_pole in range(startw, endw):
              Z = (sqrtE - self.data[i_pole, MP_EA]) * dopp
              w_val = faddeeva(Z) * dopp * invE * np.sqrt(np.pi)
              if self.formalism == FORM_MLBW:
                  sigT += ((self.data[i_pole, MLBW_RT] *
                            sigT_factor[self.l_value[i_pole]-1] +
                            self.data[i_pole, MLBW_RX]) * w_val).real
                  sigA += (self.data[i_pole, MLBW_RA] * w_val).real
                  sigF += (self.data[i_pole, MLBW_RF] * w_val).real
              elif self.formalism == FORM_RM:
                  sigT += (self.data[i_pole, RM_RT] * w_val *
                           sigT_factor[self.l_value[i_pole]-1]).real
                  sigA += (self.data[i_pole, RM_RA] * w_val).real
                  sigF += (self.data[i_pole, RM_RF] * w_val).real

        return sigT, sigA, sigF
