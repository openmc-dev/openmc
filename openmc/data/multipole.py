from numbers import Real
from math import exp, erf, pi, sqrt

import os
import h5py
import numpy as np
from scipy.signal import find_peaks

import openmc.checkvalue as cv
from ..exceptions import DataError
from ..mixin import EqualityMixin
from . import WMP_VERSION, WMP_VERSION_MAJOR
from .data import K_BOLTZMANN
from .neutron import IncidentNeutron
from .resonance import ResonanceRange

import vectfit as m

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

# Upper temperature limit
TEMPERATURE_LIMIT = 3000

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
    np.ndarray
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

def _vectfit_xs(energy, ce_xs, mts, rtol=1e-3, atol=1e-5, orders=None,
                n_iter_vf=30, log=False, path_plot=None):
    r"""Generate multipole data from point-wise cross sections.

    Parameters
    ----------
    energy : np.ndarray
        Energy array.
    cs_xs : np.ndarray
        Point-wise cross sections to be fitted.
    mts : Iterable of Integral
        Reaction list.
    rtol : Real, optional
        Relative error tolerance.
    atol : Real, optional
        Absolute error tolerance.
    orders : Iterable of Integral, optional
        A list of orders (number of poles) to be searched.
    n_iter_vf : Integral, optional
        Number of maximum VF iterations.
    log : Bool, optional
        Whether to print log.
    path_plot : str, optional
        Path to save the figures.

    Returns
    -------
    Tuple
        (poles, residues).

    """

    MIN_CROSS_SECTION = 1e-7
    N_FINER = 10

    ne = energy.size
    nmt = len(mts)
    if ce_xs.shape != (nmt, ne):
        raise ValueError('Inconsistent cross section data.')

    # construct test data: interpolate xs with finer grids
    ne_test = (ne-1)*N_FINER + 1
    test_energy = np.interp(np.arange(ne_test),
                            np.arange(ne_test, step=N_FINER), energy)
    test_energy[[0, -1]] = energy[[0, -1]] # avoid numerical issue
    test_xs_ref = np.zeros((nmt, ne_test))
    for i in range(nmt):
        test_xs_ref[i] = np.interp(test_energy, energy, ce_xs[i])

    if log:
        print("Energy: {} to {} ({} points)".format(energy[0], energy[-1], ne))
    # inputs
    f = ce_xs * energy # sigma*E
    s = np.sqrt(energy) # sqrt(E)
    test_s = np.sqrt(test_energy)
    weight = 1.0/f
    # very small cross sections can lead to huge weights which harm accuracy
    for i in range(nmt):
        if np.all(ce_xs[i]<=MIN_CROSS_SECTION):
            weight[i] = 1.0
        elif np.any(ce_xs[i]<=MIN_CROSS_SECTION):
            weight[i, ce_xs[i]<=MIN_CROSS_SECTION] = \
               max(weight[i, ce_xs[i]>MIN_CROSS_SECTION])

    # order search
    peaks, _ = find_peaks(ce_xs[0]+ce_xs[1])
    n_peaks = peaks.size
    if orders is not None:
        # make sure orders are even integers
        orders = list(set([int(i/2)*2 for i in orders if i>=2]))
    else:
        lowest_order = max(2, 2*n_peaks)
        highest_order = max(200, 4*n_peaks)
        orders = list(range(lowest_order, highest_order+1, 2))

    if log:
        print("Found {} peaks".format(n_peaks))
        print("Fitting orders from {} to {}".format(orders[0], orders[-1]))

    found_ideal = False
    n_discarded = 0 # for accelation, number of discarded searches
    best_quality = best_ratio = -np.inf
    for i, order in enumerate(orders):
        if log:
            print("Order={} {}/{}".format(order, i, len(orders)))
        # initial guessed poles
        poles = np.linspace(s[0], s[-1], order/2)
        poles = poles + poles*0.01j
        poles = np.sort(np.append(poles, np.conj(poles)))

        found_better = False
        # fitting iteration
        for i_vf in range(n_iter_vf):
            if log:
                print("VF iteration {}/{}".format(i_vf+1, n_iter_vf))

            # call vf
            poles, residues, cf, f_fit, rms = m.vectfit(f, s, poles, weight)

            # convert real pole to conjugate pairs
            n_real_poles = 0
            new_poles = []
            for p in poles:
                p_r, p_i = np.real(p), np.imag(p)
                if p_r > s[0] and p_r < s[-1] and p_i == 0.:
                    new_poles += [p_r+p_r*0.01j, p_r-p_r*0.01j]
                    n_real_poles += 1
                else:
                    new_poles += [p]
            new_poles = np.array(new_poles)
            # re-calculate residues if poles changed
            if n_real_poles > 0:
                new_poles, residues, cf, f_fit, rms = \
                      m.vectfit(f, s, new_poles, weight, skip_pole=True)

            # assess the result on test grid
            test_xs = m.evaluate(test_s, new_poles, residues)/test_energy
            abserr = np.abs(test_xs - test_xs_ref)
            relerr = abserr/test_xs_ref
            if np.any(np.isnan(abserr)):
                maxre, ratio, ratio2 = np.inf, -np.inf, -np.inf
            elif np.all(abserr <= atol):
                maxre, ratio, ratio2 = 0., 1., 1.
            else:
                maxre = np.max(relerr[abserr > atol])
                ratio = np.sum((relerr<rtol) | (abserr<atol)) / relerr.size
                ratio2 = np.sum((relerr<10*rtol) | (abserr<atol)) / relerr.size

            quality = 100*(ratio+ratio2-min(0.1*maxre, 1)-5e-4*new_poles.size)

            if np.any(test_xs < -atol):
                quality = -np.inf

            if log:
                print("  Max relative error: {:.3f}%".format(maxre*100))
                print("  Satisfaction: {:.1f}%, {:.1f}%".format(ratio*100, ratio2*100))
                print("  Quality: {:.2f}".format(quality))

            if quality > best_quality:
                if log:
                    print("  Best by far!")
                found_better = True
                best_quality, best_ratio = quality, ratio
                best_poles, best_residues = new_poles, residues
                best_test_xs, best_relerr = test_xs, relerr
                if ratio >= 1.0:
                    if log:
                        print("Found ideal results. Stop!")
                    found_ideal = True
                    break
            else:
                if log:
                    print("  Discarded!")

        if found_ideal:
            break

        # acceleration
        if found_better:
            n_discarded = 0
        else:
            if order > max(2*n_peaks, 50) and best_ratio > 0.7:
                n_discarded += 1
                if n_discarded >= 10 or (n_discarded >= 5 and best_ratio > 0.9):
                    if log:
                        print("Couldn't get better results. Stop!")
                    break

    # merge conjugate poles
    real_idx = conj_idx = []
    found_conj = False
    for i, p in enumerate(best_poles):
        if found_conj:
            found_conj = False
            continue
        if np.imag(p) == 0.:
            real_idx.append(i)
        else:
            if i < best_poles.size and np.conj(p) == best_poles[i+1]:
                found_conj = True
                conj_idx.append(i)
            else:
                raise RuntimeError("Complex poles are not conjugate!")
    if log:
        print("Found {} real poles and {} conjugate complex pairs.".format(
               len(real_idx), len(conj_idx)))
    mp_poles = best_poles[real_idx+conj_idx]
    mp_residues = np.concatenate((best_residues[:, real_idx],
                                  best_residues[:, conj_idx]*2), axis=1)/1j
    if log:
        print("Final number of poles: {}".format(mp_poles.size))

    if path_plot:
        import matplotlib
        matplotlib.use("agg")
        import matplotlib.pyplot as plt
        if not os.path.exists(path_plot):
            os.makedirs(path_plot)
        for i, mt in enumerate(mts):
            fig, ax1 = plt.subplots()
            lns1 = ax1.loglog(test_energy, test_xs_ref[i], 'g', label="ACE xs")
            lns2 = ax1.loglog(test_energy, best_test_xs[i], 'b', label="VF xs")
            ax2 = ax1.twinx()
            lns3 = ax2.loglog(test_energy, best_relerr[i], 'r',
                              label="Relative error", alpha=0.5)
            lns = lns1 + lns2 + lns3
            labels = [l.get_label() for l in lns]
            ax1.legend(lns, labels, loc='best')
            ax1.set_xlabel('energy (eV)')
            ax1.set_ylabel('cross section (b)', color='b')
            ax1.tick_params('y', colors='b')
            ax2.set_ylabel('relative error', color='r')
            ax2.tick_params('y', colors='r')

            plt.title("MT {} vectfitted with {} poles".format(mt, mp_poles.size))
            fig.tight_layout()
            figfile = os.path.join(path_plot, "{}_vf.png".format(mt))
            plt.savefig(figfile)
            plt.close()
            if log:
                print("Plot figure: {}".format(figfile))

    return (mp_poles, mp_residues)

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
            need_to_close = False
        else:
            h5file = h5py.File(str(group_or_filename), 'r')
            need_to_close = True

            # Make sure version matches
            if 'version' in h5file.attrs:
                major, minor = h5file.attrs['version']
                if major != WMP_VERSION_MAJOR:
                    raise DataError(
                        'WMP data format uses version {}. {} whereas your '
                        'installation of the OpenMC Python API expects version '
                        '{}.x.'.format(major, minor, WMP_VERSION_MAJOR))
            else:
                raise DataError(
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

        # If HDF5 file was opened here, make sure it gets closed
        if need_to_close:
            h5file.close()

        return out

    @classmethod
    def from_vectfit(cls, endf_file, error=1e-3, njoy_error=5e-4, log=False,
                     path=None, **kwargs):
        """Generate windowed multipole neutron data via Vector Fitting.

        Parameters
        ----------
        endf_file : str
            Path to ENDF evaluation
        error : float, optional
            Fractional error tolerance for data fitting
        njoy_error : float, optional
            Fractional error tolerance for processing point-wise data with NJOY
        log : bool, optional
            Whether to display log
        path : str, optional
            Path to save the figures.
        **kwargs
            Keyword arguments passed to :func:`openmc.data.multipole._vectfit_xs`

        Returns
        -------
        openmc.data.WindowedMultipole
            Resonant cross sections represented in the windowed multipole
            format.

        """

        # ======================================================================
        # PREPARE POINT-WISE XS

        # make 0K ACE data using njoy
        if log:
            print("Running NJOY to get 0K point-wise data...")
        nuc_ce = IncidentNeutron.from_njoy(endf_file, temperatures=[0.0],
                        error=njoy_error, broadr=False, heatr=False, purr=False)

        if log:
            print("Parsing cross section within resolved resonance range...")
        # RRR bound
        endf_res = IncidentNeutron.from_endf(endf_file).resonances

        rrr_bound_energy = nuc_ce.energy['0K'][-1]
        try:
            rrr = endf_res.resolved
        except:
            rrr = None
        # if resolved resonance parameters exist
        if rrr is not None and hasattr(rrr, 'energy_max') and \
             type(rrr) is not ResonanceRange:
            rrr_bound_energy = rrr.energy_max
        else:
            try:
                # set rrr bound as lower bound of unresolved
                rrr_bound_energy = endf_res.unresolved.energy_min
            except:
                pass
        rrr_bound_idx = np.searchsorted(nuc_ce.energy['0K'], rrr_bound_energy,
                                        side='right') - 1

        # first threshold
        first_threshold_idx = float("inf")
        first_threshold_mt = None
        for mt in nuc_ce.reactions:
            if hasattr(nuc_ce.reactions[mt].xs['0K'], '_threshold_idx'):
                threshold_idx = nuc_ce.reactions[mt].xs['0K']._threshold_idx
                if 0 < threshold_idx < first_threshold_idx:
                    first_threshold_idx = threshold_idx
                    first_threshold_mt = mt

        # lower of RRR bound and first threshold
        e_max_idx = min(rrr_bound_idx, first_threshold_idx)

        if log:
            print("RRR idx: {}, first threshold idx: {}".format(rrr_bound_idx,
                  first_threshold_idx))

        # parse energy and summed cross sections
        energy = nuc_ce.energy['0K'][:e_max_idx+1]
        E_min, E_max = energy[0], energy[-1]
        n_points = energy.size

        total_xs = nuc_ce[1].xs['0K'](energy)
        if 2 in nuc_ce:
            elastic_xs = nuc_ce[2].xs['0K'](energy)
        else:
            elastic_xs = np.zeros_like(total_xs)
        if 27 in nuc_ce:
            absorption_xs = nuc_ce[27].xs['0K'](energy)
        else:
            absorption_xs = np.zeros_like(total_xs)
        fissionable = False
        if 18 in nuc_ce:
            fission_xs = nuc_ce[18].xs['0K'](energy)
            fissionable = True

        # make vectors
        if fissionable:
            ce_xs = np.vstack((elastic_xs, absorption_xs, fission_xs))
            mts = [2, 27, 18]
        else:
            ce_xs = np.vstack((elastic_xs, absorption_xs))
            mts = [2, 27]

        if log:
            print("  MTs: {}, Energy range: {:e} to {:e} eV ({} points)".format(
                  mts, E_min, E_max, n_points))

        # ======================================================================
        # PERFORM VECTOR FITTING

        alpha = nuc_ce.atomic_weight_ratio/(K_BOLTZMANN*TEMPERATURE_LIMIT)

        # divide into pieces for complex nuclides
        peaks, _ = find_peaks(total_xs)
        n_peaks = peaks.size
        if n_peaks > 300 or n_points > 50000 or n_peaks * n_points > 100*30000:
            n_pieces = max(5, n_peaks // 80,  n_points // 5000)
        else:
            n_pieces = 1
        piece_width = (sqrt(E_max) - sqrt(E_min)) / n_pieces

        # VF piece by piece
        for i_piece in range(n_pieces):
            if log:
                print("Piece {}/{}".format(i_piece+1, n_pieces))
            # start E of this piece
            e_bound = (sqrt(E_min) + piece_width*(i_piece-0.5))**2
            if i_piece == 0 or sqrt(alpha*e_bound) < 4.0:
                e_start = E_min
                e_start_idx = 0
            else:
                e_start = max(E_min, (sqrt(alpha*e_bound)-4.0)**2/alpha)
                e_start_idx = np.searchsorted(energy, e_start, side='right') - 1
            # end E of this piece
            e_bound = (sqrt(E_min) + piece_width*(i_piece+1))**2
            e_end = min(E_max, (sqrt(alpha*e_bound) + 4.0)**2/alpha)
            e_end_idx = np.searchsorted(energy, e_end, side='left') + 1
            piece_range = range(e_start_idx, min(e_end_idx+1, n_points))

            # fitting xs
            if path:
                path_plot = os.path.join(path, "{}".format(i_piece+1))
            else:
                path_plot = None
            poles, residues = _vectfit_xs(energy[piece_range],
                   ce_xs[:, piece_range], mts, rtol=error, log=log,
                   path_plot=path_plot, **kwargs)

        # ======================================================================
        # WINDOWING

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
        mode : {'r+', 'w', 'x', 'a'}
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
