from numbers import Real
from math import exp, erf, pi, sqrt
from copy import deepcopy

import os
import h5py
import pickle
import numpy as np
from scipy.signal import find_peaks

import openmc.checkvalue as cv

from ..exceptions import DataError
from ..mixin import EqualityMixin
from . import WMP_VERSION, WMP_VERSION_MAJOR
from .data import K_BOLTZMANN
from .neutron import IncidentNeutron
from .resonance import ResonanceRange
from .vectfit import vectfit

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

# Upper temperature limit (K)
TEMPERATURE_LIMIT = 3000

# Logging control
DETAILED_LOGGING = 2


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
    E : float
        Energy to evaluate at.
    dopp : float
        sqrt(atomic weight ratio / kT) in units of eV.
    n : int
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
                n_vf_iter=30, log=False, path_out=None):
    """Convert point-wise cross section to multipole data via vector fitting.

    Parameters
    ----------
    energy : np.ndarray
        Energy array
    ce_xs : np.ndarray
        Point-wise cross sections to be fitted, with shape (number of reactions,
        number of energy points)
    mts : Iterable of int
        Reaction list
    rtol : float, optional
        Relative error tolerance
    atol : float, optional
        Absolute error tolerance
    orders : Iterable of int, optional
        A list of orders (number of poles) to be searched
    n_vf_iter : int, optional
        Number of maximum VF iterations
    log : bool or int, optional
        Whether to print running logs (use int for verbosity control)
    path_out : str, optional
        Path to save the figures to show discrepancies between the original and
        fitted cross sections for different reactions

    Returns
    -------
    tuple
        (poles, residues)

    """
    ne = energy.size
    nmt = len(mts)
    if ce_xs.shape != (nmt, ne):
        raise ValueError('Inconsistent cross section data.')

    # construct test data: interpolate xs with finer grids
    n_finer = 10
    ne_test = (ne - 1)*n_finer + 1
    test_energy = np.interp(np.arange(ne_test),
                            np.arange(ne_test, step=n_finer), energy)
    test_energy[[0, -1]] = energy[[0, -1]]  # avoid numerical issue
    test_xs_ref = np.zeros((nmt, ne_test))
    for i in range(nmt):
        test_xs_ref[i] = np.interp(test_energy, energy, ce_xs[i])

    if log:
        print(f"  energy: {energy[0]:.3e} to {energy[-1]:.3e} eV ({ne} points)")
        print(f"  error tolerance: rtol={rtol}, atol={atol}")

    # transform xs (sigma) and energy (E) to f (sigma*E) and s (sqrt(E)) to be
    # compatible with the multipole representation
    f = ce_xs * energy
    s = np.sqrt(energy)
    test_s = np.sqrt(test_energy)

    # inverse weighting is used for minimizing the relative deviation instead of
    # absolute deviation in vector fitting
    with np.errstate(divide='ignore'):
        weight = 1.0/f

    # avoid too large weights which will harm the fitting accuracy
    min_cross_section = 1e-7
    for i in range(nmt):
        if np.all(ce_xs[i] <= min_cross_section):
            weight[i] = 1.0
        elif np.any(ce_xs[i] <= min_cross_section):
            weight[i, ce_xs[i] <= min_cross_section] = \
               max(weight[i, ce_xs[i] > min_cross_section])

    # detect peaks (resonances) and determine VF order search range
    peaks, _ = find_peaks(ce_xs[0] + ce_xs[1])
    n_peaks = peaks.size
    if orders is not None:
        # make sure orders are even integers
        orders = list(set([int(i/2)*2 for i in orders if i >= 2]))
    else:
        lowest_order = max(2, 2*n_peaks)
        highest_order = max(200, 4*n_peaks)
        orders = list(range(lowest_order, highest_order + 1, 2))

    if log:
        print(f"Found {n_peaks} peaks")
        print(f"Fitting orders from {orders[0]} to {orders[-1]}")

    # perform VF with increasing orders
    found_ideal = False
    n_discarded = 0  # for accelation, number of discarded searches
    best_quality = best_ratio = -np.inf
    for i, order in enumerate(orders):
        if log:
            print(f"Order={order}({i}/{len(orders)})")
        # initial guessed poles
        poles_r = np.linspace(s[0], s[-1], order//2)
        poles = poles_r + poles_r*0.01j
        poles = np.sort(np.append(poles, np.conj(poles)))

        found_better = False
        # fitting iteration
        for i_vf in range(n_vf_iter):
            if log >= DETAILED_LOGGING:
                print(f"VF iteration {i_vf + 1}/{n_vf_iter}")

            # call vf
            poles, residues, cf, f_fit, rms = vectfit(f, s, poles, weight)

            # convert real pole to conjugate pairs
            n_real_poles = 0
            new_poles = []
            for p in poles:
                p_r, p_i = np.real(p), np.imag(p)
                if (s[0] <= p_r <= s[-1]) and p_i == 0.:
                    new_poles += [p_r+p_r*0.01j, p_r-p_r*0.01j]
                    n_real_poles += 1
                else:
                    new_poles += [p]
            new_poles = np.array(new_poles)
            # re-calculate residues if poles changed
            if n_real_poles > 0:
                if log >= DETAILED_LOGGING:
                    print(f"  # real poles: {n_real_poles}")
                new_poles, residues, cf, f_fit, rms = \
                      vectfit(f, s, new_poles, weight, skip_pole=True)

            # assess the result on test grid
            test_xs = vf.evaluate(test_s, new_poles, residues) / test_energy
            abserr = np.abs(test_xs - test_xs_ref)
            with np.errstate(invalid='ignore', divide='ignore'):
                relerr = abserr / test_xs_ref
                if np.any(np.isnan(abserr)):
                    maxre, ratio, ratio2 = np.inf, -np.inf, -np.inf
                elif np.all(abserr <= atol):
                    maxre, ratio, ratio2 = 0., 1., 1.
                else:
                    maxre = np.max(relerr[abserr > atol])
                    ratio = np.sum((relerr < rtol) | (abserr < atol)) / relerr.size
                    ratio2 = np.sum((relerr < 10*rtol) | (abserr < atol)) / relerr.size

            # define a metric for choosing the best fitting results
            # basically, it is preferred to have more points within accuracy
            # tolerance, smaller maximum deviation and fewer poles
            #TODO: improve the metric with clearer basis
            quality = ratio + ratio2 - min(0.1*maxre, 1) - 0.001*new_poles.size

            if np.any(test_xs < -atol):
                quality = -np.inf

            if log >= DETAILED_LOGGING:
                print(f"  # poles: {new_poles.size}")
                print(f"  Max relative error: {maxre * 100:.3f}%")
                print(f"  Satisfaction: {ratio * 100:.1f}%, {ratio2 * 100:.1f}%")
                print(f"  Quality: {quality:.2f}")

            if quality > best_quality:
                if log >= DETAILED_LOGGING:
                    print("  Best so far!")
                found_better = True
                best_quality, best_ratio = quality, ratio
                best_poles, best_residues = new_poles, residues
                best_test_xs, best_relerr = test_xs, relerr
                if best_ratio >= 1.0:
                    if log:
                        print("Found ideal results. Stop!")
                    found_ideal = True
                    break
            else:
                if log >= DETAILED_LOGGING:
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
                    if log >= DETAILED_LOGGING:
                        print("Couldn't get better results. Stop!")
                    break

    # merge conjugate poles
    real_idx = []
    conj_idx = []
    found_conj = False
    for i, p in enumerate(best_poles):
        if found_conj:
            found_conj = False
            continue
        if np.imag(p) == 0.:
            real_idx.append(i)
        else:
            if i < best_poles.size and np.conj(p) == best_poles[i + 1]:
                found_conj = True
                conj_idx.append(i)
            else:
                raise RuntimeError("Complex poles are not conjugate!")
    if log:
        print("Found {} real poles and {} conjugate complex pairs.".format(
               len(real_idx), len(conj_idx)))
    mp_poles = best_poles[real_idx + conj_idx]
    mp_residues = np.concatenate((best_residues[:, real_idx],
                                  best_residues[:, conj_idx]*2), axis=1)/1j
    if log:
        print(f"Final number of poles: {mp_poles.size}")

    if path_out:
        if not os.path.exists(path_out):
            os.makedirs(path_out)
        for i, mt in enumerate(mts):
            if not test_xs_ref[i].any():
                continue
            import matplotlib.pyplot as plt
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

            plt.title(f"MT {mt} vector fitted with {mp_poles.size} poles")
            fig.tight_layout()
            fig_file = os.path.join(path_out, "{:.0f}-{:.0f}_MT{}.png".format(
                                    energy[0], energy[-1], mt))
            plt.savefig(fig_file)
            plt.close()
            if log:
                print(f"Saved figure: {fig_file}")

    return (mp_poles, mp_residues)


def vectfit_nuclide(endf_file, njoy_error=5e-4, vf_pieces=None,
                    log=False, path_out=None, mp_filename=None, **kwargs):
    r"""Generate multipole data for a nuclide from ENDF.

    Parameters
    ----------
    endf_file : str
        Path to ENDF evaluation
    njoy_error : float, optional
        Fractional error tolerance for processing point-wise data with NJOY
    vf_pieces : integer, optional
        Number of equal-in-momentum spaced energy pieces for data fitting
    log : bool or int, optional
        Whether to print running logs (use int for verbosity control)
    path_out : str, optional
        Path to write out mutipole data file and vector fitting figures
    mp_filename : str, optional
        File name to write out multipole data
    **kwargs
        Keyword arguments passed to :func:`openmc.data.multipole._vectfit_xs`

    Returns
    -------
    mp_data
        Dictionary containing necessary multipole data of the nuclide

    """

    # ======================================================================
    # PREPARE POINT-WISE XS

    # make 0K ACE data using njoy
    if log:
        print(f"Running NJOY to get 0K point-wise data (error={njoy_error})...")

    nuc_ce = IncidentNeutron.from_njoy(endf_file, temperatures=[0.0],
             error=njoy_error, broadr=False, heatr=False, purr=False)

    if log:
        print("Parsing cross sections within resolved resonance range...")

    # Determine upper energy: the lower of RRR upper bound and first threshold
    endf_res = IncidentNeutron.from_endf(endf_file).resonances
    if hasattr(endf_res, 'resolved') and \
       hasattr(endf_res.resolved, 'energy_max') and \
       type(endf_res.resolved) is not ResonanceRange:
        E_max = endf_res.resolved.energy_max
    elif hasattr(endf_res, 'unresolved') and \
         hasattr(endf_res.unresolved, 'energy_min'):
        E_max = endf_res.unresolved.energy_min
    else:
        E_max = nuc_ce.energy['0K'][-1]
    E_max_idx = np.searchsorted(nuc_ce.energy['0K'], E_max, side='right') - 1
    for mt in nuc_ce.reactions:
        if hasattr(nuc_ce.reactions[mt].xs['0K'], '_threshold_idx'):
            threshold_idx = nuc_ce.reactions[mt].xs['0K']._threshold_idx
            if 0 < threshold_idx < E_max_idx:
                E_max_idx = threshold_idx

    # parse energy and cross sections
    energy = nuc_ce.energy['0K'][:E_max_idx + 1]
    E_min, E_max = energy[0], energy[-1]
    n_points = energy.size
    total_xs = nuc_ce[1].xs['0K'](energy)
    elastic_xs = nuc_ce[2].xs['0K'](energy)

    try:
        absorption_xs = nuc_ce[27].xs['0K'](energy)
    except KeyError:
        absorption_xs = np.zeros_like(total_xs)

    fissionable = False
    try:
        fission_xs = nuc_ce[18].xs['0K'](energy)
        fissionable = True
    except KeyError:
        pass

    # make vectors
    if fissionable:
        ce_xs = np.vstack((elastic_xs, absorption_xs, fission_xs))
        mts = [2, 27, 18]
    else:
        ce_xs = np.vstack((elastic_xs, absorption_xs))
        mts = [2, 27]

    if log:
        print(f"  MTs: {mts}")
        print(f"  Energy range: {E_min:.3e} to {E_max:.3e} eV ({n_points} points)")

    # ======================================================================
    # PERFORM VECTOR FITTING

    if vf_pieces is None:
        # divide into pieces for complex nuclides
        peaks, _ = find_peaks(total_xs)
        n_peaks = peaks.size
        if n_peaks > 200 or n_points > 30000 or n_peaks * n_points > 100*10000:
            vf_pieces = max(5, n_peaks // 50,  n_points // 2000)
        else:
            vf_pieces = 1
    piece_width = (sqrt(E_max) - sqrt(E_min)) / vf_pieces

    alpha = nuc_ce.atomic_weight_ratio/(K_BOLTZMANN*TEMPERATURE_LIMIT)

    poles, residues = [], []
    # VF piece by piece
    for i_piece in range(vf_pieces):
        if log:
            print(f"Vector fitting piece {i_piece + 1}/{vf_pieces}...")
        # start E of this piece
        e_bound = (sqrt(E_min) + piece_width*(i_piece-0.5))**2
        if i_piece == 0 or sqrt(alpha*e_bound) < 4.0:
            e_start = E_min
            e_start_idx = 0
        else:
            e_start = max(E_min, (sqrt(alpha*e_bound) - 4.0)**2/alpha)
            e_start_idx = np.searchsorted(energy, e_start, side='right') - 1
        # end E of this piece
        e_bound = (sqrt(E_min) + piece_width*(i_piece + 1))**2
        e_end = min(E_max, (sqrt(alpha*e_bound) + 4.0)**2/alpha)
        e_end_idx = np.searchsorted(energy, e_end, side='left') + 1
        e_idx = range(e_start_idx, min(e_end_idx + 1, n_points))

        p, r = _vectfit_xs(energy[e_idx], ce_xs[:, e_idx], mts, log=log,
                           path_out=path_out, **kwargs)

        poles.append(p)
        residues.append(r)

    # collect multipole data into a dictionary
    mp_data = {"name": nuc_ce.name,
               "AWR": nuc_ce.atomic_weight_ratio,
               "E_min": E_min,
               "E_max": E_max,
               "poles": poles,
               "residues": residues}

    # dump multipole data to file
    if path_out:
        if not os.path.exists(path_out):
            os.makedirs(path_out)
        if not mp_filename:
            mp_filename = f"{nuc_ce.name}_mp.pickle"
        mp_filename = os.path.join(path_out, mp_filename)
        with open(mp_filename, 'wb') as f:
            pickle.dump(mp_data, f)
        if log:
            print(f"Dumped multipole data to file: {mp_filename}")

    return mp_data


def _windowing(mp_data, n_cf, rtol=1e-3, atol=1e-5, n_win=None, spacing=None,
               log=False):
    """Generate windowed multipole library from multipole data with specific
        settings of window size, curve fit order, etc.

    Parameters
    ----------
    mp_data : dict
        Multipole data
    n_cf : int
        Curve fitting order
    rtol : float, optional
        Maximum relative error tolerance
    atol : float, optional
        Minimum absolute error tolerance
    n_win : int, optional
        Number of equal-in-mementum spaced energy windows
    spacing : float, optional
        Inner window spacing (sqrt energy space)
    log : bool or int, optional
        Whether to print running logs (use int for verbosity control)

    Returns
    -------
    openmc.data.WindowedMultipole
        Resonant cross sections represented in the windowed multipole
        format.

    """
    # unpack multipole data
    name = mp_data["name"]
    awr = mp_data["AWR"]
    E_min = mp_data["E_min"]
    E_max = mp_data["E_max"]
    mp_poles = mp_data["poles"]
    mp_residues = mp_data["residues"]

    n_pieces = len(mp_poles)
    piece_width = (sqrt(E_max) - sqrt(E_min)) / n_pieces
    alpha = awr / (K_BOLTZMANN*TEMPERATURE_LIMIT)

    # determine window size
    if n_win is None:
        if spacing is not None:
            # ensure the windows are within the multipole energy range
            n_win = int((sqrt(E_max) - sqrt(E_min)) / spacing)
            E_max = (sqrt(E_min) + n_win*spacing)**2
        else:
            n_win = 1000
    # inner window size
    spacing = (sqrt(E_max) - sqrt(E_min)) / n_win
    # make sure inner window size is smaller than energy piece size
    if spacing > piece_width:
        raise ValueError('Window spacing cannot be larger than piece spacing.')

    if log:
        print("Windowing:")
        print(f"  config: # windows={n_win}, spacing={spacing}, CF order={n_cf}")
        print(f"  error tolerance: rtol={rtol}, atol={atol}")

    # sort poles (and residues) by the real component of the pole
    for ip in range(n_pieces):
        indices = mp_poles[ip].argsort()
        mp_poles[ip] = mp_poles[ip][indices]
        mp_residues[ip] = mp_residues[ip][:, indices]

    # initialize an array to record whether each pole is used or not
    poles_unused = [np.ones_like(p, dtype=int) for p in mp_poles]

    # optimize the windows: the goal is to find the least set of significant
    # consecutive poles and curve fit coefficients to reproduce cross section
    win_data = []
    for iw in range(n_win):
        if log >= DETAILED_LOGGING:
            print(f"Processing window {iw + 1}/{n_win}...")

        # inner window boundaries
        inbegin = sqrt(E_min) + spacing * iw
        inend = inbegin + spacing
        incenter = (inbegin + inend) / 2.0
        # extend window energy range for Doppler broadening
        if iw == 0 or sqrt(alpha)*inbegin < 4.0:
            e_start = inbegin**2
        else:
            e_start = max(E_min, (sqrt(alpha)*inbegin - 4.0)**2/alpha)
        e_end = min(E_max, (sqrt(alpha)*inend + 4.0)**2/alpha)

        # locate piece and relevant poles
        i_piece = min(n_pieces - 1, int((inbegin - sqrt(E_min))/piece_width + 0.5))
        poles, residues = mp_poles[i_piece], mp_residues[i_piece]
        n_poles = poles.size

        # generate energy points for fitting: equally spaced in momentum
        n_points = min(max(100, int((e_end - e_start)*4)), 10000)
        energy_sqrt = np.linspace(np.sqrt(e_start), np.sqrt(e_end), n_points)
        energy = energy_sqrt**2

        # reference xs from multipole form, note the residue terms in the
        # multipole and vector fitting representations differ by a 1j
        xs_ref = vf.evaluate(energy_sqrt, poles, residues*1j) / energy

        # curve fit matrix
        matrix = np.vstack([energy**(0.5*i - 1) for i in range(n_cf + 1)]).T

        # start from 0 poles, initialize pointers to the center nearest pole
        center_pole_ind = np.argmin((np.fabs(poles.real - incenter)))
        lp = rp = center_pole_ind
        while True:
            if log >= DETAILED_LOGGING:
                print(f"Trying poles {lp} to {rp}")

            # calculate the cross sections contributed by the windowed poles
            if rp > lp:
                xs_wp = vf.evaluate(energy_sqrt, poles[lp:rp],
                                    residues[:, lp:rp]*1j) / energy
            else:
                xs_wp = np.zeros_like(xs_ref)

            # do least square curve fit on the remains
            coefs = np.linalg.lstsq(matrix, (xs_ref - xs_wp).T, rcond=None)[0]
            xs_fit = (matrix @ coefs).T

            # assess the result
            abserr = np.abs(xs_fit + xs_wp - xs_ref)
            with np.errstate(invalid='ignore', divide='ignore'):
                relerr = abserr / xs_ref
            if not np.any(np.isnan(abserr)):
                re = relerr[abserr > atol]
                if re.size == 0 or np.all(re <= rtol) or \
                   (re.max() <= 2*rtol and (re > rtol).sum() <= 0.01*relerr.size) or \
                   (iw == 0 and np.all(relerr.mean(axis=1) <= rtol)):
                    # meet tolerances
                    if log >= DETAILED_LOGGING:
                        print("Accuracy satisfied.")
                    break

            # we expect pure curvefit will succeed for the first window
            # TODO: find the energy boundary below which no poles are allowed
            if iw == 0:
                raise RuntimeError('Pure curvefit failed for the first window!')

            # try to include one more pole (next center nearest)
            if rp >= n_poles:
                lp -= 1
            elif lp <= 0 or poles[rp] - incenter <= incenter - poles[lp - 1]:
                rp += 1
            else:
                lp -= 1

        # save data for this window
        win_data.append((i_piece, lp, rp, coefs))

        # mark the windowed poles as used poles
        poles_unused[i_piece][lp:rp] = 0

    # flatten and shrink by removing unused poles
    data = []  # used poles and residues
    for ip in range(n_pieces):
        used = (poles_unused[ip] == 0)
        # stack poles and residues for library format
        data.append(np.vstack([mp_poles[ip][used], mp_residues[ip][:, used]]).T)
    # stack poles/residues in sequence vertically
    data = np.vstack(data)
    # new start/end pole indices
    windows = []
    curvefit = []
    for iw in range(n_win):
        ip, lp, rp, coefs = win_data[iw]
        # adjust indices and change to 1-based for the library format
        n_prev_poles = sum([poles_unused[i].size for i in range(ip)])
        n_unused = sum([(poles_unused[i] == 1).sum() for i in range(ip)]) + \
                  (poles_unused[ip][:lp] == 1).sum()
        lp += n_prev_poles - n_unused + 1
        rp += n_prev_poles - n_unused
        windows.append([lp, rp])
        curvefit.append(coefs)

    # construct the WindowedMultipole object
    wmp = WindowedMultipole(name)
    wmp.spacing = spacing
    wmp.sqrtAWR = sqrt(awr)
    wmp.E_min = E_min
    wmp.E_max = E_max
    wmp.data = data
    wmp.windows = np.asarray(windows)
    wmp.curvefit = np.asarray(curvefit)
    # TODO: check if Doppler brodening of the polynomial curvefit is negligible
    wmp.broaden_poly = np.ones((n_win,), dtype=bool)

    return wmp


class WindowedMultipole(EqualityMixin):
    """Resonant cross sections represented in the windowed multipole format.

    Parameters
    ----------
    name : str
        Name of the nuclide using the GNDS naming convention

    Attributes
    ----------
    name : str
        Name of the nuclide using the GNDS naming convention
    spacing : float
        The width of each window in sqrt(E)-space.  For example, the frst window
        will end at (sqrt(E_min) + spacing)**2 and the second window at
        (sqrt(E_min) + 2*spacing)**2.
    sqrtAWR : float
        Square root of the atomic weight ratio of the target nuclide.
    E_min : float
        Lowest energy in eV the library is valid for.
    E_max : float
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

    @name.setter
    def name(self, name):
        cv.check_type('name', name, str)
        self._name = name

    @property
    def fit_order(self):
        return self.curvefit.shape[1] - 1

    @property
    def fissionable(self):
        return self.data.shape[1] == 4

    @property
    def n_poles(self):
        return self.data.shape[0]

    @property
    def n_windows(self):
        return self.windows.shape[0]

    @property
    def poles_per_window(self):
        return (self.windows[:, 1] - self.windows[:, 0] + 1).mean()

    @property
    def spacing(self):
        return self._spacing

    @spacing.setter
    def spacing(self, spacing):
        if spacing is not None:
            cv.check_type('spacing', spacing, Real)
            cv.check_greater_than('spacing', spacing, 0.0, equality=False)
        self._spacing = spacing

    @property
    def sqrtAWR(self):
        return self._sqrtAWR

    @sqrtAWR.setter
    def sqrtAWR(self, sqrtAWR):
        if sqrtAWR is not None:
            cv.check_type('sqrtAWR', sqrtAWR, Real)
            cv.check_greater_than('sqrtAWR', sqrtAWR, 0.0, equality=False)
        self._sqrtAWR = sqrtAWR

    @property
    def E_min(self):
        return self._E_min

    @E_min.setter
    def E_min(self, E_min):
        if E_min is not None:
            cv.check_type('E_min', E_min, Real)
            cv.check_greater_than('E_min', E_min, 0.0, equality=True)
        self._E_min = E_min

    @property
    def E_max(self):
        return self._E_max

    @E_max.setter
    def E_max(self, E_max):
        if E_max is not None:
            cv.check_type('E_max', E_max, Real)
            cv.check_greater_than('E_max', E_max, 0.0, equality=False)
        self._E_max = E_max

    @property
    def data(self):
        return self._data

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

    @property
    def windows(self):
        return self._windows

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

    @property
    def broaden_poly(self):
        return self._broaden_poly

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

    @property
    def curvefit(self):
        return self._curvefit

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

        out.broaden_poly = group['broaden_poly'][...].astype(bool)
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
    def from_endf(cls, endf_file, log=False, vf_options=None, wmp_options=None):
        """Generate windowed multipole neutron data from an ENDF evaluation.

        .. versionadded:: 0.12.1

        Parameters
        ----------
        endf_file : str
            Path to ENDF evaluation
        log : bool or int, optional
            Whether to print running logs (use int for verbosity control)
        vf_options : dict, optional
            Dictionary of keyword arguments, e.g. {'njoy_error': 0.001},
            passed to :func:`openmc.data.multipole.vectfit_nuclide`
        wmp_options : dict, optional
            Dictionary of keyword arguments, e.g. {'search': True, 'rtol': 0.01},
            passed to :func:`openmc.data.WindowedMultipole.from_multipole`

        Returns
        -------
        openmc.data.WindowedMultipole
            Resonant cross sections represented in the windowed multipole
            format.

        """

        if vf_options is None:
            vf_options = {}

        if wmp_options is None:
            wmp_options = {}

        if log:
            vf_options.update(log=log)
            wmp_options.update(log=log)

        # generate multipole data from EDNF
        mp_data = vectfit_nuclide(endf_file, **vf_options)

        # windowing
        return cls.from_multipole(mp_data, **wmp_options)

    @classmethod
    def from_multipole(cls, mp_data, search=None, log=False, **kwargs):
        """Generate windowed multipole neutron data from multipole data.

        Parameters
        ----------
        mp_data : dictionary or str
            Dictionary or Path to the multipole data stored in a pickle file
        search : bool, optional
            Whether to search for optimal window size and curvefit order.
            Defaults to True if no windowing parameters are specified.
        log : bool or int, optional
            Whether to print running logs (use int for verbosity control)
        **kwargs
            Keyword arguments passed to :func:`openmc.data.multipole._windowing`

        Returns
        -------
        openmc.data.WindowedMultipole
            Resonant cross sections represented in the windowed multipole
            format.

        """

        if isinstance(mp_data, str):
            # load multipole data from file
            with open(mp_data, 'rb') as f:
                mp_data = pickle.load(f)

        if search is None:
            if 'n_cf' in kwargs and ('n_win' in kwargs or 'spacing' in kwargs):
                search = False
            else:
                search = True

        # windowing with specific options
        if not search:
            # set default value for curvefit order if not specified
            if 'n_cf' not in kwargs:
                kwargs.update(n_cf=5)
            return _windowing(mp_data, log=log, **kwargs)

        # search optimal WMP from a range of window sizes and CF orders
        if log:
            print("Start searching ...")
        n_poles = sum([p.size for p in mp_data["poles"]])
        n_win_min = max(5, n_poles // 20)
        n_win_max = 2000 if n_poles < 2000 else 8000
        best_wmp = best_metric = None
        for n_w in np.unique(np.linspace(n_win_min, n_win_max, 20, dtype=int)):
            for n_cf in range(10, 1, -1):
                if log:
                    print(f"Testing N_win={n_w} N_cf={n_cf}")

                # update arguments dictionary
                kwargs.update(n_win=n_w, n_cf=n_cf)

                # windowing
                try:
                    wmp = _windowing(mp_data, log=log, **kwargs)
                except Exception as e:
                    if log:
                        print('Failed: ' + str(e))
                    break

                # select wmp library with metric:
                # - performance: average # used poles per window and CF order
                # - memory: # windows
                metric = -(wmp.poles_per_window * 10. + wmp.fit_order * 1. +
                           wmp.n_windows * 0.01)
                if best_wmp is None or metric > best_metric:
                    if log:
                        print("Best library so far.")
                    best_wmp = deepcopy(wmp)
                    best_metric = metric

        # return the best wmp library
        if log:
            print("Final library: {} poles, {} windows, {:.2g} poles per window, "
                  "{} CF order".format(best_wmp.n_poles, best_wmp.n_windows,
                   best_wmp.poles_per_window, best_wmp.fit_order))

        return best_wmp

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
            Scattering, absorption, and fission microscopic cross sections
            at the given energy and temperature.

        """

        if E < self.E_min: return (0, 0, 0)
        if E > self.E_max: return (0, 0, 0)

        # ======================================================================
        # Bookkeeping

        # Define some frequently used variables.
        sqrtkT = sqrt(K_BOLTZMANN * T)
        sqrtE = sqrt(E)
        invE = 1.0 / E

        # Locate us.  The i_window calc omits a + 1 present from the legacy
        # Fortran version of OpenMC because of the 1-based vs. 0-based
        # indexing.  Similarly startw needs to be decreased by 1.  endw does
        # not need to be decreased because range(startw, endw) does not include
        # endw.
        i_window = min(self.n_windows - 1,
                       int(np.floor((sqrtE - sqrt(self.E_min)) / self.spacing)))
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
            for i_poly in range(self.fit_order + 1):
                sig_s += (self.curvefit[i_window, i_poly, _FIT_S]
                          * broadened_polynomials[i_poly])
                sig_a += (self.curvefit[i_window, i_poly, _FIT_A]
                          * broadened_polynomials[i_poly])
                if self.fissionable:
                    sig_f += (self.curvefit[i_window, i_poly, _FIT_F]
                              * broadened_polynomials[i_poly])
        else:
            temp = invE
            for i_poly in range(self.fit_order + 1):
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
            Scattering, absorption, and fission microscopic cross sections
            at the given energy and temperature.

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
            f.attrs['filetype'] = np.bytes_('data_wmp')
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
