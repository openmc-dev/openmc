"""Chebyshev Rational Approximation Method module

Implements two different forms of CRAM for use in openmc.deplete.
"""

from itertools import repeat
from multiprocessing import Pool
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

from .. import comm


def deplete(chain, x, op_result, dt, print_out):
    """Deplete materials using given reaction rates for a specified time

    Parameters
    ----------
    chain : openmc.deplete.Chain
        Depletion chain
    x : list of numpy.ndarray
        Atom number vectors for each material
    op_result : openmc.deplete.OperatorResult
        Result of applying transport operator (contains reaction rates)
    dt : float
        Time in [s] to deplete for
    print_out : bool
        Whether to show elapsed time

    Returns
    -------
    x_result : list of numpy.ndarray
        Updated atom number vectors for each material

    """
    t_start = time.time()

    # Set up iterators
    n_mats = len(x)
    chains = repeat(chain, n_mats)
    vecs = (x[i] for i in range(n_mats))
    rates = (op_result.rates[i, :, :] for i in range(n_mats))
    dts = repeat(dt, n_mats)

    # Use multiprocessing pool to distribute work
    with Pool() as pool:
        iters = zip(chains, vecs, rates, dts)
        x_result = list(pool.starmap(_cram_wrapper, iters))

    t_end = time.time()
    if comm.rank == 0:
        if print_out:
            print("Time to matexp: ", t_end - t_start)

    return x_result


def _cram_wrapper(chain, n0, rates, dt):
    """Wraps depletion matrix creation / CRAM solve for multiprocess execution

    Parameters
    ----------
    chain : DepletionChain
        Depletion chain used to construct the burnup matrix
    n0 : numpy.array
        Vector to operate a matrix exponent on.
    rates : numpy.ndarray
        2D array indexed by nuclide then by cell.
    dt : float
        Time to integrate to.

    Returns
    -------
    numpy.array
        Results of the matrix exponent.
    """
    A = chain.form_matrix(rates)
    return CRAM48(A, n0, dt)


def CRAM16(A, n0, dt):
    """Chebyshev Rational Approximation Method, order 16

    Algorithm is the 16th order Chebyshev Rational Approximation Method,
    implemented in the more stable `incomplete partial fraction (IPF)
    <https://doi.org/10.13182/NSE15-26>`_ form.

    Parameters
    ----------
    A : scipy.linalg.csr_matrix
        Matrix to take exponent of.
    n0 : numpy.array
        Vector to operate a matrix exponent on.
    dt : float
        Time to integrate to.

    Returns
    -------
    numpy.array
        Results of the matrix exponent.

    """

    alpha = np.array([+2.124853710495224e-16,
                      +5.464930576870210e+3 - 3.797983575308356e+4j,
                      +9.045112476907548e+1 - 1.115537522430261e+3j,
                      +2.344818070467641e+2 - 4.228020157070496e+2j,
                      +9.453304067358312e+1 - 2.951294291446048e+2j,
                      +7.283792954673409e+2 - 1.205646080220011e+5j,
                      +3.648229059594851e+1 - 1.155509621409682e+2j,
                      +2.547321630156819e+1 - 2.639500283021502e+1j,
                      +2.394538338734709e+1 - 5.650522971778156e+0j],
                     dtype=np.complex128)
    theta = np.array([+0.0,
                      +3.509103608414918 + 8.436198985884374j,
                      +5.948152268951177 + 3.587457362018322j,
                      -5.264971343442647 + 16.22022147316793j,
                      +1.419375897185666 + 10.92536348449672j,
                      +6.416177699099435 + 1.194122393370139j,
                      +4.993174737717997 + 5.996881713603942j,
                      -1.413928462488886 + 13.49772569889275j,
                      -10.84391707869699 + 19.27744616718165j],
                     dtype=np.complex128)

    n = A.shape[0]

    alpha0 = 2.124853710495224e-16

    k = 8

    y = np.array(n0, dtype=np.float64)
    for l in range(1, k+1):
        y = 2.0*np.real(alpha[l]*sla.spsolve(A*dt - theta[l]*sp.eye(n), y)) + y

    y *= alpha0
    return y


def CRAM48(A, n0, dt):
    """Chebyshev Rational Approximation Method, order 48

    Algorithm is the 48th order Chebyshev Rational Approximation Method,
    implemented in the more stable `incomplete partial fraction (IPF)
    <https://doi.org/10.13182/NSE15-26>`_ form.

    Parameters
    ----------
    A : scipy.linalg.csr_matrix
        Matrix to take exponent of.
    n0 : numpy.array
        Vector to operate a matrix exponent on.
    dt : float
        Time to integrate to.

    Returns
    -------
    numpy.array
        Results of the matrix exponent.

    """

    theta_r = np.array([-4.465731934165702e+1, -5.284616241568964e+0,
                        -8.867715667624458e+0, +3.493013124279215e+0,
                        +1.564102508858634e+1, +1.742097597385893e+1,
                        -2.834466755180654e+1, +1.661569367939544e+1,
                        +8.011836167974721e+0, -2.056267541998229e+0,
                        +1.449208170441839e+1, +1.853807176907916e+1,
                        +9.932562704505182e+0, -2.244223871767187e+1,
                        +8.590014121680897e-1, -1.286192925744479e+1,
                        +1.164596909542055e+1, +1.806076684783089e+1,
                        +5.870672154659249e+0, -3.542938819659747e+1,
                        +1.901323489060250e+1, +1.885508331552577e+1,
                        -1.734689708174982e+1, +1.316284237125190e+1])
    theta_i = np.array([+6.233225190695437e+1, +4.057499381311059e+1,
                        +4.325515754166724e+1, +3.281615453173585e+1,
                        +1.558061616372237e+1, +1.076629305714420e+1,
                        +5.492841024648724e+1, +1.316994930024688e+1,
                        +2.780232111309410e+1, +3.794824788914354e+1,
                        +1.799988210051809e+1, +5.974332563100539e+0,
                        +2.532823409972962e+1, +5.179633600312162e+1,
                        +3.536456194294350e+1, +4.600304902833652e+1,
                        +2.287153304140217e+1, +8.368200580099821e+0,
                        +3.029700159040121e+1, +5.834381701800013e+1,
                        +1.194282058271408e+0, +3.583428564427879e+0,
                        +4.883941101108207e+1, +2.042951874827759e+1])
    theta = np.array(theta_r + theta_i * 1j, dtype=np.complex128)

    alpha_r = np.array([+6.387380733878774e+2, +1.909896179065730e+2,
                        +4.236195226571914e+2, +4.645770595258726e+2,
                        +7.765163276752433e+2, +1.907115136768522e+3,
                        +2.909892685603256e+3, +1.944772206620450e+2,
                        +1.382799786972332e+5, +5.628442079602433e+3,
                        +2.151681283794220e+2, +1.324720240514420e+3,
                        +1.617548476343347e+4, +1.112729040439685e+2,
                        +1.074624783191125e+2, +8.835727765158191e+1,
                        +9.354078136054179e+1, +9.418142823531573e+1,
                        +1.040012390717851e+2, +6.861882624343235e+1,
                        +8.766654491283722e+1, +1.056007619389650e+2,
                        +7.738987569039419e+1, +1.041366366475571e+2])
    alpha_i = np.array([-6.743912502859256e+2, -3.973203432721332e+2,
                        -2.041233768918671e+3, -1.652917287299683e+3,
                        -1.783617639907328e+4, -5.887068595142284e+4,
                        -9.953255345514560e+3, -1.427131226068449e+3,
                        -3.256885197214938e+6, -2.924284515884309e+4,
                        -1.121774011188224e+3, -6.370088443140973e+4,
                        -1.008798413156542e+6, -8.837109731680418e+1,
                        -1.457246116408180e+2, -6.388286188419360e+1,
                        -2.195424319460237e+2, -6.719055740098035e+2,
                        -1.693747595553868e+2, -1.177598523430493e+1,
                        -4.596464999363902e+3, -1.738294585524067e+3,
                        -4.311715386228984e+1, -2.777743732451969e+2])
    alpha = np.array(alpha_r + alpha_i * 1j, dtype=np.complex128)
    n = A.shape[0]

    alpha0 = 2.258038182743983e-47

    k = 24

    y = np.array(n0, dtype=np.float64)
    for l in range(k):
        y = 2.0*np.real(alpha[l]*sla.spsolve(A*dt - theta[l]*sp.eye(n), y)) + y

    y *= alpha0
    return y
