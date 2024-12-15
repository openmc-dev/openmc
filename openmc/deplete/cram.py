"""Chebyshev Rational Approximation Method module

Implements two different forms of CRAM for use in openmc.deplete.
"""

import numbers

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

from openmc.checkvalue import check_type, check_length
from .abc import DepSystemSolver

__all__ = ["CRAM16", "CRAM48", "Cram16Solver", "Cram48Solver", "IPFCramSolver"]


class IPFCramSolver(DepSystemSolver):
    r"""CRAM depletion solver that uses incomplete partial factorization

    Provides a :meth:`__call__` that utilizes an incomplete
    partial factorization (IPF) for the Chebyshev Rational Approximation
    Method (CRAM), as described in the following paper: M. Pusa, "`Higher-Order
    Chebyshev Rational Approximation Method and Application to Burnup Equations
    <https://doi.org/10.13182/NSE15-26>`_," Nucl. Sci. Eng., 182:3, 297-318.

    Parameters
    ----------
    alpha : numpy.ndarray
        Complex residues of poles used in the factorization. Must be a
        vector with even number of items.
    theta : numpy.ndarray
        Complex poles. Must have an equal size as ``alpha``.
    alpha0 : float
        Limit of the approximation at infinity

    Attributes
    ----------
    alpha : numpy.ndarray
        Complex residues of poles :attr:`theta` in the incomplete partial
        factorization. Denoted as :math:`\tilde{\alpha}`
    theta : numpy.ndarray
        Complex poles :math:`\theta` of the rational approximation
    alpha0 : float
        Limit of the approximation at infinity

    """

    def __init__(self, alpha, theta, alpha0):
        check_type("alpha", alpha, np.ndarray, numbers.Complex)
        check_type("theta", theta, np.ndarray, numbers.Complex)
        check_length("theta", theta, alpha.size)
        check_type("alpha0", alpha0, numbers.Real)
        self.alpha = alpha
        self.theta = theta
        self.alpha0 = alpha0

    def __call__(self, A, n0, dt):
        """Solve depletion equations using IPF CRAM

        Parameters
        ----------
        A : scipy.sparse.csr_matrix
            Sparse transmutation matrix ``A[j, i]`` desribing rates at
            which isotope ``i`` transmutes to isotope ``j``
        n0 : numpy.ndarray
            Initial compositions, typically given in number of atoms in some
            material or an atom density
        dt : float
            Time [s] of the specific interval to be solved

        Returns
        -------
        numpy.ndarray
            Final compositions after ``dt``

        """
        A = dt * sp.csc_matrix(A, dtype=np.float64)
        y = n0.copy()
        ident = sp.eye(A.shape[0], format="csc")
        for alpha, theta in zip(self.alpha, self.theta):
            y += 2 * np.real(alpha * sla.spsolve(A - theta * ident, y))
        return y * self.alpha0


# Coefficients for IPF Cram 16
c16_alpha = np.array(
    [
        +5.464930576870210e3 - 3.797983575308356e4j,
        +9.045112476907548e1 - 1.115537522430261e3j,
        +2.344818070467641e2 - 4.228020157070496e2j,
        +9.453304067358312e1 - 2.951294291446048e2j,
        +7.283792954673409e2 - 1.205646080220011e5j,
        +3.648229059594851e1 - 1.155509621409682e2j,
        +2.547321630156819e1 - 2.639500283021502e1j,
        +2.394538338734709e1 - 5.650522971778156e0j,
    ],
    dtype=np.complex128,
)

c16_theta = np.array(
    [
        +3.509103608414918 + 8.436198985884374j,
        +5.948152268951177 + 3.587457362018322j,
        -5.264971343442647 + 16.22022147316793j,
        +1.419375897185666 + 10.92536348449672j,
        +6.416177699099435 + 1.194122393370139j,
        +4.993174737717997 + 5.996881713603942j,
        -1.413928462488886 + 13.49772569889275j,
        -10.84391707869699 + 19.27744616718165j,
    ],
    dtype=np.complex128,
)

c16_alpha0 = 2.124853710495224e-16
Cram16Solver = IPFCramSolver(c16_alpha, c16_theta, c16_alpha0)
CRAM16 = Cram16Solver.__call__

del c16_alpha, c16_alpha0, c16_theta

# Coefficients for 48th order IPF Cram

theta_r = np.array(
    [
        -4.465731934165702e1,
        -5.284616241568964e0,
        -8.867715667624458e0,
        +3.493013124279215e0,
        +1.564102508858634e1,
        +1.742097597385893e1,
        -2.834466755180654e1,
        +1.661569367939544e1,
        +8.011836167974721e0,
        -2.056267541998229e0,
        +1.449208170441839e1,
        +1.853807176907916e1,
        +9.932562704505182e0,
        -2.244223871767187e1,
        +8.590014121680897e-1,
        -1.286192925744479e1,
        +1.164596909542055e1,
        +1.806076684783089e1,
        +5.870672154659249e0,
        -3.542938819659747e1,
        +1.901323489060250e1,
        +1.885508331552577e1,
        -1.734689708174982e1,
        +1.316284237125190e1,
    ]
)

theta_i = np.array(
    [
        +6.233225190695437e1,
        +4.057499381311059e1,
        +4.325515754166724e1,
        +3.281615453173585e1,
        +1.558061616372237e1,
        +1.076629305714420e1,
        +5.492841024648724e1,
        +1.316994930024688e1,
        +2.780232111309410e1,
        +3.794824788914354e1,
        +1.799988210051809e1,
        +5.974332563100539e0,
        +2.532823409972962e1,
        +5.179633600312162e1,
        +3.536456194294350e1,
        +4.600304902833652e1,
        +2.287153304140217e1,
        +8.368200580099821e0,
        +3.029700159040121e1,
        +5.834381701800013e1,
        +1.194282058271408e0,
        +3.583428564427879e0,
        +4.883941101108207e1,
        +2.042951874827759e1,
    ]
)

c48_theta = np.array(theta_r + theta_i * 1j, dtype=np.complex128)

alpha_r = np.array(
    [
        +6.387380733878774e2,
        +1.909896179065730e2,
        +4.236195226571914e2,
        +4.645770595258726e2,
        +7.765163276752433e2,
        +1.907115136768522e3,
        +2.909892685603256e3,
        +1.944772206620450e2,
        +1.382799786972332e5,
        +5.628442079602433e3,
        +2.151681283794220e2,
        +1.324720240514420e3,
        +1.617548476343347e4,
        +1.112729040439685e2,
        +1.074624783191125e2,
        +8.835727765158191e1,
        +9.354078136054179e1,
        +9.418142823531573e1,
        +1.040012390717851e2,
        +6.861882624343235e1,
        +8.766654491283722e1,
        +1.056007619389650e2,
        +7.738987569039419e1,
        +1.041366366475571e2,
    ]
)

alpha_i = np.array(
    [
        -6.743912502859256e2,
        -3.973203432721332e2,
        -2.041233768918671e3,
        -1.652917287299683e3,
        -1.783617639907328e4,
        -5.887068595142284e4,
        -9.953255345514560e3,
        -1.427131226068449e3,
        -3.256885197214938e6,
        -2.924284515884309e4,
        -1.121774011188224e3,
        -6.370088443140973e4,
        -1.008798413156542e6,
        -8.837109731680418e1,
        -1.457246116408180e2,
        -6.388286188419360e1,
        -2.195424319460237e2,
        -6.719055740098035e2,
        -1.693747595553868e2,
        -1.177598523430493e1,
        -4.596464999363902e3,
        -1.738294585524067e3,
        -4.311715386228984e1,
        -2.777743732451969e2,
    ]
)

c48_alpha = np.array(alpha_r + alpha_i * 1j, dtype=np.complex128)

c48_alpha0 = 2.258038182743983e-47

Cram48Solver = IPFCramSolver(c48_alpha, c48_theta, c48_alpha0)

del c48_alpha, c48_alpha0, c48_theta, alpha_r, alpha_i, theta_r, theta_i

CRAM48 = Cram48Solver.__call__
