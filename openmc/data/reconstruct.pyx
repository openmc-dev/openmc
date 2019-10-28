from libc.stdlib cimport malloc, calloc, free
from libc.math cimport cos, sin, sqrt, atan, M_PI

cimport numpy as np
import numpy as np
from numpy.linalg import inv
cimport cython


cdef extern from "complex.h":
    double cabs(double complex)
    double complex conj(double complex)
    double creal(complex double)
    double cimag(complex double)
    double complex cexp(double complex)

# Physical constants are from CODATA 2014
cdef double NEUTRON_MASS_ENERGY = 939.5654133e6  # eV/c^2
cdef double HBAR_C = 197.3269788e5  # eV-b^0.5


@cython.cdivision(True)
def wave_number(double A, double E):
    r"""Neutron wave number in center-of-mass system.

    ENDF-102 defines the neutron wave number in the center-of-mass system in
    Equation D.10 as

    .. math::
        k = \frac{2m_n}{\hbar} \frac{A}{A + 1} \sqrt{|E|}

    Parameters
    ----------
    A : double
        Ratio of target mass to neutron mass
    E : double
        Energy in eV

    Returns
    -------
    double
        Neutron wave number in b^-0.5

    """
    return A/(A + 1)*sqrt(2*NEUTRON_MASS_ENERGY*abs(E))/HBAR_C

@cython.cdivision(True)
cdef double _wave_number(double A, double E):
    return A/(A + 1)*sqrt(2*NEUTRON_MASS_ENERGY*abs(E))/HBAR_C


@cython.cdivision(True)
cdef double phaseshift(int l, double rho):
    """Calculate hardsphere phase shift as given in ENDF-102, Equation D.13

    Parameters
    ----------
    l : int
        Angular momentum quantum number
    rho : float
        Product of the wave number and the channel radius

    Returns
    -------
    double
        Hardsphere phase shift

    """
    if l == 0:
        return rho
    elif l == 1:
        return rho - atan(rho)
    elif l == 2:
        return rho - atan(3*rho/(3 - rho**2))
    elif l == 3:
        return rho - atan((15*rho - rho**3)/(15 - 6*rho**2))
    elif l == 4:
        return rho - atan((105*rho - 10*rho**3)/(105 - 45*rho**2 + rho**4))


@cython.cdivision(True)
def penetration_shift(int l, double rho):
    r"""Calculate shift and penetration factors as given in ENDF-102, Equations D.11
    and D.12.

    Parameters
    ----------
    l : int
        Angular momentum quantum number
    rho : float
        Product of the wave number and the channel radius

    Returns
    -------
    double
        Penetration factor for given :math:`l`
    double
        Shift factor for given :math:`l`

    """
    cdef double den

    if l == 0:
        return rho, 0.
    elif l == 1:
        den = 1 + rho**2
        return rho**3/den, -1/den
    elif l == 2:
        den = 9 + 3*rho**2 + rho**4
        return rho**5/den, -(18 + 3*rho**2)/den
    elif l == 3:
        den = 225 + 45*rho**2 + 6*rho**4 + rho**6
        return rho**7/den, -(675 + 90*rho**2 + 6*rho**4)/den
    elif l == 4:
        den = 11025 + 1575*rho**2 + 135*rho**4 + 10*rho**6 + rho**8
        return rho**9/den, -(44100 + 4725*rho**2 + 270*rho**4 + 10*rho**6)/den


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def reconstruct_mlbw(mlbw, double E):
    """Evaluate cross section using MLBW data.

    Parameters
    ----------
    mlbw : openmc.data.MultiLevelBreitWigner
        Multi-level Breit-Wigner resonance parameters
    E : double
        Energy in eV at which to evaluate the cross section

    Returns
    -------
    elastic : double
        Elastic scattering cross section in barns
    capture : double
        Radiative capture cross section in barns
    fission : double
        Fission cross section in barns

    """
    cdef int i, nJ, ij, l, n_res, i_res
    cdef double elastic, capture, fission
    cdef double A, k, rho, rhohat, I
    cdef double P, S, phi, cos2phi, sin2phi
    cdef double Ex, Q, rhoc, rhochat, P_c, S_c
    cdef double jmin, jmax, j, Dl
    cdef double E_r, gt, gn, gg, gf, gx, P_r, S_r, P_rx
    cdef double gnE, gtE, Eprime, x, f
    cdef double *g
    cdef double (*s)[2]
    cdef double [:,:] params

    I = mlbw.target_spin
    A = mlbw.atomic_weight_ratio
    k = _wave_number(A, E)

    elastic = 0.
    capture = 0.
    fission = 0.

    for i, l in enumerate(mlbw._l_values):
        params = mlbw._parameter_matrix[l]

        rho = k*mlbw.channel_radius[l](E)
        rhohat = k*mlbw.scattering_radius[l](E)
        P, S = penetration_shift(l, rho)
        phi = phaseshift(l, rhohat)
        cos2phi = cos(2*phi)
        sin2phi = sin(2*phi)

        # Determine shift and penetration at modified energy
        if mlbw._competitive[i]:
            Ex = E + mlbw.q_value[l]*(A + 1)/A
            rhoc = mlbw.channel_radius[l](Ex)
            rhochat = mlbw.scattering_radius[l](Ex)
            P_c, S_c = penetration_shift(l, rhoc)
            if Ex < 0:
                P_c = 0

        # Determine range of total angular momentum values based on equation
        # 41 in LA-UR-12-27079
        jmin = abs(abs(I - l) - 0.5)
        jmax = I + l + 0.5
        nJ = int(jmax - jmin + 1)

        # Determine Dl factor using Equation 43 in LA-UR-12-27079
        Dl = 2*l + 1
        g = <double *> malloc(nJ*sizeof(double))
        for ij in range(nJ):
            j = jmin + ij
            g[ij] = (2*j + 1)/(4*I + 2)
            Dl -= g[ij]

        s = <double (*)[2]> calloc(2*nJ, sizeof(double))
        for i_res in range(params.shape[0]):
            # Copy resonance parameters
            E_r = params[i_res, 0]
            j = params[i_res, 2]
            ij = int(j - jmin)
            gt = params[i_res, 3]
            gn = params[i_res, 4]
            gg = params[i_res, 5]
            gf = params[i_res, 6]
            gx = params[i_res, 7]
            P_r = params[i_res, 8]
            S_r = params[i_res, 9]
            P_rx = params[i_res, 10]

            # Calculate neutron and total width at energy E
            gnE = P*gn/P_r  # ENDF-102, Equation D.7
            gtE = gnE + gg + gf
            if gx > 0:
                gtE += gx*P_c/P_rx

            Eprime = E_r + (S_r - S)/(2*P_r)*gn  # ENDF-102, Equation D.9
            x = 2*(E - Eprime)/gtE    # LA-UR-12-27079, Equation 26
            f = 2*gnE/(gtE*(1 + x*x)) # Common factor in Equation 40
            s[ij][0] += f             # First sum in Equation 40
            s[ij][1] += f*x           # Second sum in Equation 40
            capture += f*g[ij]*gg/gtE
            if gf > 0:
                fission += f*g[ij]*gf/gtE

        for ij in range(nJ):
            # Add all but last term of LA-UR-12-27079, Equation 40
            elastic += g[ij]*((1 - cos2phi - s[ij][0])**2 +
                              (sin2phi + s[ij][1])**2)

        # Add final term with Dl from Equation 40
        elastic += 2*Dl*(1 - cos2phi)

        # Free memory
        free(g)
        free(s)

    capture *= 2*M_PI/(k*k)
    fission *= 2*M_PI/(k*k)
    elastic *= M_PI/(k*k)

    return (elastic, capture, fission)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def reconstruct_slbw(slbw, double E):
    """Evaluate cross section using SLBW data.

    Parameters
    ----------
    slbw : openmc.data.SingleLevelBreitWigner
        Single-level Breit-Wigner resonance parameters
    E : double
        Energy in eV at which to evaluate the cross section

    Returns
    -------
    elastic : double
        Elastic scattering cross section in barns
    capture : double
        Radiative capture cross section in barns
    fission : double
        Fission cross section in barns

    """
    cdef int i, l, i_res
    cdef double elastic, capture, fission
    cdef double A, k, rho, rhohat, I
    cdef double P, S, phi, cos2phi, sin2phi, sinphi2
    cdef double Ex, rhoc, rhochat, P_c, S_c
    cdef double E_r, J, gt, gn, gg, gf, gx, P_r, S_r, P_rx
    cdef double gnE, gtE, Eprime, f
    cdef double x, theta, psi, chi
    cdef double [:,:] params

    I = slbw.target_spin
    A = slbw.atomic_weight_ratio
    k = _wave_number(A, E)

    elastic = 0.
    capture = 0.
    fission = 0.

    for i, l in enumerate(slbw._l_values):
        params = slbw._parameter_matrix[l]

        rho = k*slbw.channel_radius[l](E)
        rhohat = k*slbw.scattering_radius[l](E)
        P, S = penetration_shift(l, rho)
        phi = phaseshift(l, rhohat)
        cos2phi = cos(2*phi)
        sin2phi = sin(2*phi)
        sinphi2 = sin(phi)**2

        # Add potential scattering -- first term in ENDF-102, Equation D.2
        elastic += 4*M_PI/(k*k)*(2*l + 1)*sinphi2

        # Determine shift and penetration at modified energy
        if slbw._competitive[i]:
            Ex = E + slbw.q_value[l]*(A + 1)/A
            rhoc = slbw.channel_radius[l](Ex)
            rhochat = slbw.scattering_radius[l](Ex)
            P_c, S_c = penetration_shift(l, rhoc)
            if Ex < 0:
                P_c = 0

        for i_res in range(params.shape[0]):
            # Copy resonance parameters
            E_r = params[i_res, 0]
            J = params[i_res, 2]
            gt = params[i_res, 3]
            gn = params[i_res, 4]
            gg = params[i_res, 5]
            gf = params[i_res, 6]
            gx = params[i_res, 7]
            P_r = params[i_res, 8]
            S_r = params[i_res, 9]
            P_rx = params[i_res, 10]

            # Calculate neutron and total width at energy E
            gnE = P*gn/P_r  # Equation D.7
            gtE = gnE + gg + gf
            if gx > 0:
                gtE += gx*P_c/P_rx

            Eprime = E_r + (S_r - S)/(2*P_r)*gn  # Equation D.9
            gJ = (2*J + 1)/(4*I + 2)  # Mentioned in section D.1.1.4

            # Calculate common factor for elastic, capture, and fission
            # cross sections
            f = M_PI/(k*k)*gJ*gnE/((E - Eprime)**2 + gtE**2/4)

            # Add contribution to elastic per Equation D.2
            elastic += f*(gnE*cos2phi - 2*(gg + gf)*sinphi2
                          + 2*(E - Eprime)*sin2phi)

            # Add contribution to capture per Equation D.3
            capture += f*gg

            # Add contribution to fission per Equation D.6
            if gf > 0:
                fission += f*gf

    return (elastic, capture, fission)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def reconstruct_rm(rm, double E):
    """Evaluate cross section using Reich-Moore data.

    Parameters
    ----------
    rm : openmc.data.ReichMoore
        Reich-Moore resonance parameters
    E : double
        Energy in eV at which to evaluate the cross section

    Returns
    -------
    elastic : double
        Elastic scattering cross section in barns
    capture : double
        Radiative capture cross section in barns
    fission : double
        Fission cross section in barns

    """
    cdef int i, l, m, n, i_res
    cdef int i_s, num_s, i_J, num_J
    cdef double elastic, capture, fission, total
    cdef double A, k, rho, rhohat, I
    cdef double P, S, phi
    cdef double smin, smax, s, Jmin, Jmax, J, j
    cdef double E_r, gn, gg, gfa, gfb, P_r
    cdef double E_diff, abs_value, gJ
    cdef double Kr, Ki, x
    cdef double complex Ubar, U_, factor
    cdef bint hasfission
    cdef np.ndarray[double, ndim=2] one
    cdef np.ndarray[double complex, ndim=2] K, Imat, U
    cdef double [:,:] params

    # Get nuclear spin
    I = rm.target_spin

    elastic = 0.
    fission = 0.
    total = 0.
    A = rm.atomic_weight_ratio
    k = _wave_number(A, E)
    one = np.eye(3)
    K = np.zeros((3,3), dtype=complex)

    for i, l in enumerate(rm._l_values):
        # Check for l-dependent scattering radius
        rho = k*rm.channel_radius[l](E)
        rhohat = k*rm.scattering_radius[l](E)

        # Calculate shift and penetrability
        P, S = penetration_shift(l, rho)

        # Calculate phase shift
        phi = phaseshift(l, rhohat)

        # Calculate common factor on collision matrix terms (term outside curly
        # braces in ENDF-102, Eq. D.27)
        Ubar = cexp(-2j*phi)

        # The channel spin is the vector sum of the target spin, I, and the
        # neutron spin, 1/2, so can take on values of |I - 1/2| < s < I + 1/2
        smin = abs(I - 0.5)
        smax = I + 0.5
        num_s = int(smax - smin + 1)

        for i_s in range(num_s):
            s = i_s + smin

            # Total angular momentum is the vector sum of l and s and can assume
            # values between |l - s| < J < l + s
            Jmin = abs(l - s)
            Jmax = l + s
            num_J = int(Jmax - Jmin + 1)

            for i_J in range(num_J):
                J = i_J + Jmin

                # Initialize K matrix
                for m in range(3):
                    for n in range(3):
                        K[m,n] = 0.0

                hasfission = False
                if (l, J) in rm._parameter_matrix:
                    params = rm._parameter_matrix[l, J]

                    for i_res in range(params.shape[0]):
                        # Sometimes, the same (l, J) quantum numbers can occur
                        # for different values of the channel spin, s. In this
                        # case, the sign of the channel spin indicates which
                        # spin is to be used. If the spin is negative assume
                        # this resonance comes from the I - 1/2 channel and vice
                        # versa.
                        j = params[i_res, 2]
                        if l > 0:
                            if (j < 0 and s != smin) or (j > 0 and s != smax):
                                continue

                        # Copy resonance parameters
                        E_r = params[i_res, 0]
                        gn = params[i_res, 3]
                        gg = params[i_res, 4]
                        gfa = params[i_res, 5]
                        gfb = params[i_res, 6]
                        P_r = params[i_res, 7]

                        # Calculate neutron width at energy E
                        gn = sqrt(P*gn/P_r)

                        # Calculate j/2 * inverse of denominator of K matrix terms
                        factor = 0.5j/(E_r - E - 0.5j*gg)

                        # Upper triangular portion of K matrix -- see ENDF-102,
                        # Equation D.28
                        K[0,0] = K[0,0] + gn*gn*factor
                        if gfa != 0.0 or gfb != 0.0:
                            # Negate fission widths if necessary
                            gfa = (-1 if gfa < 0 else 1)*sqrt(abs(gfa))
                            gfb = (-1 if gfb < 0 else 1)*sqrt(abs(gfb))

                            K[0,1] = K[0,1] + gn*gfa*factor
                            K[0,2] = K[0,2] + gn*gfb*factor
                            K[1,1] = K[1,1] + gfa*gfa*factor
                            K[1,2] = K[1,2] + gfa*gfb*factor
                            K[2,2] = K[2,2] + gfb*gfb*factor
                            hasfission = True

                # Get collision matrix
                gJ = (2*J + 1)/(4*I + 2)
                if hasfission:
                    # Copy upper triangular portion of K to lower triangular
                    K[1,0] = K[0,1]
                    K[2,0] = K[0,2]
                    K[2,1] = K[1,2]

                    Imat = inv(one - K)
                    U = Ubar*(2*Imat - one)  # ENDF-102, Eq. D.27
                    elastic += gJ*cabs(1 - U[0,0])**2  # ENDF-102, Eq. D.24
                    total += 2*gJ*(1 - creal(U[0,0]))  # ENDF-102, Eq. D.23

                    # Calculate fission from ENDF-102, Eq. D.26
                    fission += 4*gJ*(cabs(Imat[1,0])**2 + cabs(Imat[2,0])**2)
                else:
                    U_ = Ubar*(2/(1 - K[0,0]) - 1)
                    if abs(creal(K[0,0])) < 3e-4 and abs(phi) < 3e-4:
                        # If K and phi are both very small, the calculated cross
                        # sections can lose precision because the real part of U
                        # ends up very close to unity. To get around this, we
                        # use Euler's formula to express Ubar by real and
                        # imaginary parts, expand cos(2phi) = 1 - 2phi^2 +
                        # O(phi^4), and then simplify
                        Kr = creal(K[0,0])
                        Ki = cimag(K[0,0])
                        x = 2*(-Kr + (Kr*Kr + Ki*Ki)*(1 - phi*phi) + phi*phi -
                               sin(2*phi)*Ki)/((1 - Kr)*(1 - Kr) + Ki*Ki)
                        total += 2*gJ*x
                        elastic += gJ*(x*x + cimag(U_)**2)
                    else:
                        total += 2*gJ*(1 - creal(U_))   # ENDF-102, Eq. D.23
                        elastic += gJ*cabs(1 - U_)**2   # ENDF-102, Eq. D.24

    # Calculate capture as difference of other cross sections as per ENDF-102,
    # Equation D.25
    capture = total - elastic - fission

    elastic *= M_PI/(k*k)
    capture *= M_PI/(k*k)
    fission *= M_PI/(k*k)

    return (elastic, capture, fission)
