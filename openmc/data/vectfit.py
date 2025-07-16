"""
Fast Relaxed Vector Fitting function

Approximate f(s) with a rational function:
        f(s)=R*(s*I-A)^(-1) + Polynomials*s
where f(s) is a vector of elements.

When f(s) is a vector, all elements become fitted with a common pole set.
The identification is done using the pole relocating method known as Vector
Fitting [1] with relaxed non-triviality constraint for faster convergence
and smaller fitting errors [2], and utilization of matrix structure for fast
solution of the pole identifion step [3].

[1] B. Gustavsen and A. Semlyen, "Rational approximation of frequency
    domain responses by Vector Fitting", IEEE Trans. Power Delivery,
    vol. 14, no. 3, pp. 1052-1061, July 1999.
[2] B. Gustavsen, "Improving the pole relocating properties of vector
    fitting", IEEE Trans. Power Delivery, vol. 21, no. 3, pp. 1587-1592,
    July 2006.
[3] D. Deschrijver, M. Mrozowski, T. Dhaene, and D. De Zutter,
    "Macromodeling of Multiport Systems Using a Fast Implementation of
    the Vector Fitting Method", IEEE Microwave and Wireless Components
    Letters, vol. 18, no. 6, pp. 383-385, June 2008.

All credit goes to:
- Bjorn Gustavsen for his MATLAB implementation.
(http://www.sintef.no/Projectweb/VECTFIT/)
- Jingang Liang for his C++ implementation.
(https://github.com/mit-crpg/vectfit.git)

"""
import numpy as np
from scipy.linalg import lstsq, eigvals, norm, qr
from typing import Tuple


def evaluate(
    eval_points: np.ndarray,
    pole_values: np.ndarray,
    residue_matrix: np.ndarray,
    poly_coefficients: np.ndarray = None,
) -> np.ndarray:
    """
    Evaluate the rational function approximation:
        f(s) ≈ sum(residue / (s - pole)) + sum(poly_coefficients * s^j)

    Parameters
    ----------
    eval_points : np.ndarray
        1D array of real scalar frequency values (s).
    pole_values : np.ndarray
        1D array of complex poles.
    residue_matrix : np.ndarray
        2D or 1D array of complex residues (shape: [num_vectors, num_poles] or [num_poles]).
    poly_coefficients : np.ndarray, optional
        2D or 1D array of real polynomial coefficients (shape: [num_vectors, num_polys]).

    Returns
    -------
    np.ndarray
        2D array of evaluated real function values (shape: [num_vectors, num_samples]).
    """
    eval_points = np.asarray(eval_points)
    pole_values = np.asarray(pole_values)
    residue_matrix = np.asarray(residue_matrix)

    if residue_matrix.ndim == 1:
        residue_matrix = residue_matrix.reshape((1, -1))
    num_vectors, _ = residue_matrix.shape
    num_samples = len(eval_points)

    if poly_coefficients is not None and isinstance(poly_coefficients, list):
        poly_coefficients = np.array(poly_coefficients)
    if poly_coefficients is None or poly_coefficients.size == 0:
        poly_coefficients = np.zeros((num_vectors, 0))
    else:
        poly_coefficients = np.asarray(poly_coefficients)
        if poly_coefficients.ndim == 1:
            poly_coefficients = poly_coefficients.reshape((1, -1))
        elif poly_coefficients.shape[0] != num_vectors:
            raise ValueError("Mismatch in residues and poly_coefficients shapes")

    num_coeffs = poly_coefficients.shape[1]
    result = np.zeros((num_vectors, num_samples))

    for vec_idx in range(num_vectors):
        term = np.sum(
            residue_matrix[vec_idx][:, None] / (eval_points - pole_values[:, None]),
            axis=0,
        )
        result[vec_idx] = np.real(term)
        for poly_idx in range(num_coeffs):
            result[vec_idx] += (
                poly_coefficients[vec_idx, poly_idx] * eval_points**poly_idx
            )

    return result


def vectfit(
    response_matrix: np.ndarray,
    eval_points: np.ndarray,
    initial_poles: np.ndarray,
    weights: np.ndarray,
    n_polys: int = 0,
    skip_pole_update: bool = False,
    skip_residue_update: bool = False,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    """
    Perform vector fitting using the Fast Relaxed Vector Fitting algorithm.

    Parameters
    ----------
    response_matrix : np.ndarray
        Complex matrix of frequency responses (shape: [num_vectors, num_samples]).
    eval_points : np.ndarray
        Real frequency samples (s), shape (num_samples,).
    initial_poles : np.ndarray
        Initial guess for poles (complex), shape (num_poles,).
    weights : np.ndarray
        Weighting matrix for fitting (same shape as response_matrix).
    n_polys : int, optional
        Number of real polynomial terms to include.
    skip_pole_update : bool, optional
        Whether to skip pole relocation step.
    skip_residue_update : bool, optional
        Whether to skip residue fitting step.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]
        - Updated poles (np.ndarray)
        - Residues (np.ndarray)
        - Polynomial coefficients (np.ndarray)
        - Fitted response matrix (np.ndarray)
        - Root-mean-square error (float)
    """
    tol_low = 1e-18
    tol_high = 1e18
    response_matrix = np.asarray(response_matrix)
    eval_points = np.asarray(eval_points)
    initial_poles = np.asarray(initial_poles)
    weights = np.asarray(weights)

    num_vectors, num_samples = response_matrix.shape
    num_poles = len(initial_poles)

    if n_polys < 0 or n_polys > 11:
        raise ValueError("num_polynomials must be in [0, 11]")

    residue_matrix = np.zeros((num_vectors, num_poles), dtype=complex)
    poly_coefficients = np.zeros((num_vectors, n_polys))
    fit_result = np.zeros_like(response_matrix)
    rms_error = 0.0

    if num_poles == 0 and n_polys == 0:
        rms_error = norm(response_matrix) / np.sqrt(num_vectors * num_samples)
        return initial_poles, residue_matrix, poly_coefficients, fit_result, rms_error

    if not skip_pole_update and num_poles > 0:
        updated_poles = _identify_poles(
            num_poles,
            num_samples,
            n_polys,
            initial_poles,
            eval_points,
            tol_high,
            weights,
            response_matrix,
            num_vectors,
            tol_low,
        )
    else:
        updated_poles = initial_poles

    if not skip_residue_update:
        fit_result, rms_error = _identify_residues(
            num_poles,
            updated_poles,
            num_samples,
            n_polys,
            eval_points,
            num_vectors,
            weights,
            response_matrix,
            poly_coefficients,
            residue_matrix,
        )

    return updated_poles, residue_matrix, poly_coefficients, fit_result, rms_error


def _identify_poles(
    num_poles: int,
    num_samples: int,
    num_polys: int,
    poles: np.ndarray,
    eval_points: np.ndarray,
    tol_high: float,
    weights: np.ndarray,
    response_matrix: np.ndarray,
    num_vectors: int,
    tol_low: float,
) -> np.ndarray:
    """
    Internal routine to update poles via relaxed vector fitting.

    Parameters
    ----------
    num_poles : int
        Number of poles.
    num_samples : int
        Number of frequency samples.
    num_polys : int
        Number of polynomial terms.
    poles : np.ndarray
        Initial poles (complex), shape (num_poles,).
    eval_points : np.ndarray
        Real frequency values, shape (num_samples,).
    tol_high : float
        Upper tolerance threshold for denominator.
    weights : np.ndarray
        Weighting matrix, shape (num_vectors, num_samples).
    response_matrix : np.ndarray
        Complex frequency responses, shape (num_vectors, num_samples).
    num_vectors : int
        Number of response vectors.
    tol_low : float
        Lower tolerance threshold for denominator.

    Returns
    -------
    np.ndarray
        Updated poles as eigenvalues (shape: [num_poles]).
    """
    conj_index = _label_conjugate_poles(poles)
    
    dk_matrix = np.zeros((num_samples, num_poles + max(num_polys, 1)), dtype=complex)
    for m in range(num_poles):
        p = poles[m]
        if conj_index[m] == 0:
            dk_matrix[:, m] = 1.0 / (eval_points - p)
        elif conj_index[m] == 1:
            dk_matrix[:, m] = 1.0 / (eval_points - p) + 1.0 / (eval_points - np.conj(p))
        elif conj_index[m] == 2:
            dk_matrix[:, m] = 1j / (eval_points - np.conj(p)) - 1j / (eval_points - p)

    dk_matrix[np.isinf(dk_matrix)] = tol_high + 0j
    dk_matrix[:, num_poles] = 1.0 + 0j
    for m in range(1, num_polys):
        dk_matrix[:, num_poles + m] = eval_points**m + 0j

    scale_factor = (
        np.sqrt(
            sum(norm(weights[m] * response_matrix[m]) ** 2 for m in range(num_vectors))
        )
        / num_samples
    )
    lhs_matrix = np.zeros((num_vectors * (num_poles + 1), num_poles + 1))
    rhs_vector = np.zeros(num_vectors * (num_poles + 1))

    for vec_idx in range(num_vectors):
        A1 = np.zeros(
            (num_samples, num_poles + num_polys + num_poles + 1), dtype=complex
        )
        for m in range(num_poles + num_polys):
            A1[:, m] = weights[vec_idx] * dk_matrix[:, m]
        for m in range(num_poles + 1):
            A1[:, num_poles + num_polys + m] = (
                -weights[vec_idx] * dk_matrix[:, m] * response_matrix[vec_idx]
            )

        A = np.zeros((2 * num_samples + 1, num_poles + num_polys + num_poles + 1))
        A[:num_samples] = A1.real
        A[num_samples : 2 * num_samples] = A1.imag

        if vec_idx == num_vectors - 1:
            for m in range(num_poles + 1):
                A[2 * num_samples, num_poles + num_polys + m] = scale_factor * np.real(
                    dk_matrix[:, m].sum()
                )

        Q, R = qr(A, mode="economic")
        lhs_matrix[vec_idx * (num_poles + 1) : (vec_idx + 1) * (num_poles + 1)] = R[
            num_poles + num_polys : num_poles + num_polys + num_poles + 1,
            num_poles + num_polys : num_poles + num_polys + num_poles + 1,
        ]
        if vec_idx == num_vectors - 1:
            rhs_vector[vec_idx * (num_poles + 1) : (vec_idx + 1) * (num_poles + 1)] = (
                num_samples
                * scale_factor
                * Q[-1, num_poles + num_polys : num_poles + num_polys + num_poles + 1]
            )

    column_scale = 1.0 / np.linalg.norm(lhs_matrix, axis=0)
    lhs_matrix *= column_scale
    solution, *_ = lstsq(lhs_matrix, rhs_vector)
    solution *= column_scale
    coeffs = solution[:-1]
    denom = solution[-1]

    if abs(denom) < tol_low or abs(denom) > tol_high:
        lhs_matrix = np.zeros((num_vectors * num_poles, num_poles))
        rhs_vector = np.zeros(num_vectors * num_poles)
        if denom == 0.0:
            denom = 1.0
        elif abs(denom) < tol_low:
            denom = np.sign(denom) * tol_low
        elif abs(denom) > tol_high:
            denom = np.sign(denom) * tol_high

        for vec_idx in range(num_vectors):
            A1 = np.zeros(
                (num_samples, num_poles + num_polys + num_poles), dtype=complex
            )
            for m in range(num_poles + num_polys):
                A1[:, m] = weights[vec_idx] * dk_matrix[:, m]
            for m in range(num_poles):
                A1[:, num_poles + num_polys + m] = (
                    -weights[vec_idx] * dk_matrix[:, m] * response_matrix[vec_idx]
                )
            A = np.vstack((A1.real, A1.imag))
            b1 = denom * weights[vec_idx] * response_matrix[vec_idx]
            b = np.concatenate((b1.real, b1.imag))
            Q, R = qr(A, mode="economic")
            lhs_matrix[vec_idx * num_poles : (vec_idx + 1) * num_poles] = R[
                num_poles + num_polys : num_poles + num_polys + num_poles,
                num_poles + num_polys : num_poles + num_polys + num_poles,
            ]
            rhs_vector[vec_idx * num_poles : (vec_idx + 1) * num_poles] = (
                Q[:, num_poles + num_polys : num_poles + num_polys + num_poles].T @ b
            )

        column_scale = 1.0 / np.linalg.norm(lhs_matrix, axis=0)
        lhs_matrix *= column_scale
        coeffs, *_ = lstsq(lhs_matrix, rhs_vector)
        coeffs *= column_scale

    lambda_matrix = np.zeros((num_poles, num_poles))
    scale_vector = np.ones((num_poles, 1))
    for m in range(num_poles):
        if conj_index[m] == 0:
            lambda_matrix[m, m] = np.real(poles[m])
        elif conj_index[m] == 1:
            real_part = np.real(poles[m])
            imag_part = np.imag(poles[m])
            lambda_matrix[m, m] = real_part
            lambda_matrix[m + 1, m + 1] = real_part
            lambda_matrix[m + 1, m] = -imag_part
            lambda_matrix[m, m + 1] = imag_part
            scale_vector[m, 0] = 2.0
            scale_vector[m + 1, 0] = 0.0

    residue_matrix = lambda_matrix - (scale_vector @ coeffs.reshape(1, -1)) / denom
    return eigvals(residue_matrix)


def _identify_residues(
    num_poles: int,
    poles: np.ndarray,
    num_samples: int,
    num_polys: int,
    eval_points: np.ndarray,
    num_vectors: int,
    weights: np.ndarray,
    response_matrix: np.ndarray,
    poly_coefficients: np.ndarray,
    residue_matrix: np.ndarray,
) -> Tuple[np.ndarray, float]:
    """
    Internal routine to compute residues and polynomial coefficients.

    Parameters
    ----------
    num_poles : int
        Number of poles.
    poles : np.ndarray
        Current poles (complex), shape (num_poles,).
    num_samples : int
        Number of frequency samples.
    num_polys : int
        Number of polynomial terms.
    eval_points : np.ndarray
        Real frequency values, shape (num_samples,).
    num_vectors : int
        Number of response vectors.
    weights : np.ndarray
        Weighting matrix, shape (num_vectors, num_samples).
    response_matrix : np.ndarray
        Complex frequency responses, shape (num_vectors, num_samples).
    poly_coefficients : np.ndarray
        Array to store output polynomial coefficients (in-place).
    residue_matrix : np.ndarray
        Array to store output residues (in-place).

    Returns
    -------
    Tuple[np.ndarray, float]
        - Fitted response matrix (np.ndarray)
        - Root-mean-square fitting error (float)
    """
    conj_index = _label_conjugate_poles(poles)

    dk_matrix = np.zeros((num_samples, num_poles + num_polys), dtype=complex)
    for m in range(num_poles):
        p = poles[m]
        if conj_index[m] == 0:
            dk_matrix[:, m] = 1.0 / (eval_points - p)
        elif conj_index[m] == 1:
            dk_matrix[:, m] = 1.0 / (eval_points - p) + 1.0 / (eval_points - np.conj(p))
        elif conj_index[m] == 2:
            dk_matrix[:, m] = 1j / (eval_points - np.conj(p)) - 1j / (eval_points - p)

    for m in range(num_polys):
        dk_matrix[:, num_poles + m] = eval_points**m + 0j

    real_residues = np.zeros((num_vectors, num_poles))
    for vec_idx in range(num_vectors):
        A = dk_matrix * weights[vec_idx][:, None]
        b = weights[vec_idx] * response_matrix[vec_idx]
        lhs_matrix = np.vstack((A.real, A.imag))
        rhs_vector = np.concatenate((b.real, b.imag))
        scale_column = 1.0 / np.linalg.norm(lhs_matrix, axis=0)
        lhs_matrix *= scale_column
        x, *_ = lstsq(lhs_matrix, rhs_vector)
        x *= scale_column
        real_residues[vec_idx] = x[:num_poles]
        if num_polys > 0:
            poly_coefficients[vec_idx] = x[num_poles : num_poles + num_polys]

    for m in range(num_poles):
        if conj_index[m] == 0:
            residue_matrix[:, m] = real_residues[:, m]
        elif conj_index[m] == 1:
            residue_matrix[:, m] = real_residues[:, m] + 1j * real_residues[:, m + 1]
            residue_matrix[:, m + 1] = (
                real_residues[:, m] - 1j * real_residues[:, m + 1]
            )

    fit_result = evaluate(eval_points, poles, residue_matrix, poly_coefficients)
    rms_error = norm(fit_result - response_matrix) / np.sqrt(num_vectors * num_samples)
    return fit_result, rms_error

def _label_conjugate_poles(poles: np.ndarray) -> np.ndarray:
    """
    Ensure complex poles appear in conjugate pairs and label them accordingly.

    Parameters
    ----------
    poles : np.ndarray
        1D array of complex poles.

    Returns
    -------
    np.ndarray
        Array of integers indicating pole type:
        - 0: real pole
        - 1: first in a complex-conjugate pair
        - 2: second in a complex-conjugate pair

    Raises
    ------
    ValueError
        If any complex pole does not have a valid conjugate pair.
    """
    num_poles = len(poles)
    conj_index = np.zeros(num_poles, dtype=int)

    m = 0
    while m < num_poles:
        if np.imag(poles[m]) != 0.0:
            if m == 0 or conj_index[m - 1] in [0, 2]:
                if m >= num_poles - 1 or not np.isclose(np.conj(poles[m]), poles[m + 1]):
                    raise ValueError("Complex poles must appear in conjugate pairs")
                conj_index[m] = 1
                conj_index[m + 1] = 2
                m += 1
        m += 1

    return conj_index

