"""
Fast Relaxed Vector Fitting function

Approximate f(s) with a rational function:
        f(s)=R*(s*I-A)^(-1) + Polynomials*s
where f(s) is a vector of elements.

When f(s) is a vector, all elements become fitted with a common pole set. The
identification is done using the pole relocating method known as Vector Fitting
[1] with relaxed non-triviality constraint for faster convergence and smaller
fitting errors [2], and utilization of matrix structure for fast solution of the
pole identifion step [3].

[1] B. Gustavsen and A. Semlyen, "Rational approximation of frequency
    domain responses by Vector Fitting", IEEE Trans. Power Delivery, vol. 14,
    no. 3, pp. 1052-1061, July 1999.
[2] B. Gustavsen, "Improving the pole relocating properties of vector
    fitting", IEEE Trans. Power Delivery, vol. 21, no. 3, pp. 1587-1592, July
    2006.
[3] D. Deschrijver, M. Mrozowski, T. Dhaene, and D. De Zutter,
    "Macromodeling of Multiport Systems Using a Fast Implementation of the
    Vector Fitting Method", IEEE Microwave and Wireless Components Letters, vol.
    18, no. 6, pp. 383-385, June 2008.

All credit goes to:
 - Bjorn Gustavsen for his MATLAB implementation.
   (http://www.sintef.no/Projectweb/VECTFIT/)
 - Jingang Liang for his C++ implementation.
   (https://github.com/mit-crpg/vectfit.git)

"""

from typing import Tuple

import numpy as np
from scipy.linalg import eigvals, lstsq, norm, qr


def evaluate(
    eval_points: np.ndarray,
    pole_values: np.ndarray,
    residue_matrix: np.ndarray,
    poly_coefficients: np.ndarray | None = None,
) -> np.ndarray:
    """Evaluate the rational function approximation:
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

    Raises
    ------
    ValueError
        If input arrays have incompatible shapes.
    """
    eval_points = np.asarray(eval_points)
    pole_values = np.asarray(pole_values)
    residue_matrix = np.asarray(residue_matrix)

    if eval_points.ndim != 1:
        raise ValueError("eval_points must be a 1D array")
    if pole_values.ndim != 1:
        raise ValueError("pole_values must be a 1D array")

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

    # term: sum over poles of (residues / (eval_points - poles))
    denominator = (
        eval_points[np.newaxis, :] - pole_values[:, np.newaxis]
    )  # shape: (num_poles, num_eval)
    pole_terms = residue_matrix @ (1.0 / denominator)  # shape: (num_vectors, num_eval)
    result = np.real(pole_terms)

    # polynomial part: sum over poly_idx of (coeff * eval_points**poly_idx)
    if num_coeffs > 0:
        powers = (
            eval_points[np.newaxis, :] ** np.arange(num_coeffs)[:, np.newaxis]
        )  # shape: (num_coeffs, num_eval)
        result += poly_coefficients @ powers  # shape: (num_vectors, num_eval)

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
    """Perform vector fitting using the Fast Relaxed Vector Fitting algorithm.

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
        raise ValueError("n_polys must be in [0, 11]")

    residue_matrix = np.zeros((num_vectors, num_poles), dtype=np.complex128)
    poly_coefficients = np.zeros((num_vectors, n_polys))
    fit_result = np.zeros_like(response_matrix)
    rms_error = 0.0

    if num_poles == 0 and n_polys == 0:
        rms_error = norm(response_matrix) / np.sqrt(num_vectors * num_samples)
        return initial_poles, residue_matrix, poly_coefficients, fit_result, rms_error

    if not skip_pole_update and num_poles > 0:
        updated_poles = identify_poles(
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
        fit_result, rms_error = identify_residues(
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


def compute_dk_matrix(
    dk_matrix: np.ndarray,
    eval_points: np.ndarray,
    poles: np.ndarray,
    conj_index: np.ndarray,
    num_poles: int,
    num_polys: int,
    tol_high: float = None,
):
    """Compute the dk_matrix used in windowed multipole evaluations.

    Parameters
    ----------
    dk_matrix : ndarray of shape (len(eval_points), M)
        The full matrix used in least-squares fitting or evaluation.
    eval_points : ndarray of shape (N,)
        Energy points at which to evaluate.
    poles : ndarray of shape (num_poles,)
        Complex poles used in the resonance model.
    conj_index : ndarray of shape (num_poles,)
        Index array indicating pole conjugacy behavior: 0 (normal), 1 (add conjugate), 2 (imaginary part).
    num_poles : int
        Number of complex poles.
    num_polys : int
        Number of polynomial terms (including constant term).
    tol_high : float
        Replacement value for infinities.
    """
    # Broadcast shapes
    eval_points_col = eval_points[:, np.newaxis]
    poles_row = poles[np.newaxis, :]

    # Compute base terms
    term1 = 1.0 / (eval_points_col - poles_row)
    term2 = 1.0 / (eval_points_col - np.conj(poles_row))
    term3 = 1j / (eval_points_col - np.conj(poles_row)) - 1j / (
        eval_points_col - poles_row
    )

    # Masks for different conjugacy types
    mask0 = conj_index == 0
    mask1 = conj_index == 1
    mask2 = conj_index == 2

    # Fill dk_matrix with pole terms
    dk_matrix[:, :num_poles][:, mask0] = term1[:, mask0]
    dk_matrix[:, :num_poles][:, mask1] = term1[:, mask1] + term2[:, mask1]
    dk_matrix[:, :num_poles][:, mask2] = term3[:, mask2]

    # Replace infinities with high tolerance value
    if tol_high is not None:
        inf_mask = np.isinf(dk_matrix)
        dk_matrix[inf_mask] = tol_high + 0j

    # Add polynomial basis (Chebyshev-like, just powers here)
    powers = np.arange(num_polys)
    dk_matrix[:, num_poles : num_poles + num_polys] = eval_points_col**powers + 0j
    return dk_matrix


def row_block_matrix(
    dk_matrix: np.ndarray,
    weights: np.ndarray,
    response_matrix: np.ndarray,
    vec_idx: int,
    num_poles: int,
    num_polys: int,
) -> np.ndarray:
    """
    Construct a single matrix row block for the given vector index.

    Parameters
    ----------
    dk_matrix : ndarray of shape (num_samples, num_poles + num_polys)
        Basis function evaluations at each sample point.
    weights : ndarray of shape (num_vectors, num_samples)
        Sample weights for each vector.
    response_matrix : ndarray of shape (num_vectors, num_samples)
        Response values at each sample point.
    vec_idx : int
        Index of the vector to construct the A1 block for.
    num_poles : int
        Number of poles used in the model.
    num_polys : int
        Number of polynomial basis terms.

    Returns
    -------
    A : ndarray of shape (num_samples, num_poles + num_polys + num_poles + 1)
        Weighted and assembled matrix block for the current vector.
    """
    num_samples = dk_matrix.shape[0]
    A = np.zeros(
        (num_samples, num_poles + num_polys + num_poles + 1), dtype=np.complex128
    )

    # Weighted basis terms
    A[:, : num_poles + num_polys] = (
        weights[vec_idx][:, np.newaxis] * dk_matrix[:, : num_poles + num_polys]
    )

    # Weighted response terms (includes poles + 1)
    A[:, num_poles + num_polys : num_poles + num_polys + num_poles + 1] = (
        -weights[vec_idx][:, np.newaxis]
        * dk_matrix[:, : num_poles + 1]
        * response_matrix[vec_idx][:, np.newaxis]
    )

    return A


def process_constrained_block(
    vec_idx: int,
    dk_matrix: np.ndarray,
    weights: np.ndarray,
    response_matrix: np.ndarray,
    num_samples: int,
    num_poles: int,
    num_polys: int,
    scale_factor: float,
    num_vectors: int,
) -> Tuple[int, np.ndarray, np.ndarray]:
    """
    Construct a constrained least-squares system block for the given vector index.

    This function computes the A matrix using weighted evaluations of the basis functions
    and response terms. It appends a constraint row to enforce physical properties
    (e.g., normalization) **only for the final vector index**. The full matrix A is
    decomposed via QR, and the resulting triangular block is returned.

    This routine is intended for use in the main vector fitting loop when the denominator
    is well-conditioned but requires an additional constraint row for physical consistency.

    Parameters
    ----------
    vec_idx : int
        Index of the vector to process.
    dk_matrix : ndarray of shape (num_samples, num_poles + num_polys)
        Evaluated basis functions at sample points.
    weights : ndarray of shape (num_vectors, num_samples)
        Weight matrix per vector.
    response_matrix : ndarray of shape (num_vectors, num_samples)
        Response function values for each vector.
    num_samples : int
        Number of sample points.
    num_poles : int
        Number of poles in the model.
    num_polys : int
        Number of polynomial terms in the model.
    scale_factor : float
        Scaling factor applied to the final constraint row.
    num_vectors : int
        Total number of vectors to process.

    Returns
    -------
    vec_idx : int
        Index of the processed vector.
    lhs_block : ndarray of shape (num_poles + 1, num_poles + 1)
        Triangular matrix block from QR decomposition.
    rhs_block : ndarray of shape (num_poles + 1,) or None
        Right-hand side vector block (only returned for final vec_idx), else None.
    """
    A1 = row_block_matrix(
        dk_matrix, weights, response_matrix, vec_idx, num_poles, num_polys
    )
    A = np.zeros((2 * num_samples + 1, num_poles + num_polys + num_poles + 1))
    A[:num_samples] = A1.real
    A[num_samples : 2 * num_samples] = A1.imag
    # Handle final row only if vec_idx is last
    if vec_idx == num_vectors - 1:
        A[
            2 * num_samples,
            num_poles + num_polys : num_poles + num_polys + num_poles + 1,
        ] = scale_factor * np.real(dk_matrix[:, : num_poles + 1].sum(axis=0))

    Q, R = qr(A, mode="economic")

    lhs_block = R[
        num_poles + num_polys : num_poles + num_polys + num_poles + 1,
        num_poles + num_polys : num_poles + num_polys + num_poles + 1,
    ]

    if vec_idx == num_vectors - 1:
        rhs_block = (
            num_samples
            * scale_factor
            * Q[-1, num_poles + num_polys : num_poles + num_polys + num_poles + 1]
        )
    else:
        rhs_block = np.zeros_like(
            Q[-1, num_poles + num_polys : num_poles + num_polys + num_poles + 1]
        )

    return vec_idx, lhs_block, rhs_block


def process_unconstrained_block(
    vec_idx: int,
    dk_matrix: np.ndarray,
    weights: np.ndarray,
    response_matrix: np.ndarray,
    denom: float,
    num_poles: int,
    num_polys: int,
) -> Tuple[int, np.ndarray, np.ndarray]:
    """
    Construct an unconstrained least-squares system block for the given vector index.

    This function is used when the fitting denominator becomes ill-conditioned
    (too small or too large), and the original constrained system is replaced by
    an alternative regularized least-squares problem. The A matrix is built by stacking
    the real and imaginary parts of the basis evaluations, and the RHS vector b is
    scaled by `denom`.

    A standard QR decomposition is used to extract the square block of the system,
    which can be solved independently from the constrained system.

    Parameters
    ----------
    vec_idx : int
        Index of the vector to process.
    dk_matrix : ndarray of shape (num_samples, num_poles + num_polys)
        Evaluated basis functions at sample points.
    weights : ndarray of shape (num_vectors, num_samples)
        Weight matrix per vector.
    response_matrix : ndarray of shape (num_vectors, num_samples)
        Response function values for each vector.
    denom : float
        Scaling factor applied to the right-hand side vector b.
    num_poles : int
        Number of poles in the model.
    num_polys : int
        Number of polynomial terms in the model.

    Returns
    -------
    vec_idx : int
        Index of the processed vector.
    lhs_block : ndarray of shape (num_poles, num_poles)
        Triangular matrix block from QR decomposition.
    rhs_block : ndarray of shape (num_poles,)
        Right-hand side vector block for this vector.
    """
    A1 = row_block_matrix(
        dk_matrix, weights, response_matrix, vec_idx, num_poles, num_polys
    )
    A = np.vstack((A1.real, A1.imag))

    b1 = denom * weights[vec_idx] * response_matrix[vec_idx]
    b = np.concatenate((b1.real, b1.imag))

    Q, R = qr(A, mode="economic")

    lhs_block = R[
        num_poles + num_polys : num_poles + num_polys + num_poles,
        num_poles + num_polys : num_poles + num_polys + num_poles,
    ]
    rhs_block = Q[:, num_poles + num_polys : num_poles + num_polys + num_poles].T @ b
    return vec_idx, lhs_block, rhs_block


def identify_poles(
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
    conj_index = label_conjugate_poles(poles)
    dk_matrix = np.zeros(
        (num_samples, num_poles + max(num_polys, 1)), dtype=np.complex128
    )
    compute_dk_matrix(
        dk_matrix, eval_points, poles, conj_index, num_poles, num_polys, tol_high
    )

    scale_factor = (
        np.sqrt(
            sum(norm(weights[m] * response_matrix[m]) ** 2 for m in range(num_vectors))
        )
        / num_samples
    )
    lhs_matrix = np.zeros((num_vectors * (num_poles + 1), num_poles + 1))
    rhs_vector = np.zeros(num_vectors * (num_poles + 1))

    for vec_idx in range(num_vectors):
        vec_idx, lhs_block, rhs_block = process_constrained_block(
            vec_idx,
            dk_matrix,
            weights,
            response_matrix,
            num_samples,
            num_poles,
            num_polys,
            scale_factor,
            num_vectors,
        )
        i0 = vec_idx * (num_poles + 1)
        i1 = (vec_idx + 1) * (num_poles + 1)
        lhs_matrix[i0:i1] = lhs_block
        rhs_vector[i0:i1] = rhs_block

    column_scale = np.linalg.norm(lhs_matrix, axis=0)
    column_scale[column_scale == 0] = 1.0
    lhs_matrix *= column_scale

    solution, *_ = lstsq(lhs_matrix, rhs_vector)
    solution *= column_scale
    coeffs = solution[:-1]
    denom = solution[-1]

    if abs(denom) < tol_low or abs(denom) > tol_high:
        lhs_matrix = np.zeros((num_vectors * num_poles, num_poles))
        rhs_vector = np.zeros(num_vectors * num_poles)
        # Adjust denom
        if denom == 0.0:
            denom = 1.0
        elif abs(denom) < tol_low:
            denom = np.sign(denom) * tol_low
        elif abs(denom) > tol_high:
            denom = np.sign(denom) * tol_high

        # Allocate output
        lhs_matrix = np.zeros((num_vectors * num_poles, num_poles))
        rhs_vector = np.zeros(num_vectors * num_poles)

        for vec_idx in range(num_vectors):
            vec_idx, lhs_block, rhs_block = process_unconstrained_block(
                vec_idx,
                dk_matrix,
                weights,
                response_matrix,
                denom,
                num_poles,
                num_polys,
            )
            i0 = vec_idx * num_poles
            i1 = (vec_idx + 1) * num_poles
            lhs_matrix[i0:i1] = lhs_block
            rhs_vector[i0:i1] = rhs_block

        column_scale = np.linalg.norm(lhs_matrix, axis=0)
        column_scale[column_scale == 0] = 1.0
        lhs_matrix *= column_scale
        coeffs, *_ = lstsq(lhs_matrix, rhs_vector)
        coeffs *= column_scale

    lambda_matrix = np.zeros((num_poles, num_poles))
    scale_vector = np.ones((num_poles, 1))

    # Mask for real poles (conj_index == 0)
    mask_real = conj_index == 0
    real_indices = np.where(mask_real)[0]
    lambda_matrix[real_indices, real_indices] = np.real(poles[real_indices])

    # Mask for start of complex conjugate pairs (conj_index == 1)
    mask_cplx_start = conj_index == 1
    cplx_indices = np.where(mask_cplx_start)[0]

    # Extract real and imaginary parts of complex conjugate poles
    real_parts = np.real(poles[cplx_indices])
    imag_parts = np.imag(poles[cplx_indices])

    # Diagonal assignments
    lambda_matrix[cplx_indices, cplx_indices] = real_parts
    lambda_matrix[cplx_indices + 1, cplx_indices + 1] = real_parts

    # Off-diagonal assignments
    lambda_matrix[cplx_indices, cplx_indices + 1] = imag_parts
    lambda_matrix[cplx_indices + 1, cplx_indices] = -imag_parts

    # Scaling vector adjustments
    scale_vector[cplx_indices, 0] = 2.0
    scale_vector[cplx_indices + 1, 0] = 0.0

    residue_matrix = lambda_matrix - np.outer(scale_vector.squeeze(), coeffs) / denom
    return eigvals(residue_matrix)


def solve_vector_block(
    vec_idx: int,
    dk_matrix: np.ndarray,
    weights: np.ndarray,
    response_matrix: np.ndarray,
    num_poles: int,
    num_polys: int,
) -> Tuple[int, np.ndarray, np.ndarray]:
    """
    Solve the least-squares system for a single vector index.

    Parameters
    ----------
    vec_idx : int
        Index of the vector to solve.
    dk_matrix : ndarray
        Basis function evaluations of shape (num_samples, num_poles + num_polys).
    weights : ndarray
        Weight array of shape (num_vectors, num_samples).
    response_matrix : ndarray
        Response array of shape (num_vectors, num_samples).
    num_poles : int
        Number of poles.
    num_polys : int
        Number of polynomial coefficients.

    Returns
    -------
    vec_idx : int
        The index of the solved vector.
    residues : ndarray
        Solution vector for the residues (length = num_poles).
    poly_coeffs : ndarray or None
        Solution vector for polynomial coefficients (length = num_polys), or None if num_polys == 0.
    """
    A = dk_matrix * weights[vec_idx][:, np.newaxis]
    b = weights[vec_idx] * response_matrix[vec_idx]

    lhs_matrix = np.vstack((A.real, A.imag))
    rhs_vector = np.concatenate((b.real, b.imag))

    scale_column = np.linalg.norm(lhs_matrix, axis=0)
    scale_column[scale_column == 0] = 1.0
    lhs_matrix *= scale_column
    x = lstsq(lhs_matrix, rhs_vector)[0]
    x *= scale_column

    residues = x[:num_poles]
    poly_coeffs = x[num_poles : num_poles + num_polys] if num_polys > 0 else None

    return vec_idx, residues, poly_coeffs


def identify_residues(
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
    conj_index = label_conjugate_poles(poles)
    dk_matrix = np.zeros((num_samples, num_poles + num_polys), dtype=np.complex128)

    compute_dk_matrix(dk_matrix, eval_points, poles, conj_index, num_poles, num_polys)

    real_residues = np.zeros((num_vectors, num_poles), dtype=np.float64)

    for vec_idx in range(num_vectors):
        vec_idx, residues, poly_coeffs = solve_vector_block(
            vec_idx,
            dk_matrix,
            weights,
            response_matrix,
            num_poles,
            num_polys,
        )
        real_residues[vec_idx] = residues
        if poly_coeffs is not None:
            poly_coefficients[vec_idx] = poly_coeffs

    # Mask for real poles
    mask_real = conj_index == 0
    real_indices = np.where(mask_real)[0]
    residue_matrix[:, real_indices] = real_residues[:, real_indices]

    # Mask for first of complex conjugate pairs
    mask_cplx_start = conj_index == 1
    cplx_indices = np.where(mask_cplx_start)[0]

    # Compute complex residues using vectorized operations
    residue_matrix[:, cplx_indices] = (
        real_residues[:, cplx_indices] + 1j * real_residues[:, cplx_indices + 1]
    )
    residue_matrix[:, cplx_indices + 1] = (
        real_residues[:, cplx_indices] - 1j * real_residues[:, cplx_indices + 1]
    )

    fit_result = evaluate(eval_points, poles, residue_matrix, poly_coefficients)
    rms_error = norm(fit_result - response_matrix) / np.sqrt(num_vectors * num_samples)
    return fit_result, rms_error


def label_conjugate_poles(poles: np.ndarray) -> np.ndarray:
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

    # Identify complex poles (nonzero imaginary part)
    is_complex = np.imag(poles) != 0.0

    # Find conjugate pairs: poles[i+1] ≈ conj(poles[i])
    is_pair_start = is_complex[:-1] & np.isclose(np.conj(poles[:-1]), poles[1:])

    # Mark valid conjugate pair entries
    conj_index[:-1][is_pair_start] = 1  # mark i with 1
    conj_index[1:][is_pair_start] = 2  # mark i+1 with 2

    # Now validate: all complex poles must be part of valid conjugate pairs
    unmatched_complex = is_complex & (conj_index == 0)
    if np.any(unmatched_complex):
        raise ValueError("Complex poles must appear in conjugate pairs")

    return conj_index
