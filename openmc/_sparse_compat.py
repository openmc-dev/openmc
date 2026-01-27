"""Compatibility module for scipy.sparse arrays

This module provides a compatibility layer for working with scipy.sparse arrays
across different scipy versions. Sparse arrays were introduced gradually in
scipy, with full support arriving in scipy 1.15. This module provides a unified
API that uses sparse arrays when available and falls back to sparse matrices for
older scipy versions.

For more information on the migration from sparse matrices to sparse arrays,
see: https://docs.scipy.org/doc/scipy/reference/sparse.migration_to_sparray.html
"""

import scipy
from scipy import sparse as sp

# Check scipy version for feature availability
_SCIPY_VERSION = tuple(map(int, scipy.__version__.split('.')[:2]))

if _SCIPY_VERSION >= (1, 15):
    # Use sparse arrays
    csr_array = sp.csr_array
    csc_array = sp.csc_array
    dok_array = sp.dok_array
    lil_array = sp.lil_array
    eye_array = sp.eye_array
    block_array = sp.block_array
else:
    # Fall back to sparse matrices
    csr_array = sp.csr_matrix
    csc_array = sp.csc_matrix
    dok_array = sp.dok_matrix
    lil_array = sp.lil_matrix
    eye_array = sp.eye
    block_array = sp.bmat

__all__ = [
    'csr_array',
    'csc_array',
    'dok_array',
    'lil_array',
    'eye_array',
    'block_array',
]
