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
_HAS_SPARSE_ARRAYS = _SCIPY_VERSION >= (1, 15)

# Define which constructors and format checkers to use
if _HAS_SPARSE_ARRAYS:
    # Use sparse arrays (scipy 1.12+)
    csc_array = sp.csc_array
    csc_matrix = sp.csc_array  # Alias for backward compatibility
    dok_array = sp.dok_array
    dok_matrix = sp.dok_array  # Alias for backward compatibility

    def eye_array(n, **kwargs):
        """Create an identity sparse array"""
        return sp.eye_array(n, **kwargs)

    def block_array(blocks, **kwargs):
        """Create a block sparse array from a 2D array of sparse arrays/matrices"""
        return sp.block_array(blocks, **kwargs)

    def hstack(blocks, **kwargs):
        """Stack sparse arrays horizontally (column wise)"""
        # Ensure at least one input is a sparse array to get array output
        if blocks and not any(isinstance(b, sp.sparray) for b in blocks):
            blocks = list(blocks)
            if len(blocks) > 0 and hasattr(sp, 'csc_array'):
                blocks[0] = sp.csc_array(blocks[0])
        return sp.hstack(blocks, **kwargs)

    def vstack(blocks, **kwargs):
        """Stack sparse arrays vertically (row wise)"""
        # Ensure at least one input is a sparse array to get array output
        if blocks and not any(isinstance(b, sp.sparray) for b in blocks):
            blocks = list(blocks)
            if len(blocks) > 0 and hasattr(sp, 'csc_array'):
                blocks[0] = sp.csc_array(blocks[0])
        return sp.vstack(blocks, **kwargs)

    def is_sparse(A):
        """Check if A is a sparse array or matrix"""
        return sp.issparse(A)

    def is_sparse_format(A, format_str):
        """Check if A is a sparse array/matrix with the given format"""
        return sp.issparse(A) and A.format == format_str

else:
    # Fall back to sparse matrices (scipy < 1.12)
    csc_array = sp.csc_matrix
    csc_matrix = sp.csc_matrix
    dok_array = sp.dok_matrix
    dok_matrix = sp.dok_matrix

    def eye_array(n, **kwargs):
        """Create an identity sparse matrix"""
        return sp.eye(n, **kwargs)

    def block_array(blocks, **kwargs):
        """Create a block sparse matrix from a 2D array of sparse matrices"""
        return sp.bmat(blocks, **kwargs)

    def hstack(blocks, **kwargs):
        """Stack sparse matrices horizontally (column wise)"""
        return sp.hstack(blocks, **kwargs)

    def vstack(blocks, **kwargs):
        """Stack sparse matrices vertically (row wise)"""
        return sp.vstack(blocks, **kwargs)

    def is_sparse(A):
        """Check if A is a sparse matrix"""
        return sp.issparse(A)

    def is_sparse_format(A, format_str):
        """Check if A is a sparse matrix with the given format"""
        return sp.issparse(A) and A.format == format_str


__all__ = [
    'csc_array',
    'csc_matrix',
    'dok_array',
    'dok_matrix',
    'eye_array',
    'block_array',
    'hstack',
    'vstack',
    'is_sparse',
    'is_sparse_format',
]
