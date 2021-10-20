from math import gcd
import scipy.sparse as sps
from scipy.sparse.linalg import spsolve

import numpy as np
from numpy.lib.stride_tricks import as_strided as ast

def lcm(a, b):
    """Find the least common multiple of two numbers.

    Parameters
    ----------
    a : int
        First number for the least common multiple calculation.
    b : int
        Second number for the least common multiple calculation.

    """

    return a * b // gcd(a, b)

def block_diag(array):
    """Provide a sparse block diagonal matrix.

    Parameters
    ----------
    array : np.array
        Array to be changed to block diagonal.

    """

    ni, ng, ng = array.shape
    diags = np.zeros((ng*2-1, ni * ng))
    ndiag = [0] + [g for g in range(1,ng)] + [-g for g in range(1,ng)]

    for i in range(ni):
        for r in range(ng):
            for c in range(ng):
                diags[c-r][i*ng+r] = array[i, r, c]

    diags2 = [diags[0]]
    for g in range(1,ng):
        diags2.append(diags[g][:-g])
    for g in range(1,ng):
        diags2.append(diags[-g][g:])

    return sps.diags(diags2, ndiag)

def map_array(array, from_shape, to_shape, normalize=True, lcm_applied=False):
    """Map an array from one shape to another.

    Parameters
    ----------
    array : np.array
        Array to be mapped.
    from_shape : tuple
        Starting shape.
    to_shape : tuple
        Final shape.
    normalize : bool
        Whether to normalize the array.
    lcm_applied : bool
        Whether the least common multiple is applied.

    """

    array.shape = from_shape
    num_dims = len(from_shape)
    if len(to_shape) != num_dims:
        msg = 'from_shape and to_shape have different dimension!'
        raise ValueError(msg)

    if not lcm_applied:
        lcm_shape = []
        for d in range(num_dims):
            lcm_shape.append(lcm(from_shape[d], to_shape[d]))
        lcm_shape = tuple(lcm_shape)

        # Map the array to the lcm mesh
        array = map_array(array, from_shape, lcm_shape, normalize, True)
        from_shape = lcm_shape

    # loop over dimensions
    for d in range(num_dims):

        # condense along dimension
        if from_shape[d] > to_shape[d]:
            ratio = int(from_shape[d] / to_shape[d])
            block = (1,) * d + (ratio,) + (1,) * (num_dims-1-d)
            axes = tuple(range(num_dims,2*num_dims))
            if normalize:
                array = block_view(array, block).mean(axis=axes)
            else:
                array = block_view(array, block).sum(axis=axes)

        # expand along dimension
        else:
            ratio = to_shape[d] / from_shape[d]
            if normalize:
                array = np.repeat(array, ratio, axis=d)
            else:
                array = np.repeat(array, ratio, axis=d) / float(ratio)

    return array

def surface_integral(array, from_shape, to_shape):
    """Perform a surface integral.

    Parameters
    ----------
    array : np.array
        Array to be integrated.
    from_shape : tuple
        Starting shape.
    to_shape : tuple
        Final shape.

    """

    array.shape = from_shape
    num_dims = len(from_shape)
    if len(to_shape) != num_dims:
        msg = 'from_shape and to_shape have different dimension!'
        raise ValueError(msg)

    to_array = np.zeros(to_shape)
    ratios = [int(a/b) for a,b in zip(from_shape,to_shape)]

    # Sum x-min currents
    view = array[:,:,0:from_shape[2]:ratios[2],:,0:2]
    vs = view.shape
    view.shape = (to_shape[0], ratios[0], to_shape[1], ratios[1],vs[2],vs[3],vs[4])
    to_array[...,0:2] = np.sum(view, axis=(1,3))

    # Sum x-max currents
    view = array[:,:,ratios[2]-1:from_shape[2]:ratios[2],:,2:4]
    vs = view.shape
    view.shape = (to_shape[0], ratios[0], to_shape[1], ratios[1], vs[2],vs[3],vs[4])
    to_array[...,2:4] = np.sum(view, axis=(1,3))

    # Sum y-min currents
    view = array[:,0:from_shape[1]:ratios[1],:,:,4:6]
    vs = view.shape
    view.shape = (to_shape[0], ratios[0], vs[1], to_shape[2], ratios[2],vs[3],vs[4])
    to_array[...,4:6] = np.sum(view, axis=(1,4))

    # Sum y-max currents
    view = array[:,ratios[1]-1:from_shape[1]:ratios[1],:,:,6:8]
    vs = view.shape
    view.shape = (to_shape[0], ratios[0], vs[1], to_shape[2], ratios[2],vs[3],vs[4])
    to_array[...,6:8] = np.sum(view, axis=(1,4))

    # Sum z-min currents
    view = array[0:from_shape[0]:ratios[0],:,:,:,8:10]
    vs = view.shape
    view.shape = (vs[0], to_shape[1], ratios[1], to_shape[2], ratios[2],vs[3],vs[4])
    to_array[...,8:10] = np.sum(view, axis=(2,4))

    # Sum z-max currents
    view = array[ratios[0]-1:from_shape[0]:ratios[0],:,:,:,10:12]
    vs = view.shape
    view.shape = (vs[0], to_shape[1], ratios[1], to_shape[2], ratios[2],vs[3],vs[4])
    to_array[...,10:12] = np.sum(view, axis=(2,4))

    to_array.shape = to_shape
    return to_array


def block_view(A, block):
    """Provide a block view to 3+D array. No error checking made.
    Therefore meaningful (as implemented) only for blocks strictly
    compatible with the shape of A.

    Parameters
    ----------
    A : np.array
        Array to be viewed.
    block : tuple
        Block size.

    """

    block_dim = len(block)

    # Check that blocks fit nicely into first dimensions of A
    if block_dim > len(A.shape):
        msg = 'Cannot access block_view of {}-D array with block size {}-D'\
              ' that is greater than array dimensions'.format(len(A.shape),
                                                              block_dim)
        raise ValueError(msg)
    else:
        block = block + A.shape[block_dim:]
        block_dim = len(A.shape)

    for i,b in enumerate(block):
        if A.shape[i] % b != 0:
            msg = 'Cannot access block_view of {}-D array with block {}'\
                  ' since block dimension {} does not fit evenly in'\
                  ' Array dimension {}'.format(len(A.shape), block,
                                               b, A.shape[i])
            raise ValueError(msg)

    # simple shape and strides computations may seem at first strange
    # unless one is able to recognize the 'tuple additions' involved
    shape = ()
    strides = ()

    for i,b in enumerate(block):
        shape = shape + (int(A.shape[i] / block[i]),)
        strides = strides + (b * A.strides[i],)

    shape = shape + block
    strides = strides + A.strides

    return ast(A, shape=shape, strides=strides)

def nan_inf_to_zero(array):
    """Replace inf with zeros in an array.

    Parameters
    ----------
    array : np.array
        Array to undergo inf to zero replacement.

    """

    array[array == -np.inf] = 0.
    array[array ==  np.inf] = 0.
    return np.nan_to_num(array)

def compute_eigenvalue(A, M, flux, tolerance=1.e-6, max_eig_iterations=10000):
    """Compute k eigenvalue using source iteration.

    Parameters
    ----------
    A : np.array
        Destruction matrix.
    M : np.array
        Production matrix.
    flux : np.array
        Initial guess of flux.
    tolerance : float
        Maximum amount of error allowed before ceasing iteration.
    max_eig_iterations : int
        Maximum number of iterations.

    """

    # Ensure flux is a 1D array
    flux = flux.flatten()

    # Compute the initial source
    old_source = M * flux
    norm = old_source.mean()
    old_source  = old_source / norm
    flux  = flux / norm
    k_eff = 1.0

    print('initial k {}'.format((M*flux).sum() / (A*flux).sum()))

    for i in range(max_eig_iterations):

        # Solve linear system
        flux = spsolve(A, old_source)

        # Compute new source
        new_source = M * flux

        # Compute and set k-eff
        k_eff = new_source.mean()

        # Scale the new source by 1 / k-eff
        new_source  = new_source / k_eff

        # Compute the fractional solution change
        nonzero = new_source != 0
        frac_change_array = (new_source[nonzero]-old_source[nonzero])/new_source[nonzero]
        frac_change_array = np.square(frac_change_array)
        frac_change = np.sqrt(frac_change_array.mean())

        # Copy new source to old source
        old_source = np.copy(new_source)

        print('eigen solve iter {:03d} change {:1.5e} k-eff {:1.6f}'\
                  .format(i, frac_change, k_eff))

        if frac_change < tolerance and i > 2:
            break

    return flux, k_eff

