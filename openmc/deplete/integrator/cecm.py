""" The CE/CM integrator."""

import copy
from itertools import repeat
import os
from multiprocessing import Pool
import time

from .. import comm
from .cram import CRAM48, cram_wrapper
from .save_results import save_results


def cecm(operator, print_out=True):
    """The CE/CM integrator.

    Implements the second order CE/CM Predictor-Corrector algorithm [ref]_.
    This algorithm is mathematically defined as:

    .. math::
        y' &= A(y, t) y(t)

        A_p &= A(y_n, t_n)

        y_m &= \\text{expm}(A_p h/2) y_n

        A_c &= A(y_m, t_n + h/2)

        y_{n+1} &= \\text{expm}(A_c h) y_n

    .. [ref]
        Isotalo, Aarno. "Comparison of Neutronics-Depletion Coupling Schemes
        for Burnup Calculations-Continued Study." Nuclear Science and
        Engineering 180.3 (2015): 286-300.

    Parameters
    ----------
    operator : Operator
        The operator object to simulate on.
    print_out : bool, optional
        Whether or not to print out time.
    """

    # Generate initial conditions
    with operator as vec:
        n_mats = len(vec)

        t = 0.0

        for i, dt in enumerate(operator.settings.dt_vec):
            # Create vectors
            x = [copy.deepcopy(vec)]
            seeds = []
            eigvls = []
            rates_array = []

            eigvl, rates, seed = operator(x[0])

            eigvls.append(eigvl)
            seeds.append(seed)
            rates_array.append(rates)

            t_start = time.time()

            chains = repeat(operator.chain, n_mats)
            vecs = (x[0][i] for i in range(n_mats))
            rates = (rates_array[0][i, :, :] for i in range(n_mats))
            dts = repeat(dt/2, n_mats)

            with Pool() as pool:
                iters = zip(chains, vecs, rates, dts)
                x_result = list(pool.starmap(cram_wrapper, iters))

            t_end = time.time()
            if comm.rank == 0:
                if print_out:
                    print("Time to matexp: ", t_end - t_start)

            x.append(x_result)

            eigvl, rates, seed = operator(x[1])

            eigvls.append(eigvl)
            seeds.append(seed)
            rates_array.append(rates)

            t_start = time.time()

            chains = repeat(operator.chain, n_mats)
            vecs = (x[0][i] for i in range(n_mats))
            rates = (rates_array[1][i, :, :] for i in range(n_mats))
            dts = repeat(dt, n_mats)

            with Pool() as pool:
                iters = zip(chains, vecs, rates, dts)
                x_result = list(pool.starmap(cram_wrapper, iters))

            t_end = time.time()
            if comm.rank == 0:
                if print_out:
                    print("Time to matexp: ", t_end - t_start)

            # Create results, write to disk
            save_results(operator, x, rates_array, eigvls, seeds, [t, t + dt], i)

            t += dt
            vec = copy.deepcopy(x_result)

        # Perform one last simulation
        x = [copy.deepcopy(vec)]
        seeds = []
        eigvls = []
        rates_array = []
        eigvl, rates, seed = operator(x[0])

        eigvls.append(eigvl)
        seeds.append(seed)
        rates_array.append(rates)

        # Create results, write to disk
        save_results(operator, x, rates_array, eigvls, seeds, [t, t],
                     len(operator.settings.dt_vec))
