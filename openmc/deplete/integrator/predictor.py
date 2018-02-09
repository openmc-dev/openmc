""" The Predictor algorithm."""

import copy
from itertools import repeat
import os
from multiprocessing import Pool
import time

from .. import comm
from .cram import CRAM48, cram_wrapper
from .save_results import save_results


def predictor(operator, print_out=True):
    """The basic predictor integrator.

    Implements the first order predictor algorithm. This algorithm is
    mathematically defined as:

    .. math::
        y' &= A(y, t) y(t)

        A_p &= A(y_n, t_n)

        y_{n+1} &= \\text{expm}(A_p h) y_n

    Parameters
    ----------
    operator : Operator
        The operator object to simulate on.
    print_out : bool, optional
        Whether or not to print out time.
    """

    # Save current directory
    dir_home = os.getcwd()

    # Move to folder
    os.makedirs(operator.settings.output_dir, exist_ok=True)
    os.chdir(operator.settings.output_dir)

    # Generate initial conditions
    vec = operator.initial_condition()

    n_mats = len(vec)

    t = 0.0

    for i, dt in enumerate(operator.settings.dt_vec):
        # Create vectors
        x = [copy.deepcopy(vec)]
        seeds = []
        eigvls = []
        rates_array = []

        eigvl, rates, seed = operator.eval(x[0])

        eigvls.append(eigvl)
        seeds.append(seed)
        rates_array.append(rates)

        # Create results, write to disk
        save_results(operator, x, rates_array, eigvls, seeds, [t, t + dt], i)

        t_start = time.time()

        chains = repeat(operator.chain, n_mats)
        vecs = (x[0][i] for i in range(n_mats))
        rates = (rates_array[0][i, :, :] for i in range(n_mats))
        dts = repeat(dt, n_mats)

        with Pool() as pool:
            iters = zip(chains, vecs, rates, dts)
            x_result = list(pool.starmap(cram_wrapper, iters))

        t_end = time.time()
        if comm.rank == 0:
            if print_out:
                print("Time to matexp: ", t_end - t_start)

        t += dt
        vec = copy.deepcopy(x_result)

    # Perform one last simulation
    x = [copy.deepcopy(vec)]
    seeds = []
    eigvls = []
    rates_array = []
    eigvl, rates, seed = operator.eval(x[0])

    eigvls.append(eigvl)
    seeds.append(seed)
    rates_array.append(rates)

    # Create results, write to disk
    save_results(operator, x, rates_array, eigvls, seeds, [t, t],
                 len(operator.settings.dt_vec))

    # Return to origin
    os.chdir(dir_home)
