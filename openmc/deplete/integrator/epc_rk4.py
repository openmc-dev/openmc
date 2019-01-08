""" The EPC-RK4 integrator.

Implements the EPC-RK4 algorithm.
"""

import copy
import os
import time

from mpi4py import MPI

from .cram import CRAM48
from .save_results import save_results

def epc_rk4(operator, print_out=True):
    """ Performs integration of an operator using the EPC-RK4 algorithm.

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
        rates_array.append(copy.deepcopy(rates))

        x_result = []

        t_start = time.time()
        for mat in range(n_mats):
            # Form matrix
            f0 = operator.form_matrix(rates_array[0], mat)

            x_new = CRAM48(1/2 * f0, x[0][mat], dt)

            x_result.append(x_new)

        t_end = time.time()
        if MPI.COMM_WORLD.rank == 0:
            if print_out:
                print("Time to matexp: ", t_end - t_start)

        x.append(x_result)

        eigvl, rates, seed = operator.eval(x[1])

        eigvls.append(eigvl)
        seeds.append(seed)
        rates_array.append(copy.deepcopy(rates))

        x_result = []

        t_start = time.time()
        for mat in range(n_mats):
            # Form matrix
            f1 = operator.form_matrix(rates_array[1], mat)

            x_new = CRAM48(1/2 * f1, x[0][mat], dt)

            x_result.append(x_new)

        t_end = time.time()
        if MPI.COMM_WORLD.rank == 0:
            if print_out:
                print("Time to matexp: ", t_end - t_start)

        x.append(x_result)

        eigvl, rates, seed = operator.eval(x[2])

        eigvls.append(eigvl)
        seeds.append(seed)
        rates_array.append(copy.deepcopy(rates))

        x_result = []

        t_start = time.time()
        for mat in range(n_mats):
            # Form matrix
            f2 = operator.form_matrix(rates_array[2], mat)

            x_new = CRAM48(f2, x[0][mat], dt)

            x_result.append(x_new)

        t_end = time.time()
        if MPI.COMM_WORLD.rank == 0:
            if print_out:
                print("Time to matexp: ", t_end - t_start)

        x.append(x_result)

        eigvl, rates, seed = operator.eval(x[3])

        eigvls.append(eigvl)
        seeds.append(seed)
        rates_array.append(copy.deepcopy(rates))

        x_result = []

        t_start = time.time()
        for mat in range(n_mats):
            # Form matrix
            f0 = operator.form_matrix(rates_array[0], mat)
            f1 = operator.form_matrix(rates_array[1], mat)
            f2 = operator.form_matrix(rates_array[2], mat)
            f3 = operator.form_matrix(rates_array[3], mat)

            x_new = CRAM48(1/6*f0 + 1/3*f1 + 1/3*f2 + 1/6*f3, x[0][mat], dt)

            x_result.append(x_new)

        t_end = time.time()
        if MPI.COMM_WORLD.rank == 0:
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
    eigvl, rates, seed = operator.eval(x[0])

    eigvls.append(eigvl)
    seeds.append(seed)
    rates_array.append(rates)

    # Create results, write to disk
    save_results(operator, x, rates_array, eigvls, seeds, [t, t],
                 len(operator.settings.dt_vec))

    # Return to origin
    os.chdir(dir_home)
