""" The CE/LI CFQ4 integrator.

Implements the CE/LI Predictor-Corrector algorithm using commutator free
high order integrators.

This algorithm is mathematically defined as:

.. math:
    y' = A(y, t) y(t)
    A_p = A(y_n, t_n)
    y_p = expm(A_p h) y_n
    A_c = A(y_p, t_n)
    A(t) = t/dt * A_c + (dt - t)/dt * A_p

Here, A(t) is integrated using the fourth order algorithm described below.

From
----
    Thalhammer, Mechthild. "A fourth-order commutator-free exponential
    integrator for nonautonomous differential equations." SIAM journal on
    numerical analysis 44.2 (2006): 851-864.
"""

import copy
import os
import time

from mpi4py import MPI

from .cram import CRAM48
from .save_results import save_results

def celi_cfq4(operator, print_out=True):
    """ Performs integration of an operator using the CE/LI CFQ4 algorithm.

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

    t = 0.0

    for i, dt in enumerate(operator.settings.dt_vec):
        vec, t, _ = celi_cfq4_inner(operator, vec, i, t, dt, print_out)

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

def celi_cfq4_inner(operator, vec, i, t, dt, print_out):
    """ The inner loop of CE/LI CFQ4.

    Parameters
    ----------
    operator : Operator
        The operator object to simulate on.
    vec : list of numpy.array
        Nuclide vector, beginning of time.
    i : Int
        Current iteration number.
    t : Float
        Time at start of step.
    dt : Float
        Time step.
    print_out : bool
        Whether or not to print out time.

    Returns
    -------
    x_result : list of numpy.array
        Nuclide vector, end of time.
    Float
        Next time
    ReactionRates
        Reaction rates from beginning of step.
    """

    n_mats = len(vec)

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
        f = operator.form_matrix(rates_array[0], mat)

        x_new = CRAM48(f, x[0][mat], dt)

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
        # Form matrices
        f1 = dt * operator.form_matrix(rates_array[0], mat)
        f2 = dt * operator.form_matrix(rates_array[1], mat)

        # Perform commutator-free integral
        x_new = copy.deepcopy(x[0][mat])

        # Compute linearly interpolated f at points
        # A{1,2} = f(1/2 -/+ sqrt(3)/6)
        # Then
        # a{1,2} = 1/4 +/- sqrt(3)/6
        # m1 = a2 * A1 + a1 * A2
        # m2 = a1 * A1 + a2 * A2
        m1 = 1/12 * (f1 + 5 * f2)
        m2 = 1/12 * (5 * f1 + f2)

        x_new = CRAM48(m2, x_new, 1.0)
        x_new = CRAM48(m1, x_new, 1.0)

        x_result.append(x_new)

    t_end = time.time()
    if MPI.COMM_WORLD.rank == 0:
        if print_out:
            print("Time to matexp: ", t_end - t_start)

    # Create results, write to disk
    save_results(operator, x, rates_array, eigvls, seeds, [t, t + dt], i)

    return x_result, t + dt, rates_array[0]