""" The LE/QI CFQ4 integrator.

Implements the LE/QI Predictor-Corrector algorithm using commutator free
high order integrators.

This algorithm is mathematically defined as:

.. math:
    y' = A(y, t) y(t)
    A_m1 = A(y_n-1, t_n-1)
    A_0 = A(y_n, t_n)
    A_l(t) linear extrapolation of A_m1, A_0
    Integrate to t_n+1 to get y_p
    A_c = A(y_p, y_n+1)
    A_q(t) quadratic interpolation of A_m1, A_0, A_c

Here, A(t) is integrated using the fourth order algorithm described below.

From
----
    Thalhammer, Mechthild. "A fourth-order commutator-free exponential
    integrator for nonautonomous differential equations." SIAM journal on
    numerical analysis 44.2 (2006): 851-864.

It is initialized using the CE/LI algorithm.
"""

import copy
import os
import time

from mpi4py import MPI

from .celi_cfq4_imp import celi_cfq4_imp_inner
from .cram import CRAM48
from .save_results import save_results

def leqi_cfq4_imp(operator, print_out=True):
    """ Performs integration of an operator using the LE/QI CFQ4 algorithm.

    Parameters
    ----------
    operator : Operator
        The operator object to simulate on.
    print_out : bool, optional
        Whether or not to print out time.
    """

    m = 10

    # Save current directory
    dir_home = os.getcwd()

    # Move to folder
    os.makedirs(operator.settings.output_dir, exist_ok=True)
    os.chdir(operator.settings.output_dir)

    # Generate initial conditions
    vec = operator.initial_condition()

    n_mats = len(vec)

    t = 0.0

    # Compute initial rates
    operator.settings.particles *= 10
    eigvl_last, rates, seed = operator.eval(vec)
    rates_last = copy.deepcopy(rates)
    operator.settings.particles = int(operator.settings.particles / 10)

    # Perform single step of CE/LI CFQ4 Implicit
    dt_l = operator.settings.dt_vec[0]
    vec, t, rates_bos, eigvl_bos = celi_cfq4_imp_inner(operator, vec, rates_last, eigvl_last, 0, t, dt_l, print_out)

    rates_bar = []

    # Perform remaining LE/QI
    for i, dt in enumerate(operator.settings.dt_vec[1::]):
        # Create vectors
        x = [copy.deepcopy(vec)]

        seeds = [0]
        eigvls = [eigvl_bos]
        rates_array = [copy.deepcopy(rates_bos)]

        # Perform extrapolation
        x_result = []

        t_start = time.time()
        for mat in range(n_mats):
            # Form matrices
            f1 = dt * operator.form_matrix(rates_last, mat)
            f2 = dt * operator.form_matrix(rates_bos, mat)

            # Perform commutator-free integral
            x_new = copy.deepcopy(x[0][mat])

            # Compute linearly extrapolated f at points
            # A{1,2} = f(1/2 -/+ sqrt(3)/6)
            # Then
            # a{1,2} = 1/4 +/- sqrt(3)/6
            # m1 = a2 * A1 + a1 * A2
            # m2 = a1 * A1 + a2 * A2
            m1 = -5 * dt / (12 * dt_l) * f1 + (5 * dt + 6 * dt_l) / (12 * dt_l) * f2
            m2 = -dt / (12 * dt_l) * f1 + (dt + 6 * dt_l) / (12 * dt_l) * f2

            x_new = CRAM48(m2, x_new, 1.0)
            x_new = CRAM48(m1, x_new, 1.0)

            x_result.append(x_new)

        t_end = time.time()
        if MPI.COMM_WORLD.rank == 0:
            if print_out:
                print("Time to matexp: ", t_end - t_start)

        x.append(copy.deepcopy(x_result))
        eigvl_bar = 0.0
        rates_bar = []

        # Loop on inner
        for j in range(0, m + 1):
            eigvl, rates, seed = operator.eval(x_result)

            if j <= 1:
                rates_bar = copy.deepcopy(rates)
                eigvl_bar = eigvl
            else:
                rates_bar.rates = 1/j * rates.rates + (1 - 1/j) * rates_bar.rates
                eigvl_bar = 1/j * eigvl + (1 - 1/j) * eigvl_bar

            x_result = []

            t_start = time.time()
            for mat in range(n_mats):
                # Form matrices
                f1 = dt * operator.form_matrix(rates_last, mat)
                f2 = dt * operator.form_matrix(rates_bos, mat)
                f3 = dt * operator.form_matrix(rates_bar, mat)

                # Perform commutator-free integral
                x_new = copy.deepcopy(x[0][mat])

                # Compute quadratically interpolated f at points
                # A{1,2} = f(1/2 -/+ sqrt(3)/6)
                # Then
                # a{1,2} = 1/4 +/- sqrt(3)/6
                # m1 = a2 * A1 + a1 * A2
                # m2 = a1 * A1 + a2 * A2
                m1 = (-dt**2 / (12 * dt_l * (dt + dt_l)) * f1 +
                     (dt**2 + 2 * dt * dt_l + dt_l**2) / (12 * dt_l * (dt + dt_l)) * f2 +
                     (4 * dt * dt_l + 5 * dt_l**2) / (12 * dt_l * (dt + dt_l)) * f3)
                m2 = (-dt**2/(12 * dt_l * (dt + dt_l)) * f1 +
                     (dt**2 + 6 * dt * dt_l + 5 * dt_l**2) / (12 * dt_l * (dt + dt_l)) * f2 + 
                     dt_l / (12 * (dt + dt_l)) * f3)

                x_new = CRAM48(m2, x_new, 1.0)
                x_new = CRAM48(m1, x_new, 1.0)

                x_result.append(x_new)

            t_end = time.time()
            if MPI.COMM_WORLD.rank == 0:
                if print_out:
                    print("Time to matexp: ", t_end - t_start)

        eigvls.append(eigvl_bar)
        seeds.append(0)
        rates_array.append(copy.deepcopy(rates_bar))

        save_results(operator, x, rates_array, eigvls, seeds, [t, t + dt], i+1)

        rates_last = copy.deepcopy(rates_bos)
        rates_bos = copy.deepcopy(rates_bar)
        eigvl_bos = eigvl_bar
        t += dt
        dt_l = dt
        vec = copy.deepcopy(x_result)

    # Perform one last simulation
    x = [copy.deepcopy(vec)]
    seeds = [0]
    eigvls = [eigvl_bos]
    rates_array = [copy.deepcopy(rates_bos)]

    # Create results, write to disk
    save_results(operator, x, rates_array, eigvls, seeds, [t, t],
                 len(operator.settings.dt_vec))

    # Return to origin
    os.chdir(dir_home)
