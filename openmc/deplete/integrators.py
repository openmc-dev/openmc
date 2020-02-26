import copy
from itertools import repeat

from .abc import Integrator, SIIntegrator, OperatorResult
from .cram import timed_deplete
from ._matrix_funcs import (
    cf4_f1, cf4_f2, cf4_f3, cf4_f4, celi_f1, celi_f2,
    leqi_f1, leqi_f2, leqi_f3, leqi_f4, rk4_f1, rk4_f4
)

__all__ = [
    "PredictorIntegrator", "CECMIntegrator", "CF4Integrator",
    "CELIIntegrator", "EPCRK4Integrator", "LEQIIntegrator",
    "SICELIIntegrator", "SILEQIIntegrator"]


class PredictorIntegrator(Integrator):
    r"""Deplete using a first-order predictor algorithm.

    Implements the first-order predictor algorithm. This algorithm is
    mathematically defined as:

    .. math::
        \begin{aligned}
        y' &= A(y, t) y(t) \\
        A_p &= A(y_n, t_n) \\
        y_{n+1} &= \text{expm}(A_p h) y_n
        \end{aligned}

    Parameters
    ----------
    operator : openmc.deplete.TransportOperator
        Operator to perform transport simulations
    timesteps : iterable of float or iterable of tuple
        Array of timesteps. Note that values are not cumulative. The units are
        specified by the `timestep_units` argument when `timesteps` is an
        iterable of float. Alternatively, units can be specified for each step
        by passing an iterable of (value, unit) tuples.
    power : float or iterable of float, optional
        Power of the reactor in [W]. A single value indicates that
        the power is constant over all timesteps. An iterable
        indicates potentially different power levels for each timestep.
        For a 2D problem, the power can be given in [W/cm] as long
        as the "volume" assigned to a depletion material is actually
        an area in [cm^2]. Either ``power`` or ``power_density`` must be
        specified.
    power_density : float or iterable of float, optional
        Power density of the reactor in [W/gHM]. It is multiplied by
        initial heavy metal inventory to get total power if ``power``
        is not speficied.
    timestep_units : {'s', 'min', 'h', 'd', 'MWd/kg'}
        Units for values specified in the `timesteps` argument. 's' means
        seconds, 'min' means minutes, 'h' means hours, and 'MWd/kg' indicates
        that the values are given in burnup (MW-d of energy deposited per
        kilogram of initial heavy metal).

    Attributes
    ----------
    operator : openmc.deplete.TransportOperator
        Operator to perform transport simulations
    chain : openmc.deplete.Chain
        Depletion chain
    timesteps : iterable of float
        Size of each depletion interval in [s]
    power : iterable of float
        Power of the reactor in [W] for each interval in :attr:`timesteps`
    """
    _num_stages = 1

    def __call__(self, conc, rates, dt, power, _i=None):
        """Perform the integration across one time step

        Parameters
        ----------
        conc : numpy.ndarray
            Initial concentrations for all nuclides in [atom]
        rates : openmc.deplete.ReactionRates
            Reaction rates from operator
        dt : float
            Time in [s] for the entire depletion interval
        power : float
            Power of the system in [W]
        _i : int or None
            Iteration index. Not used

        Returns
        -------
        proc_time : float
            Time spent in CRAM routines for all materials in [s]
        conc_list : list of numpy.ndarray
            Concentrations at end of interval
        op_results : empty list
            Kept for consistency with API. No intermediate calls to
            operator with predictor

        """
        proc_time, conc_end = timed_deplete(self.chain, conc, rates, dt)
        return proc_time, [conc_end], []


class CECMIntegrator(Integrator):
    r"""Deplete using the CE/CM algorithm.

    Implements the second order `CE/CM predictor-corrector algorithm
    <https://doi.org/10.13182/NSE14-92>`_.

    "CE/CM" stands for constant extrapolation on predictor and constant
    midpoint on corrector. This algorithm is mathematically defined as:

    .. math::
        \begin{aligned}
        y' &= A(y, t) y(t) \\
        A_p &= A(y_n, t_n) \\
        y_m &= \text{expm}(A_p h/2) y_n \\
        A_c &= A(y_m, t_n + h/2) \\
        y_{n+1} &= \text{expm}(A_c h) y_n
        \end{aligned}

    Parameters
    ----------
    operator : openmc.deplete.TransportOperator
        Operator to perform transport simulations
    timesteps : iterable of float or iterable of tuple
        Array of timesteps. Note that values are not cumulative. The units are
        specified by the `timestep_units` argument when `timesteps` is an
        iterable of float. Alternatively, units can be specified for each step
        by passing an iterable of (value, unit) tuples.
    power : float or iterable of float, optional
        Power of the reactor in [W]. A single value indicates that
        the power is constant over all timesteps. An iterable
        indicates potentially different power levels for each timestep.
        For a 2D problem, the power can be given in [W/cm] as long
        as the "volume" assigned to a depletion material is actually
        an area in [cm^2]. Either ``power`` or ``power_density`` must be
        specified.
    power_density : float or iterable of float, optional
        Power density of the reactor in [W/gHM]. It is multiplied by
        initial heavy metal inventory to get total power if ``power``
        is not speficied.
    timestep_units : {'s', 'min', 'h', 'd', 'MWd/kg'}
        Units for values specified in the `timesteps` argument. 's' means
        seconds, 'min' means minutes, 'h' means hours, and 'MWd/kg' indicates
        that the values are given in burnup (MW-d of energy deposited per
        kilogram of initial heavy metal).

    Attributes
    ----------
    operator : openmc.deplete.TransportOperator
        Operator to perform transport simulations
    chain : openmc.deplete.Chain
        Depletion chain
    timesteps : iterable of float
        Size of each depletion interval in [s]
    power : iterable of float
        Power of the reactor in [W] for each interval in :attr:`timesteps`
    """
    _num_stages = 2

    def __call__(self, conc, rates, dt, power, _i=None):
        """Integrate using CE/CM

        Parameters
        ----------
        conc : numpy.ndarray
            Initial concentrations for all nuclides in [atom]
        rates : openmc.deplete.ReactionRates
            Reaction rates from operator
        dt : float
            Time in [s] for the entire depletion interval
        power : float
            Power of the system [W]
        _i : int, optional
            Current iteration count. Not used

        Returns
        -------
        proc_time : float
            Time spent in CRAM routines for all materials in [s]
        conc_list : list of numpy.ndarray
            Concentrations at each of the intermediate points with
            the final concentration as the last element
        op_results : list of openmc.deplete.OperatorResult
            Eigenvalue and reaction rates from transport simulations
        """
        # deplete across first half of inteval
        time0, x_middle = timed_deplete(self.chain, conc, rates, dt / 2)
        res_middle = self.operator(x_middle, power)

        # deplete across entire interval with BOS concentrations,
        # MOS reaction rates
        time1, x_end = timed_deplete(self.chain, conc, res_middle.rates, dt)

        return time0 + time1, [x_middle, x_end], [res_middle]


class CF4Integrator(Integrator):
    r"""Deplete using the CF4 algorithm.

    Implements the fourth order `commutator-free Lie algorithm
    <https://doi.org/10.1016/S0167-739X(02)00161-9>`_.
    This algorithm is mathematically defined as:

    .. math::
        \begin{aligned}
        F_1 &= h A(y_0) \\
        y_1 &= \text{expm}(1/2 F_1) y_0 \\
        F_2 &= h A(y_1) \\
        y_2 &= \text{expm}(1/2 F_2) y_0 \\
        F_3 &= h A(y_2) \\
        y_3 &= \text{expm}(-1/2 F_1 + F_3) y_1 \\
        F_4 &= h A(y_3) \\
        y_4 &= \text{expm}( 1/4  F_1 + 1/6 F_2 + 1/6 F_3 - 1/12 F_4)
               \text{expm}(-1/12 F_1 + 1/6 F_2 + 1/6 F_3 + 1/4  F_4) y_0
        \end{aligned}

    Parameters
    ----------
    operator : openmc.deplete.TransportOperator
        Operator to perform transport simulations
    timesteps : iterable of float or iterable of tuple
        Array of timesteps. Note that values are not cumulative. The units are
        specified by the `timestep_units` argument when `timesteps` is an
        iterable of float. Alternatively, units can be specified for each step
        by passing an iterable of (value, unit) tuples.
    power : float or iterable of float, optional
        Power of the reactor in [W]. A single value indicates that
        the power is constant over all timesteps. An iterable
        indicates potentially different power levels for each timestep.
        For a 2D problem, the power can be given in [W/cm] as long
        as the "volume" assigned to a depletion material is actually
        an area in [cm^2]. Either ``power`` or ``power_density`` must be
        specified.
    power_density : float or iterable of float, optional
        Power density of the reactor in [W/gHM]. It is multiplied by
        initial heavy metal inventory to get total power if ``power``
        is not speficied.
    timestep_units : {'s', 'min', 'h', 'd', 'MWd/kg'}
        Units for values specified in the `timesteps` argument. 's' means
        seconds, 'min' means minutes, 'h' means hours, and 'MWd/kg' indicates
        that the values are given in burnup (MW-d of energy deposited per
        kilogram of initial heavy metal).

    Attributes
    ----------
    operator : openmc.deplete.TransportOperator
        Operator to perform transport simulations
    chain : openmc.deplete.Chain
        Depletion chain
    timesteps : iterable of float
        Size of each depletion interval in [s]
    power : iterable of float
        Power of the reactor in [W] for each interval in :attr:`timesteps`
    """
    _num_stages = 4

    def __call__(self, bos_conc, bos_rates, dt, power, _i=None):
        """Perform the integration across one time step

        Parameters
        ----------
        bos_conc : numpy.ndarray
            Initial concentrations for all nuclides in [atom]
        bos_rates : openmc.deplete.ReactionRates
            Reaction rates from operator
        dt : float
            Time in [s] for the entire depletion interval
        power : float
            Power of the system in [W]
        _i : int, optional
            Current depletion step index. Not used

        Returns
        -------
        proc_time : float
            Time spent in CRAM routines for all materials in [s]
        conc_list : list of numpy.ndarray
            Concentrations at each of the intermediate points with
            the final concentration as the last element
        op_results : list of openmc.deplete.OperatorResult
            Eigenvalue and reaction rates from intermediate transport
            simulations
        """
        # Step 1: deplete with matrix 1/2*A(y0)
        time1, conc_eos1 = timed_deplete(
            self.chain, bos_conc, bos_rates, dt, matrix_func=cf4_f1)
        res1 = self.operator(conc_eos1, power)

        # Step 2: deplete with matrix 1/2*A(y1)
        time2, conc_eos2 = timed_deplete(
            self.chain, bos_conc, res1.rates, dt, matrix_func=cf4_f1)
        res2 = self.operator(conc_eos2, power)

        # Step 3: deplete with matrix -1/2*A(y0)+A(y2)
        list_rates = list(zip(bos_rates, res2.rates))
        time3, conc_eos3 = timed_deplete(
            self.chain, conc_eos1, list_rates, dt, matrix_func=cf4_f2)
        res3 = self.operator(conc_eos3, power)

        # Step 4: deplete with two matrix exponentials
        list_rates = list(zip(bos_rates, res1.rates, res2.rates, res3.rates))
        time4, conc_inter = timed_deplete(
            self.chain, bos_conc, list_rates, dt, matrix_func=cf4_f3)
        time5, conc_eos5 = timed_deplete(
            self.chain, conc_inter, list_rates, dt, matrix_func=cf4_f4)

        return (time1 + time2 + time3 + time4 + time5,
                [conc_eos1, conc_eos2, conc_eos3, conc_eos5],
                [res1, res2, res3])


class CELIIntegrator(Integrator):
    r"""Deplete using the CE/LI CFQ4 algorithm.

    Implements the CE/LI Predictor-Corrector algorithm using the `fourth order
    commutator-free integrator <https://doi.org/10.1137/05063042>`_.

    "CE/LI" stands for constant extrapolation on predictor and linear
    interpolation on corrector. This algorithm is mathematically defined as:

    .. math::
        \begin{aligned}
        y' &= A(y, t) y(t) \\
        A_0 &= A(y_n, t_n) \\
        y_p &= \text{expm}(h A_0) y_n \\
        A_1 &= A(y_p, t_n + h) \\
        y_{n+1} &= \text{expm}(\frac{h}{12} A_0 + \frac{5h}{12} A1)
                   \text{expm}(\frac{5h}{12} A_0 + \frac{h}{12} A1) y_n
        \end{aligned}

    Parameters
    ----------
    operator : openmc.deplete.TransportOperator
        Operator to perform transport simulations
    timesteps : iterable of float or iterable of tuple
        Array of timesteps. Note that values are not cumulative. The units are
        specified by the `timestep_units` argument when `timesteps` is an
        iterable of float. Alternatively, units can be specified for each step
        by passing an iterable of (value, unit) tuples.
    power : float or iterable of float, optional
        Power of the reactor in [W]. A single value indicates that
        the power is constant over all timesteps. An iterable
        indicates potentially different power levels for each timestep.
        For a 2D problem, the power can be given in [W/cm] as long
        as the "volume" assigned to a depletion material is actually
        an area in [cm^2]. Either ``power`` or ``power_density`` must be
        specified.
    power_density : float or iterable of float, optional
        Power density of the reactor in [W/gHM]. It is multiplied by
        initial heavy metal inventory to get total power if ``power``
        is not speficied.
    timestep_units : {'s', 'min', 'h', 'd', 'MWd/kg'}
        Units for values specified in the `timesteps` argument. 's' means
        seconds, 'min' means minutes, 'h' means hours, and 'MWd/kg' indicates
        that the values are given in burnup (MW-d of energy deposited per
        kilogram of initial heavy metal).

    Attributes
    ----------
    operator : openmc.deplete.TransportOperator
        Operator to perform transport simulations
    chain : openmc.deplete.Chain
        Depletion chain
    timesteps : iterable of float
        Size of each depletion interval in [s]
    power : iterable of float
        Power of the reactor in [W] for each interval in :attr:`timesteps`
    """
    _num_stages = 2

    def __call__(self, bos_conc, rates, dt, power, _i=None):
        """Perform the integration across one time step

        Parameters
        ----------
        bos_conc : numpy.ndarray
            Initial concentrations for all nuclides in [atom]
        rates : openmc.deplete.ReactionRates
            Reaction rates from operator
        dt : float
            Time in [s] for the entire depletion interval
        power : float
            Power of the system in [W]
        _i : int, optional
            Current iteration count. Not used

        Returns
        -------
        proc_time : float
            Time spent in CRAM routines for all materials in [s]
        conc_list : list of numpy.ndarray
            Concentrations at each of the intermediate points with
            the final concentration as the last element
        op_results : list of openmc.deplete.OperatorResult
            Eigenvalue and reaction rates from intermediate transport
            simulation
        """
        # deplete to end using BOS rates
        proc_time, conc_ce = timed_deplete(self.chain, bos_conc, rates, dt)
        res_ce = self.operator(conc_ce, power)

        # deplete using two matrix exponentials
        list_rates = list(zip(rates, res_ce.rates))

        time_le1, conc_inter = timed_deplete(
            self.chain, bos_conc, list_rates, dt, matrix_func=celi_f1)

        time_le2, conc_end = timed_deplete(
            self.chain, conc_inter, list_rates, dt, matrix_func=celi_f2)

        return proc_time + time_le1 + time_le1, [conc_ce, conc_end], [res_ce]


class EPCRK4Integrator(Integrator):
    r"""Deplete using the EPC-RK4 algorithm.

    Implements an extended predictor-corrector algorithm with traditional
    Runge-Kutta 4 method. This algorithm is mathematically defined as:

    .. math::
        \begin{aligned}
        F_1 &= h A(y_0) \\
        y_1 &= \text{expm}(1/2 F_1) y_0 \\
        F_2 &= h A(y_1) \\
        y_2 &= \text{expm}(1/2 F_2) y_0 \\
        F_3 &= h A(y_2) \\
        y_3 &= \text{expm}(F_3) y_0 \\
        F_4 &= h A(y_3) \\
        y_4 &= \text{expm}(1/6 F_1 + 1/3 F_2 + 1/3 F_3 + 1/6 F_4) y_0
        \end{aligned}

    Parameters
    ----------
    operator : openmc.deplete.TransportOperator
        Operator to perform transport simulations
    timesteps : iterable of float or iterable of tuple
        Array of timesteps. Note that values are not cumulative. The units are
        specified by the `timestep_units` argument when `timesteps` is an
        iterable of float. Alternatively, units can be specified for each step
        by passing an iterable of (value, unit) tuples.
    power : float or iterable of float, optional
        Power of the reactor in [W]. A single value indicates that
        the power is constant over all timesteps. An iterable
        indicates potentially different power levels for each timestep.
        For a 2D problem, the power can be given in [W/cm] as long
        as the "volume" assigned to a depletion material is actually
        an area in [cm^2]. Either ``power`` or ``power_density`` must be
        specified.
    power_density : float or iterable of float, optional
        Power density of the reactor in [W/gHM]. It is multiplied by
        initial heavy metal inventory to get total power if ``power``
        is not speficied.
    timestep_units : {'s', 'min', 'h', 'd', 'MWd/kg'}
        Units for values specified in the `timesteps` argument. 's' means
        seconds, 'min' means minutes, 'h' means hours, and 'MWd/kg' indicates
        that the values are given in burnup (MW-d of energy deposited per
        kilogram of initial heavy metal).

    Attributes
    ----------
    operator : openmc.deplete.TransportOperator
        Operator to perform transport simulations
    chain : openmc.deplete.Chain
        Depletion chain
    timesteps : iterable of float
        Size of each depletion interval in [s]
    power : iterable of float
        Power of the reactor in [W] for each interval in :attr:`timesteps`
    """
    _num_stages = 4

    def __call__(self, conc, rates, dt, power, _i=None):
        """Perform the integration across one time step

        Parameters
        ----------
        conc : numpy.ndarray
            Initial concentrations for all nuclides in [atom]
        rates : openmc.deplete.ReactionRates
            Reaction rates from operator
        dt : float
            Time in [s] for the entire depletion interval
        power : float
            Power of the system in [W]
        _i : int, optional
            Current depletion step index, unused.

        Returns
        -------
        proc_time : float
            Time spent in CRAM routines for all materials in [s]
        conc_list : list of numpy.ndarray
            Concentrations at each of the intermediate points with
            the final concentration as the last element
        op_results : list of openmc.deplete.OperatorResult
            Eigenvalue and reaction rates from intermediate transport
            simulations
        """

        # Step 1: deplete with matrix A(y0) / 2
        time1, conc1 = timed_deplete(
            self.chain, conc, rates, dt, matrix_func=rk4_f1)
        res1 = self.operator(conc1, power)

        # Step 2: deplete with matrix A(y1) / 2
        time2, conc2 = timed_deplete(
            self.chain, conc, res1.rates, dt, matrix_func=rk4_f1)
        res2 = self.operator(conc2, power)

        # Step 3: deplete with matrix A(y2)
        time3, conc3 = timed_deplete(
            self.chain, conc, res2.rates, dt)
        res3 = self.operator(conc3, power)

        # Step 4: deplete with matrix built from weighted rates
        list_rates = list(zip(rates, res1.rates, res2.rates, res3.rates))
        time4, conc4 = timed_deplete(
            self.chain, conc, list_rates, dt, matrix_func=rk4_f4)

        return (time1 + time2 + time3 + time4, [conc1, conc2, conc3, conc4],
                [res1, res2, res3])


class LEQIIntegrator(Integrator):
    r"""Deplete using the LE/QI CFQ4 algorithm.

    Implements the LE/QI Predictor-Corrector algorithm using the `fourth order
    commutator-free integrator <https://doi.org/10.1137/05063042>`_.

    "LE/QI" stands for linear extrapolation on predictor and quadratic
    interpolation on corrector. This algorithm is mathematically defined as:

    .. math::
        \begin{aligned}
        y' &= A(y, t) y(t) \\
        A_{last} &= A(y_{n-1}, t_n - h_1) \\
        A_0 &= A(y_n, t_n) \\
        F_1 &= \frac{-h_2^2}{12h_1} A_{last} + \frac{h_2(6h_1+h_2)}{12h_1} A_0 \\
        F_2 &= \frac{-5h_2^2}{12h_1} A_{last} + \frac{h_2(6h_1+5h_2)}{12h_1} A_0 \\
        y_p &= \text{expm}(F_2) \text{expm}(F_1) y_n \\
        A_1 &= A(y_p, t_n + h_2) \\
        F_3 &= \frac{-h_2^3}{12 h_1 (h_1 + h_2)} A_{last} +
              \frac{h_2 (5 h_1^2 + 6 h_2 h_1 + h_2^2)}{12 h_1 (h_1 + h_2)} A_0 +
              \frac{h_2 h_1)}{12 (h_1 + h_2)} A_1 \\
        F_4 &= \frac{-h_2^3}{12 h_1 (h_1 + h_2)} A_{last} +
              \frac{h_2 (h_1^2 + 2 h_2 h_1 + h_2^2)}{12 h_1 (h_1 + h_2)} A_0 +
              \frac{h_2 (5 h_1^2 + 4 h_2 h_1)}{12 h_1 (h_1 + h_2)} A_1 \\
        y_{n+1} &= \text{expm}(F_4) \text{expm}(F_3) y_n
        \end{aligned}

    It is initialized using the CE/LI algorithm.

    Parameters
    ----------
    operator : openmc.deplete.TransportOperator
        Operator to perform transport simulations
    timesteps : iterable of float or iterable of tuple
        Array of timesteps. Note that values are not cumulative. The units are
        specified by the `timestep_units` argument when `timesteps` is an
        iterable of float. Alternatively, units can be specified for each step
        by passing an iterable of (value, unit) tuples.
    power : float or iterable of float, optional
        Power of the reactor in [W]. A single value indicates that
        the power is constant over all timesteps. An iterable
        indicates potentially different power levels for each timestep.
        For a 2D problem, the power can be given in [W/cm] as long
        as the "volume" assigned to a depletion material is actually
        an area in [cm^2]. Either ``power`` or ``power_density`` must be
        specified.
    power_density : float or iterable of float, optional
        Power density of the reactor in [W/gHM]. It is multiplied by
        initial heavy metal inventory to get total power if ``power``
        is not speficied.
    timestep_units : {'s', 'min', 'h', 'd', 'MWd/kg'}
        Units for values specified in the `timesteps` argument. 's' means
        seconds, 'min' means minutes, 'h' means hours, and 'MWd/kg' indicates
        that the values are given in burnup (MW-d of energy deposited per
        kilogram of initial heavy metal).

    Attributes
    ----------
    operator : openmc.deplete.TransportOperator
        Operator to perform transport simulations
    chain : openmc.deplete.Chain
        Depletion chain
    timesteps : iterable of float
        Size of each depletion interval in [s]
    power : iterable of float
        Power of the reactor in [W] for each interval in :attr:`timesteps`
    """
    _num_stages = 2

    def __call__(self, bos_conc, bos_rates, dt, power, i):
        """Perform the integration across one time step

        Parameters
        ----------
        conc : numpy.ndarray
            Initial concentrations for all nuclides in [atom]
        rates : openmc.deplete.ReactionRates
            Reaction rates from operator
        dt : float
            Time in [s] for the entire depletion interval
        power : float
            Power of the system in [W]
        i : int
            Current depletion step index

        Returns
        -------
        proc_time : float
            Time spent in CRAM routines for all materials in [s]
        conc_list : list of numpy.ndarray
            Concentrations at each of the intermediate points with
            the final concentration as the last element
        op_results : list of openmc.deplete.OperatorResult
            Eigenvalue and reaction rates from intermediate transport
            simulation
        """
        if i == 0:
            if self._i_res < 1:  # need at least previous transport solution
                self._prev_rates = bos_rates
                return CELIIntegrator.__call__(
                    self, bos_conc, bos_rates, dt, power, i)
            prev_res = self.operator.prev_res[-2]
            prev_dt = self.timesteps[i] - prev_res.time[0]
            self._prev_rates = prev_res.rates[0]
        else:
            prev_dt = self.timesteps[i - 1]

        # Remaining LE/QI
        bos_res = self.operator(bos_conc, power)

        le_inputs = list(zip(
            self._prev_rates, bos_res.rates, repeat(prev_dt), repeat(dt)))

        time1, conc_inter = timed_deplete(
            self.chain, bos_conc, le_inputs, dt, matrix_func=leqi_f1)
        time2, conc_eos0 = timed_deplete(
            self.chain, conc_inter, le_inputs, dt, matrix_func=leqi_f2)

        res_inter = self.operator(conc_eos0, power)

        qi_inputs = list(zip(
            self._prev_rates, bos_res.rates, res_inter.rates,
            repeat(prev_dt), repeat(dt)))

        time3, conc_inter = timed_deplete(
            self.chain, bos_conc, qi_inputs, dt, matrix_func=leqi_f3)
        time4, conc_eos1 = timed_deplete(
            self.chain, conc_inter, qi_inputs, dt, matrix_func=leqi_f4)

        # store updated rates
        self._prev_rates = copy.deepcopy(bos_res.rates)

        return (
            time1 + time2 + time3 + time4, [conc_eos0, conc_eos1],
            [bos_res, res_inter])


class SICELIIntegrator(SIIntegrator):
    r"""Deplete using the SI-CE/LI CFQ4 algorithm.

    Implements the stochastic implicit CE/LI predictor-corrector algorithm
    using the `fourth order commutator-free integrator
    <https://doi.org/10.1137/05063042>`_.

    Detailed algorithm can be found in section 3.2 in `Colin Josey's thesis
    <http://hdl.handle.net/1721.1/113721>`_.

    Parameters
    ----------
    operator : openmc.deplete.TransportOperator
        The operator object to simulate on.
    timesteps : iterable of float or iterable of tuple
        Array of timesteps. Note that values are not cumulative. The units are
        specified by the `timestep_units` argument when `timesteps` is an
        iterable of float. Alternatively, units can be specified for each step
        by passing an iterable of (value, unit) tuples.
    power : float or iterable of float, optional
        Power of the reactor in [W]. A single value indicates that
        the power is constant over all timesteps. An iterable
        indicates potentially different power levels for each timestep.
        For a 2D problem, the power can be given in [W/cm] as long
        as the "volume" assigned to a depletion material is actually
        an area in [cm^2]. Either ``power`` or ``power_density`` must be
        specified.
    power_density : float or iterable of float, optional
        Power density of the reactor in [W/gHM]. It is multiplied by
        initial heavy metal inventory to get total power if ``power``
        is not speficied.
    timestep_units : {'s', 'min', 'h', 'd', 'MWd/kg'}
        Units for values specified in the `timesteps` argument. 's' means
        seconds, 'min' means minutes, 'h' means hours, and 'MWd/kg' indicates
        that the values are given in burnup (MW-d of energy deposited per
        kilogram of initial heavy metal).
    n_steps : int, optional
        Number of stochastic iterations per depletion interval.
        Must be greater than zero. Default : 10

    Attributes
    ----------
    operator : openmc.deplete.TransportOperator
        Operator to perform transport simulations
    chain : openmc.deplete.Chain
        Depletion chain
    timesteps : iterable of float
        Size of each depletion interval in [s]
    power : iterable of float
        Power of the reactor in [W] for each interval in :attr:`timesteps`
    n_steps : int
        Number of stochastic iterations per depletion interval
    """
    _num_stages = 2

    def __call__(self, bos_conc, bos_rates, dt, power, _i=None):
        """Perform the integration across one time step

        Parameters
        ----------
        bos_conc : numpy.ndarray
            Initial bos_concentrations for all nuclides in [atom]
        bos_rates : openmc.deplete.ReactionRates
            Reaction rates from operator
        dt : float
            Time in [s] for the entire depletion interval
        power : float
            Power of the system in [W]
        _i : int, optional
            Current depletion step index. Not used

        Returns
        -------
        proc_time : float
            Time spent in CRAM routines for all materials in [s]
        bos_conc_list : list of numpy.ndarray
            Concentrations at each of the intermediate points with
            the final bos_concentration as the last element
        op_results : list of openmc.deplete.OperatorResult
            Eigenvalue and reaction rates from intermediate transport
            simulations
        """
        proc_time, eos_conc = timed_deplete(
            self.chain, bos_conc, bos_rates, dt)
        inter_conc = copy.deepcopy(eos_conc)

        # Begin iteration
        for j in range(self.n_steps + 1):
            inter_res = self.operator(inter_conc, power)

            if j <= 1:
                res_bar = copy.deepcopy(inter_res)
            else:
                rates = 1/j * inter_res.rates + (1 - 1 / j) * res_bar.rates
                k = 1/j * inter_res.k + (1 - 1 / j) * res_bar.k
                res_bar = OperatorResult(k, rates)

            list_rates = list(zip(bos_rates, res_bar.rates))
            time1, inter_conc = timed_deplete(
                self.chain, bos_conc, list_rates, dt, matrix_func=celi_f1)
            time2, inter_conc = timed_deplete(
                self.chain, inter_conc, list_rates, dt, matrix_func=celi_f2)
            proc_time += time1 + time2

        # end iteration
        return proc_time, [eos_conc, inter_conc], [res_bar]


class SILEQIIntegrator(SIIntegrator):
    r"""Deplete using the SI-LE/QI CFQ4 algorithm.

    Implements the Stochastic Implicit LE/QI Predictor-Corrector algorithm
    using the `fourth order commutator-free integrator
    <https://doi.org/10.1137/05063042>`_.

    Detailed algorithm can be found in Section 3.2 in `Colin Josey's thesis
    <http://hdl.handle.net/1721.1/113721>`_.

    Parameters
    ----------
    operator : openmc.deplete.TransportOperator
        The operator object to simulate on.
    timesteps : iterable of float or iterable of tuple
        Array of timesteps. Note that values are not cumulative. The units are
        specified by the `timestep_units` argument when `timesteps` is an
        iterable of float. Alternatively, units can be specified for each step
        by passing an iterable of (value, unit) tuples.
    power : float or iterable of float, optional
        Power of the reactor in [W]. A single value indicates that
        the power is constant over all timesteps. An iterable
        indicates potentially different power levels for each timestep.
        For a 2D problem, the power can be given in [W/cm] as long
        as the "volume" assigned to a depletion material is actually
        an area in [cm^2]. Either ``power`` or ``power_density`` must be
        specified.
    power_density : float or iterable of float, optional
        Power density of the reactor in [W/gHM]. It is multiplied by
        initial heavy metal inventory to get total power if ``power``
        is not speficied.
    timestep_units : {'s', 'min', 'h', 'd', 'MWd/kg'}
        Units for values specified in the `timesteps` argument. 's' means
        seconds, 'min' means minutes, 'h' means hours, and 'MWd/kg' indicates
        that the values are given in burnup (MW-d of energy deposited per
        kilogram of initial heavy metal).
    n_steps : int, optional
        Number of stochastic iterations per depletion interval.
        Must be greater than zero. Default : 10

    Attributes
    ----------
    operator : openmc.deplete.TransportOperator
        Operator to perform transport simulations
    chain : openmc.deplete.Chain
        Depletion chain
    timesteps : iterable of float
        Size of each depletion interval in [s]
    power : iterable of float
        Power of the reactor in [W] for each interval in :attr:`timesteps`
    n_steps : int
        Number of stochastic iterations per depletion interval
    """
    _num_stages = 2

    def __call__(self, bos_conc, bos_rates, dt, power, i):
        """Perform the integration across one time step

        Parameters
        ----------
        bos_conc : list of numpy.ndarray
            Initial concentrations for all nuclides in [atom] for
            all depletable materials
        bos_rates : list of openmc.deplete.ReactionRates
            Reaction rates from operator for all depletable materials
        dt : float
            Time in [s] for the entire depletion interval
        power : float
            Power of the system in [W]
        i : int
            Current depletion step index

        Returns
        -------
        proc_time : float
            Time spent in CRAM routines for all materials in [s]
        conc_list : list of numpy.ndarray
            Concentrations at each of the intermediate points with
            the final concentration as the last element
        op_results : list of openmc.deplete.OperatorResult
            Eigenvalue and reaction rates from intermediate transport
            simulation
        """
        if i == 0:
            if self._i_res < 1:
                self._prev_rates = bos_rates
                # Perform CELI for initial steps
                return SICELIIntegrator.__call__(
                    self, bos_conc, bos_rates, dt, power, i)
            prev_res = self.operator.prev_res[-2]
            prev_dt = self.timesteps[i] - prev_res.time[0]
            self._prev_rates = prev_res.rates[0]
        else:
            prev_dt = self.timesteps[i - 1]

        # Perform remaining LE/QI
        inputs = list(zip(self._prev_rates, bos_rates,
                          repeat(prev_dt), repeat(dt)))
        proc_time, inter_conc = timed_deplete(
            self.chain, bos_conc, inputs, dt, matrix_func=leqi_f1)
        time1, eos_conc = timed_deplete(
            self.chain, inter_conc, inputs, dt, matrix_func=leqi_f2)

        proc_time += time1
        inter_conc = copy.deepcopy(eos_conc)

        for j in range(self.n_steps + 1):
            inter_res = self.operator(inter_conc, power)

            if j <= 1:
                res_bar = copy.deepcopy(inter_res)
            else:
                rates = 1 / j * inter_res.rates + (1 - 1 / j) * res_bar.rates
                k = 1 / j * inter_res.k + (1 - 1 / j) * res_bar.k
                res_bar = OperatorResult(k, rates)

            inputs = list(zip(self._prev_rates, bos_rates, res_bar.rates,
                              repeat(prev_dt), repeat(dt)))
            time1, inter_conc = timed_deplete(
                self.chain, bos_conc, inputs, dt, matrix_func=leqi_f3)
            time2, inter_conc = timed_deplete(
                self.chain, inter_conc, inputs, dt, matrix_func=leqi_f4)
            proc_time += time1 + time2

        return proc_time, [eos_conc, inter_conc], [res_bar]
