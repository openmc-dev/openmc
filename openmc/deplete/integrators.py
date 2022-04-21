import copy
from itertools import repeat

from .abc import Integrator, SIIntegrator, OperatorResult, add_params
from ._matrix_funcs import (
    cf4_f1, cf4_f2, cf4_f3, cf4_f4, celi_f1, celi_f2,
    leqi_f1, leqi_f2, leqi_f3, leqi_f4, rk4_f1, rk4_f4
)

__all__ = [
    "PredictorIntegrator", "CECMIntegrator", "CF4Integrator",
    "CELIIntegrator", "EPCRK4Integrator", "LEQIIntegrator",
    "SICELIIntegrator", "SILEQIIntegrator"]


@add_params
class PredictorIntegrator(Integrator):
    r"""Deplete using a first-order predictor algorithm.

    Implements the first-order predictor algorithm. This algorithm is
    mathematically defined as:

    .. math::
        \mathbf{n}_{i+1} = \exp\left(h\mathbf{A}(\mathbf{n}_i) \right) \mathbf{n}_i

    """
    _num_stages = 1

    def __call__(self, conc, rates, dt, source_rate, _i=None):
        """Perform the integration across one time step

        Parameters
        ----------
        conc : numpy.ndarray
            Initial concentrations for all nuclides in [atom]
        rates : openmc.deplete.ReactionRates
            Reaction rates from operator
        dt : float
            Time in [s] for the entire depletion interval
        source_rate : float
            Power in [W] or source rate in [neutron/sec]
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
        proc_time, conc_end = self._timed_deplete(conc, rates, dt)
        return proc_time, [conc_end], []


@add_params
class CECMIntegrator(Integrator):
    r"""Deplete using the CE/CM algorithm.

    Implements the second order `CE/CM predictor-corrector algorithm
    <https://doi.org/10.13182/NSE14-92>`_.

    "CE/CM" stands for constant extrapolation on predictor and constant
    midpoint on corrector. This algorithm is mathematically defined as:

    .. math::
        \begin{aligned}
        \mathbf{n}_{i+1/2} &= \exp \left (\frac{h}{2}\mathbf{A}(\mathbf{n}_i)
            \right) \mathbf{n}_i \\
        \mathbf{n}_{i+1} &= \exp \left(h \mathbf{A}(\mathbf{n}_{i+1/2}) \right)
            \mathbf{n}_i.
        \end{aligned}
    """
    _num_stages = 2

    def __call__(self, conc, rates, dt, source_rate, _i=None):
        """Integrate using CE/CM

        Parameters
        ----------
        conc : numpy.ndarray
            Initial concentrations for all nuclides in [atom]
        rates : openmc.deplete.ReactionRates
            Reaction rates from operator
        dt : float
            Time in [s] for the entire depletion interval
        source_rate : float
            Power in [W] or source rate in [neutron/sec]
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
        time0, x_middle = self._timed_deplete(conc, rates, dt / 2)
        res_middle = self.operator(x_middle, source_rate)

        # deplete across entire interval with BOS concentrations,
        # MOS reaction rates
        time1, x_end = self._timed_deplete(conc, res_middle.rates, dt)

        return time0 + time1, [x_middle, x_end], [res_middle]


@add_params
class CF4Integrator(Integrator):
    r"""Deplete using the CF4 algorithm.

    Implements the fourth order `commutator-free Lie algorithm
    <https://doi.org/10.1016/S0167-739X(02)00161-9>`_.
    This algorithm is mathematically defined as:

    .. math::
        \begin{aligned}
        \mathbf{A}_1 &= h\mathbf{A}(\mathbf{n}_i) \\
        \hat{\mathbf{n}}_1 &= \exp \left ( \frac{\mathbf{A}_1}{2} \right ) \mathbf{n}_i \\
        \mathbf{A}_2 &= h\mathbf{A}(\hat{\mathbf{n}}_1) \\
        \hat{\mathbf{n}}_2 &= \exp \left ( \frac{\mathbf{A}_2}{2} \right ) \mathbf{n}_i \\
        \mathbf{A}_3 &= h \mathbf{A}(\hat{\mathbf{n}}_2) \\
        \hat{\mathbf{n}}_3 &= \exp \left ( -\frac{\mathbf{A}_1}{2} + \mathbf{A}_3
            \right ) \hat{\mathbf{n}}_1 \\
        \mathbf{A}_4 &= h\mathbf{A}(\hat{\mathbf{n}}_3) \\
        \mathbf{n}_{i+1} &= \exp \left ( \frac{\mathbf{A}_1}{4} + \frac{\mathbf{A}_2}{6}
            + \frac{\mathbf{A}_3}{6} - \frac{\mathbf{A}_4}{12} \right )
        \exp \left ( -\frac{\mathbf{A}_1}{12} + \frac{\mathbf{A}_2}{6} +
            \frac{\mathbf{A}_3}{6} - \frac{\mathbf{A}_4}{4} \right ) \mathbf{n}_i.
        \end{aligned}
    """
    _num_stages = 4

    def __call__(self, bos_conc, bos_rates, dt, source_rate, _i=None):
        """Perform the integration across one time step

        Parameters
        ----------
        bos_conc : numpy.ndarray
            Initial concentrations for all nuclides in [atom]
        bos_rates : openmc.deplete.ReactionRates
            Reaction rates from operator
        dt : float
            Time in [s] for the entire depletion interval
        source_rate : float
            Power in [W] or source rate in [neutron/sec]
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
        time1, conc_eos1 = self._timed_deplete(
            bos_conc, bos_rates, dt, matrix_func=cf4_f1)
        res1 = self.operator(conc_eos1, source_rate)

        # Step 2: deplete with matrix 1/2*A(y1)
        time2, conc_eos2 = self._timed_deplete(
            bos_conc, res1.rates, dt, matrix_func=cf4_f1)
        res2 = self.operator(conc_eos2, source_rate)

        # Step 3: deplete with matrix -1/2*A(y0)+A(y2)
        list_rates = list(zip(bos_rates, res2.rates))
        time3, conc_eos3 = self._timed_deplete(
            conc_eos1, list_rates, dt, matrix_func=cf4_f2)
        res3 = self.operator(conc_eos3, source_rate)

        # Step 4: deplete with two matrix exponentials
        list_rates = list(zip(bos_rates, res1.rates, res2.rates, res3.rates))
        time4, conc_inter = self._timed_deplete(
            bos_conc, list_rates, dt, matrix_func=cf4_f3)
        time5, conc_eos5 = self._timed_deplete(
            conc_inter, list_rates, dt, matrix_func=cf4_f4)

        return (time1 + time2 + time3 + time4 + time5,
                [conc_eos1, conc_eos2, conc_eos3, conc_eos5],
                [res1, res2, res3])


@add_params
class CELIIntegrator(Integrator):
    r"""Deplete using the CE/LI CFQ4 algorithm.

    Implements the CE/LI Predictor-Corrector algorithm using the `fourth order
    commutator-free integrator <https://doi.org/10.1137/05063042>`_.

    "CE/LI" stands for constant extrapolation on predictor and linear
    interpolation on corrector. This algorithm is mathematically defined as:

    .. math::
        \begin{aligned}
        \mathbf{n}_{i+1}^p &= \exp \left ( h \mathbf{A}(\mathbf{n}_i ) \right )
        \mathbf{n}_i \\
        \mathbf{n}_{i+1} &= \exp \left( \frac{h}{12} \mathbf{A}(\mathbf{n}_i) +
            \frac{5h}{12} \mathbf{A}(\mathbf{n}_{i+1}^p) \right)
        \exp \left( \frac{5h}{12} \mathbf{A}(\mathbf{n}_i) +
        \frac{h}{12} \mathbf{A}(\mathbf{n}_{i+1}^p) \right) \mathbf{n}_i.
        \end{aligned}
    """
    _num_stages = 2

    def __call__(self, bos_conc, rates, dt, source_rate, _i=None):
        """Perform the integration across one time step

        Parameters
        ----------
        bos_conc : numpy.ndarray
            Initial concentrations for all nuclides in [atom]
        rates : openmc.deplete.ReactionRates
            Reaction rates from operator
        dt : float
            Time in [s] for the entire depletion interval
        source_rate : float
            Power in [W] or source rate in [neutron/sec]
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
        proc_time, conc_ce = self._timed_deplete(bos_conc, rates, dt)
        res_ce = self.operator(conc_ce, source_rate)

        # deplete using two matrix exponentials
        list_rates = list(zip(rates, res_ce.rates))

        time_le1, conc_inter = self._timed_deplete(
            bos_conc, list_rates, dt, matrix_func=celi_f1)

        time_le2, conc_end = self._timed_deplete(
            conc_inter, list_rates, dt, matrix_func=celi_f2)

        return proc_time + time_le1 + time_le1, [conc_ce, conc_end], [res_ce]


@add_params
class EPCRK4Integrator(Integrator):
    r"""Deplete using the EPC-RK4 algorithm.

    Implements an extended predictor-corrector algorithm with traditional
    Runge-Kutta 4 method. This algorithm is mathematically defined as:

    .. math::
        \begin{aligned}
        \mathbf{A}_1 &= h\mathbf{A}(\mathbf{n}_i) \\
        \hat{\mathbf{n}}_1 &= \exp \left ( \frac{\mathbf{A}_1}{2} \right ) \mathbf{n}_i \\
        \mathbf{A}_2 &= h\mathbf{A}(\hat{\mathbf{n}}_1) \\
        \hat{\mathbf{n}}_2 &= \exp \left ( \frac{\mathbf{A}_2}{2} \right ) \mathbf{n}_i \\
        \mathbf{A}_3 &= h \mathbf{A}(\hat{\mathbf{n}}_2) \\
        \hat{\mathbf{n}}_3 &= \exp \left ( \mathbf{A}_3 \right ) \mathbf{n}_i \\
        \mathbf{A}_4 &= h\mathbf{A}(\hat{\mathbf{n}}_3) \\
        \mathbf{n}_{i+1} &= \exp \left ( \frac{\mathbf{A}_1}{6} + \frac{\mathbf{A}_2}{3}
        + \frac{\mathbf{A}_3}{3} + \frac{\mathbf{A}_4}{6} \right ) \mathbf{n}_i.
        \end{aligned}
    """
    _num_stages = 4

    def __call__(self, conc, rates, dt, source_rate, _i=None):
        """Perform the integration across one time step

        Parameters
        ----------
        conc : numpy.ndarray
            Initial concentrations for all nuclides in [atom]
        rates : openmc.deplete.ReactionRates
            Reaction rates from operator
        dt : float
            Time in [s] for the entire depletion interval
        source_rate : float
            Power in [W] or source rate in [neutron/sec]
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
        time1, conc1 = self._timed_deplete(
            conc, rates, dt, matrix_func=rk4_f1)
        res1 = self.operator(conc1, source_rate)

        # Step 2: deplete with matrix A(y1) / 2
        time2, conc2 = self._timed_deplete(
            conc, res1.rates, dt, matrix_func=rk4_f1)
        res2 = self.operator(conc2, source_rate)

        # Step 3: deplete with matrix A(y2)
        time3, conc3 = self._timed_deplete(conc, res2.rates, dt)
        res3 = self.operator(conc3, source_rate)

        # Step 4: deplete with matrix built from weighted rates
        list_rates = list(zip(rates, res1.rates, res2.rates, res3.rates))
        time4, conc4 = self._timed_deplete(
            conc, list_rates, dt, matrix_func=rk4_f4)

        return (time1 + time2 + time3 + time4, [conc1, conc2, conc3, conc4],
                [res1, res2, res3])


@add_params
class LEQIIntegrator(Integrator):
    r"""Deplete using the LE/QI CFQ4 algorithm.

    Implements the LE/QI Predictor-Corrector algorithm using the `fourth order
    commutator-free integrator <https://doi.org/10.1137/05063042>`_.

    "LE/QI" stands for linear extrapolation on predictor and quadratic
    interpolation on corrector. This algorithm is mathematically defined as:

    .. math::
        \begin{aligned}
        \mathbf{A}_{-1} &= \mathbf{A}(\mathbf{n}_{i-1}) \\
        \mathbf{A}_0 &= \mathbf{A}(\mathbf{n}_i) \\
        \mathbf{F}_1 &= \frac{-h_i}{12h_{i-1}} \mathbf{A}_{-1} + \frac{6h_{i-1}
            + h_i}{12h_{i-1}} \mathbf{A}_0 \\
        \mathbf{F}_2 &= \frac{-5h_i}{12h_{i-1}} \mathbf{A}_{-1} + \frac{6h_{i-1}
            + 5h_i}{12h_{i-1}} \mathbf{A}_0 \\
        \mathbf{n}_{i+1}^p &= \exp (h_i \mathbf{F}_1) \exp(h_i \mathbf{F}_2)
            \mathbf{n}_i \\
        \mathbf{A}_1 &= \mathbf{A}(\mathbf{n}_{i+1}^p) \\
        \mathbf{F}_3 &= \frac{-h_i^2}{12 h_{i-1} (h_{i-1} + h_i)} \mathbf{A}_{-1} +
              \frac{5 h_{i-1}^2 + 6 h_i h_{i-1} + h_i^2}{12 h_{i-1} (h_{i-1} +
              h_i)} \mathbf{A}_0 + \frac{h_{i-1}}{12 (h_{i-1} + h_i)} \mathbf{A}_1 \\
        \mathbf{F}_4 &= \frac{-h_i^2}{12 h_{i-1} (h_{i-1} + h_i)} \mathbf{A}_{-1} +
              \frac{h_{i-1}^2 + 2 h_i h_{i-1} + h_i^2}{12 h_{i-1} (h_{i-1} + h_i)}
              \mathbf{A}_0 + \frac{5 h_{i-1}^2 + 4 h_i h_{i-1}}{12 h_{i-1}
              (h_{i-1} + h_i)} \mathbf{A}_1 \\
        \mathbf{n}_{i+1} &= \exp(h_i \mathbf{F}_4) \exp(h_i \mathbf{F}_3) \mathbf{n}_i
        \end{aligned}

    It is initialized using the CE/LI algorithm.
    """
    _num_stages = 2

    def __call__(self, bos_conc, bos_rates, dt, source_rate, i):
        """Perform the integration across one time step

        Parameters
        ----------
        conc : numpy.ndarray
            Initial concentrations for all nuclides in [atom]
        rates : openmc.deplete.ReactionRates
            Reaction rates from operator
        dt : float
            Time in [s] for the entire depletion interval
        source_rate : float
            Power in [W] or source rate in [neutron/sec]
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
                    self, bos_conc, bos_rates, dt, source_rate, i)
            prev_res = self.operator.prev_res[-2]
            prev_dt = self.timesteps[i] - prev_res.time[0]
            self._prev_rates = prev_res.rates[0]
        else:
            prev_dt = self.timesteps[i - 1]

        # Remaining LE/QI
        bos_res = self.operator(bos_conc, source_rate)

        le_inputs = list(zip(
            self._prev_rates, bos_res.rates, repeat(prev_dt), repeat(dt)))

        time1, conc_inter = self._timed_deplete(
            bos_conc, le_inputs, dt, matrix_func=leqi_f1)
        time2, conc_eos0 = self._timed_deplete(
            conc_inter, le_inputs, dt, matrix_func=leqi_f2)

        res_inter = self.operator(conc_eos0, source_rate)

        qi_inputs = list(zip(
            self._prev_rates, bos_res.rates, res_inter.rates,
            repeat(prev_dt), repeat(dt)))

        time3, conc_inter = self._timed_deplete(
            bos_conc, qi_inputs, dt, matrix_func=leqi_f3)
        time4, conc_eos1 = self._timed_deplete(
            conc_inter, qi_inputs, dt, matrix_func=leqi_f4)

        # store updated rates
        self._prev_rates = copy.deepcopy(bos_res.rates)

        return (
            time1 + time2 + time3 + time4, [conc_eos0, conc_eos1],
            [bos_res, res_inter])


@add_params
class SICELIIntegrator(SIIntegrator):
    r"""Deplete using the SI-CE/LI CFQ4 algorithm.

    Implements the stochastic implicit CE/LI predictor-corrector algorithm
    using the `fourth order commutator-free integrator
    <https://doi.org/10.1137/05063042>`_.

    Detailed algorithm can be found in section 3.2 in `Colin Josey's thesis
    <http://hdl.handle.net/1721.1/113721>`_.
    """
    _num_stages = 2

    def __call__(self, bos_conc, bos_rates, dt, source_rate, _i=None):
        """Perform the integration across one time step

        Parameters
        ----------
        bos_conc : numpy.ndarray
            Initial bos_concentrations for all nuclides in [atom]
        bos_rates : openmc.deplete.ReactionRates
            Reaction rates from operator
        dt : float
            Time in [s] for the entire depletion interval
        source_rate : float
            Power in [W] or source rate in [neutron/sec]
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
        proc_time, eos_conc = self._timed_deplete(bos_conc, bos_rates, dt)
        inter_conc = copy.deepcopy(eos_conc)

        # Begin iteration
        for j in range(self.n_steps + 1):
            inter_res = self.operator(inter_conc, source_rate)

            if j <= 1:
                res_bar = copy.deepcopy(inter_res)
            else:
                rates = 1/j * inter_res.rates + (1 - 1 / j) * res_bar.rates
                k = 1/j * inter_res.k + (1 - 1 / j) * res_bar.k
                res_bar = OperatorResult(k, rates)

            list_rates = list(zip(bos_rates, res_bar.rates))
            time1, inter_conc = self._timed_deplete(
                bos_conc, list_rates, dt, matrix_func=celi_f1)
            time2, inter_conc = self._timed_deplete(
                inter_conc, list_rates, dt, matrix_func=celi_f2)
            proc_time += time1 + time2

        # end iteration
        return proc_time, [eos_conc, inter_conc], [res_bar]


@add_params
class SILEQIIntegrator(SIIntegrator):
    r"""Deplete using the SI-LE/QI CFQ4 algorithm.

    Implements the Stochastic Implicit LE/QI Predictor-Corrector algorithm
    using the `fourth order commutator-free integrator
    <https://doi.org/10.1137/05063042>`_.

    Detailed algorithm can be found in Section 3.2 in `Colin Josey's thesis
    <http://hdl.handle.net/1721.1/113721>`_.
    """
    _num_stages = 2

    def __call__(self, bos_conc, bos_rates, dt, source_rate, i):
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
        source_rate : float
            Power in [W] or source rate in [neutron/sec]
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
                    self, bos_conc, bos_rates, dt, source_rate, i)
            prev_res = self.operator.prev_res[-2]
            prev_dt = self.timesteps[i] - prev_res.time[0]
            self._prev_rates = prev_res.rates[0]
        else:
            prev_dt = self.timesteps[i - 1]

        # Perform remaining LE/QI
        inputs = list(zip(self._prev_rates, bos_rates,
                          repeat(prev_dt), repeat(dt)))
        proc_time, inter_conc = self._timed_deplete(
            bos_conc, inputs, dt, matrix_func=leqi_f1)
        time1, eos_conc = self._timed_deplete(
            inter_conc, inputs, dt, matrix_func=leqi_f2)

        proc_time += time1
        inter_conc = copy.deepcopy(eos_conc)

        for j in range(self.n_steps + 1):
            inter_res = self.operator(inter_conc, source_rate)

            if j <= 1:
                res_bar = copy.deepcopy(inter_res)
            else:
                rates = 1 / j * inter_res.rates + (1 - 1 / j) * res_bar.rates
                k = 1 / j * inter_res.k + (1 - 1 / j) * res_bar.k
                res_bar = OperatorResult(k, rates)

            inputs = list(zip(self._prev_rates, bos_rates, res_bar.rates,
                              repeat(prev_dt), repeat(dt)))
            time1, inter_conc = self._timed_deplete(
                bos_conc, inputs, dt, matrix_func=leqi_f3)
            time2, inter_conc = self._timed_deplete(
                inter_conc, inputs, dt, matrix_func=leqi_f4)
            proc_time += time1 + time2

        return proc_time, [eos_conc, inter_conc], [res_bar]


integrator_by_name = {
    'cecm': CECMIntegrator,
    'predictor': PredictorIntegrator,
    'cf4': CF4Integrator,
    'epc_rk4': EPCRK4Integrator,
    'si_celi': SICELIIntegrator,
    'si_leqi': SILEQIIntegrator,
    'celi': CELIIntegrator,
    'leqi': LEQIIntegrator
}
