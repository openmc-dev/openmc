
import copy
import numpy as np

TIME_POINTS = ['START',
               'PREVIOUS_OUTER',
               'FORWARD_OUTER',
               'PREVIOUS_INNER',
               'FORWARD_INNER',
               'END']


class Clock(object):

    def __init__(self, start=0., end=3., dt_inner=None, t_outer=None):

        # Initialize coordinates
        self.dt_inner = dt_inner
        self.t_outer = t_outer

        # Create a dictionary of clock times
        self._times = {}
        for t in TIME_POINTS:
            self._times[t] = start

        # Reset the end time
        self._times['END'] = end
        self._outer_step = 0

    def __repr__(self):

        string = 'Clock\n'
        string += '{0: <24}{1}{2}\n'.format('\tdt inner', '=\t', self.dt_inner)
        string += '{0: <24}{1}{2}\n'.format('\tt outer', '=\t', self.t_outer)

        for t in TIME_POINTS:
            string += '{0: <24}{1}{2}\n'.format('\tTime ' + t, '=\t', self.times[t])

        return string

    @property
    def t_outer(self):
        return self._t_outer

    @property
    def times(self):
        return self._times

    @property
    def outer_step(self):
        return self._outer_step

    @property
    def num_outer_steps(self):
        return len(self.t_outer) - 1

    @property
    def dt_inner(self):
        return self._dt_inner

    @dt_inner.setter
    def dt_inner(self, dt_inner):
        self._dt_inner = np.float64(dt_inner)

    @t_outer.setter
    def t_outer(self, t_outer):
        self._t_outer = np.float64(t_outer)

    @property
    def dt_outer(self):
        return self.t_outer[self.outer_step + 1] - self.t_outer[self.outer_step]

    @times.setter
    def times(self, times):
        self._times = np.float64(times)

    @outer_step.setter
    def outer_step(self, outer_step):
        self._outer_step = outer_step
