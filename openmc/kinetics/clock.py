
import copy
import numpy as np

TIME_POINTS = ['START',
               'PREVIOUS_OUT',
               'PREVIOUS_IN',
               'CURRENT',
               'FORWARD_IN',
               'FORWARD_OUT',
               'END']

class Clock(object):

    def __init__(self, start=0., end=3., dt_outer=1.e-1, dt_inner=1.e-2):

        # Initialize coordinates
        self.dt_outer = dt_outer
        self.dt_inner = dt_inner

        # Create a dictionary of clock times
        self._times = {}
        for t in TIME_POINTS:
            self._times[t] = start

        # Reset the end time
        self._times['END'] = end


    def __deepcopy__(self, memo):

        existing = memo.get(id(self))

        # If this is the first time we have tried to copy this object, create a copy
        if existing is None:

            clone = type(self).__new__(type(self))

            memo[id(self)] = clone

            return clone

        # If this object has been copied before, return the first copy made
        else:
            return existing

    def __repr__(self):

        string = 'Clock\n'
        string += '{0: <24}{1}{2}\n'.format('\tdt inner', '=\t', self.dt_inner)
        string += '{0: <24}{1}{2}\n'.format('\tdt outer', '=\t', self.dt_outer)

        for t in TIME_POINTS:
            string += '{0: <24}{1}{2}\n'.format('\tTime ' + t, '=\t', self.times[t])

        return string

    @property
    def dt_inner(self):
        return self._dt_inner

    @property
    def dt_outer(self):
        return self._dt_outer

    @property
    def times(self):
        return self._times

    @dt_inner.setter
    def dt_inner(self, dt_inner):
        self._dt_inner = np.float64(dt_inner)

    @dt_outer.setter
    def dt_outer(self, dt_outer):
        self._dt_outer = np.float64(dt_outer)

    @times.setter
    def times(self, times):
        self._times = np.float64(times)

    def take_outer_step(self):
        """Take an outer time step and reset all the inner time step values
        to the starting point for the outer time step.

        """

        self.times['PREVIOUS_OUT'] = self.times['FORWARD_OUT']
        self.times['PREVIOUS_IN']  = self.times['FORWARD_OUT']
        self.times['FORWARD_IN']   = self.times['FORWARD_OUT']
        self.times['CURRENT']      = self.times['FORWARD_OUT']

        if (self.times['END'] > self.times['FORWARD_OUT'] + self.dt_outer):
            self.times['FORWARD_OUT'] = self.times['END']
        else:
            self.times['FORWARD_OUT'] = self.times['FORWARD_OUT'] + self.dt_outer

    def take_inner_step(self):
        """Take an inner time step.

        """
        self.times['PREVIOUS_IN']  = self.times['FORWARD_IN']
        self.times['CURRENT']      = self.times['FORWARD_IN']

        if (self.times['FORWARD_OUT'] > self.times['FORWARD_IN'] + self.dt_inner):
            self.times['FORWARD_IN'] = self.times['FORWARD_OUT']
        else:
            self.times['FORWARD_IN'] = self.times['FORWARD_IN'] + self.dt_inner

    def reset_to_previous_outer(self):
        """Reset the time values to the previous outer time.

        """

        self.times['PREVIOUS_IN']  = self.times['PREVIOUS_IN']
        self.times['FORWARD_IN']   = self.times['PREVIOUS_IN']
        self.times['CURRENT']      = self.times['PREVIOUS_IN']
