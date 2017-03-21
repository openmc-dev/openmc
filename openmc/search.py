from collections import Callable
from numbers import Real

import openmc
import openmc.model
import openmc.checkvalue as cv


_SCALAR_BRACKETED_METHODS = ['brentq', 'brenth', 'ridder', 'bisect']


class KeffSearch(object):
    """Class to perform a keff search by modifying a model parametrized by a
    single independent variable.

    Parameters
    ----------
    model_builder : collections.Callable
        Callable function which builds a model according to a passed
        parameter. This function must return an openmc.model.Model object.
    guess : Real, optional
        Initial guess for the parameter to be searched in
        :param:`model_builder`. One of :param:`guess` or :param`bracket` must
        be provided.
    target_keff : Real, optional
        keff value to search for, defaults to 1.0.
    bracket : None or Iterable of Real, optional
        Bracketing interval to search for the solution; if not provided,
        a generic non-bracketing method is used. If provided, the brackets
        are used. Defaults to no brackets provided. One of :param:`guess` or
        :param`bracket` must be provided. If both are provided, the bracket
        will be preferentially used.

    Attributes
    ----------
    model_builder : FunctionType
        Callable function which builds a model according to parameters passed.
        This function must return an openmc.model.Model
        object.
    initial_guess : Iterable of Real
        Initial guess for the parameter to be searched in
        :param:`model_builder`.
    target_keff : Real
        keff value to search for.
    bracket : None or Iterable of Real
        Bracketing interval to search for the solution; if not provided,
        a generic non-bracketing method is used. If provided, the brackets
        are used.
    model_args : dict
        Keyword-based arguments to pass to the :param:`model_builder` method.
    print_iterations : bool
        Whether or not to print the guess and the resultant keff during the
        iteration process. Defaults to False.
    print_output : bool
        Whether or not to print the OpenMC output during the iterations.
        Defaults to False.
    guesses : List of Real
        List of guesses attempted by the search
    keffs : List of Real
        List of keffs corresponding to the guess attempted by the search
    keff_uncs : List of Real
        List of keff uncertainties corresponding to the guess attempted by the
        search

    """

    def __init__(self, model_builder, guess=None, bracket=None,
                 target_keff=1.0, model_args={}):
        self.initial_guess = guess
        self.bracket = bracket
        self.model_args = model_args
        self.model_builder = model_builder
        self.target_keff = target_keff
        self.guesses = []
        self.keffs = []
        self.keff_uncs = []
        self.print_iterations = False
        self.print_output = False

    @property
    def model_builder(self):
        return self._model_builder

    @property
    def initial_guess(self):
        return self._initial_guess

    @property
    def bracket(self):
        return self._bracket

    @property
    def target_keff(self):
        return self._target_keff

    @property
    def print_iterations(self):
        return self._print_iterations

    @property
    def print_output(self):
        return self._print_output

    @property
    def model_args(self):
        return self._model_args

    @model_builder.setter
    def model_builder(self, model_builder):
        # Make sure model_builder is a function
        cv.check_type('model_builder', model_builder, Callable)

        # Run the model builder function once to make sure it provides the
        # correct output type
        if self.bracket is not None:
            model = model_builder(self.bracket[0], **self.model_args)
        elif self.initial_guess is not None:
            model = model_builder(self.initial_guess, **self.model_args)
        cv.check_type('model_builder return', model, openmc.model.Model)
        self._model_builder = model_builder

    @initial_guess.setter
    def initial_guess(self, initial_guess):
        if initial_guess is not None:
            cv.check_type('initial_guess', initial_guess, Real)
        self._initial_guess = initial_guess

    @bracket.setter
    def bracket(self, bracket):
        if bracket is not None:
            cv.check_iterable_type('bracket', bracket, Real)
            cv.check_length('bracket', bracket, 2)
        self._bracket = bracket

    @target_keff.setter
    def target_keff(self, target_keff):
        cv.check_type('target_keff', target_keff, Real)
        self._target_keff = target_keff

    @print_iterations.setter
    def print_iterations(self, print_iterations):
        cv.check_type('print_iterations', print_iterations, bool)
        self._print_iterations = print_iterations

    @print_output.setter
    def print_output(self, print_output):
        cv.check_type('print_output', print_output, bool)
        self._print_output = print_output

    @model_args.setter
    def model_args(self, model_args):
        cv.check_type('model_args', model_args, dict)
        self._model_args = model_args

    def _search_function(self, guess):
        # Build the model
        model = self.model_builder(guess, **self.model_args)

        # Run the model
        keff = model.run(output=self._print_output)

        # Close the model to ensure HDF5 will allow access during the next
        # OpenMC execution
        model.close()

        # Record the history
        self.guesses.append(guess)
        self.keffs.append(keff[0])
        self.keff_uncs.append(keff[1])

        if self._print_iterations:
            text = 'Iteration: {}; Guess of {:.2e} produced a keff of ' + \
                '{:1.5f} +/- {:1.5f}'
            print(text.format(self._i, guess, keff[0], keff[1]))
        self._i += 1

        return (keff[0] - self.target_keff)

    def search(self, tol=None, bracketed_method='brentq', **kwargs):
        """Apply root-finding algorithm to search for the target k-eigenvalue

        Parameters
        ----------
        tol : float
            Tolerance to pass to the search method
        bracketed_method : {'brentq', 'brenth', 'ridder', 'bisect'}, optional
            Solution method to use; only applies if
            :param:`bracket` is set, otherwise the Secant method is used.
            Defaults to 'brentq'.
        **kwargs
            All remaining keyword arguments are passed to the root-finding
            method.

        Returns
        -------
        zero_value : float
            Estimated value of the variable parameter where keff is the
            targeted value

        """

        cv.check_value('bracketed_method', bracketed_method,
                       _SCALAR_BRACKETED_METHODS)

        import scipy.optimize as sopt

        if self.bracket is not None:
            # Generate our arguments
            args = {'f': self._search_function, 'a': self.bracket[0],
                    'b': self.bracket[1]}
            if tol is not None:
                args['rtol'] = tol

            # Set the root finding method
            if bracketed_method == 'brentq':
                root_finder = sopt.brentq
            elif bracketed_method == 'brenth':
                root_finder = sopt.brenth
            elif bracketed_method == 'ridder':
                root_finder = sopt.ridder
            elif bracketed_method == 'bisect':
                root_finder = sopt.bisect

        elif self.initial_guess is not None:

            # Generate our arguments
            args = {'func': self._search_function, 'x0': self.initial_guess}
            if tol is not None:
                args['tol'] = tol

            # Set the root finding method
            root_finder = sopt.newton

        else:
            raise ValueError("One of the 'bracket' or 'initial_guess' "
                             "parameters must be set")

        # Set the iteration counter
        self._i = 0

        # Create a new dictionary with the arguments from args and kwargs
        args.update(kwargs)

        # Perform the search
        zero_value = root_finder(**args)

        return zero_value
