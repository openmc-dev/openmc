import numpy as np


def linearize(x, f, tolerance=0.001):
    """Return a tabulated representation of a one-variable function

    Parameters
    ----------
    x : Iterable of float
        Initial x values at which the function should be evaluated
    f : Callable
        Function of a single variable
    tolerance : float
        Tolerance on the interpolation error

    Returns
    -------
    numpy.ndarray
        Tabulated values of the independent variable
    numpy.ndarray
        Tabulated values of the dependent variable

    """
    # Make sure x is a numpy array
    x = np.asarray(x)

    # Initialize output arrays
    x_out = []
    y_out = []

    # Initialize stack
    x_stack = [x[0]]
    y_stack = [f(x[0])]

    for i in range(x.shape[0] - 1):
        x_stack.insert(0, x[i + 1])
        y_stack.insert(0, f(x[i + 1]))

        while True:
            x_high, x_low = x_stack[-2:]
            y_high, y_low = y_stack[-2:]
            x_mid = 0.5*(x_low + x_high)
            y_mid = f(x_mid)

            y_interp = y_low + (y_high - y_low)/(x_high - x_low)*(x_mid - x_low)
            error = abs((y_interp - y_mid)/y_mid)
            if error > tolerance:
                x_stack.insert(-1, x_mid)
                y_stack.insert(-1, y_mid)
            else:
                x_out.append(x_stack.pop())
                y_out.append(y_stack.pop())
                if len(x_stack) == 1:
                    break

    x_out.append(x_stack.pop())
    y_out.append(y_stack.pop())

    return np.array(x_out), np.array(y_out)

def thin(x, y, tolerance=0.001):
    """Check for (x,y) points that can be removed.

    Parameters
    ----------
    x : numpy.ndarray
        Independent variable
    y : numpy.ndarray
        Dependent variable
    tolerance : float
        Tolerance on interpolation error

    Returns
    -------
    numpy.ndarray
        Tabulated values of the independent variable
    numpy.ndarray
        Tabulated values of the dependent variable

    """
    # Initialize output arrays
    x_out = x.copy()
    y_out = y.copy()

    N = x.shape[0]
    i_left = 0
    i_right = 2

    while i_left < N - 2 and i_right < N:
        m = (y[i_right] - y[i_left])/(x[i_right] - x[i_left])

        for i in range(i_left + 1, i_right):
            # Determine error in interpolated point
            y_interp = y[i_left] + m*(x[i] - x[i_left])
            if abs(y[i]) > 0.:
                error = abs((y_interp - y[i])/y[i])
            else:
                error = 2*tolerance

            if error > tolerance:
                for i_remove in range(i_left + 1, i_right - 1):
                    x_out[i_remove] = np.nan
                    y_out[i_remove] = np.nan
                i_left = i_right - 1
                i_right = i_left + 1
                break

        i_right += 1

    for i_remove in range(i_left + 1, i_right - 1):
        x_out[i_remove] = np.nan
        y_out[i_remove] = np.nan

    return x_out[np.isfinite(x_out)], y_out[np.isfinite(y_out)]
