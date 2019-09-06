"""Functions to form the special matrix for depletion"""


def celi_f1(chain, rates, fission_yields=None):
    return (5 / 12 * chain.form_matrix(rates[0], fission_yields)
            + 1 / 12 * chain.form_matrix(rates[1], fission_yields))


def celi_f2(chain, rates, fission_yields=None):
    return (1 / 12 * chain.form_matrix(rates[0], fission_yields)
            + 5 / 12 * chain.form_matrix(rates[1], fission_yields))


def cf4_f1(chain, rates, fission_yields=None):
    return 1 / 2 * chain.form_matrix(rates, fission_yields)


def cf4_f2(chain, rates, fission_yields=None):
    return (-1 / 2 * chain.form_matrix(rates[0], fission_yields)
            + chain.form_matrix(rates[1], fission_yields))


def cf4_f3(chain, rates, fission_yields=None):
    return (1 / 4 * chain.form_matrix(rates[0], fission_yields)
            + 1 / 6 * chain.form_matrix(rates[1], fission_yields)
            + 1 / 6 * chain.form_matrix(rates[2], fission_yields)
            - 1 / 12 * chain.form_matrix(rates[3], fission_yields))


def cf4_f4(chain, rates, fission_yields=None):
    return (-1 / 12 * chain.form_matrix(rates[0], fission_yields)
            + 1 / 6 * chain.form_matrix(rates[1], fission_yields)
            + 1 / 6 * chain.form_matrix(rates[2], fission_yields)
            + 1 / 4 * chain.form_matrix(rates[3], fission_yields))


def rk4_f1(chain, rates, fission_yields=None):
    return 1 / 2 * chain.form_matrix(rates, fission_yields)


def rk4_f4(chain, rates, fission_yields=None):
    return (1 / 6 * chain.form_matrix(rates[0], fission_yields)
            + 1 / 3 * chain.form_matrix(rates[1], fission_yields)
            + 1 / 3 * chain.form_matrix(rates[2], fission_yields)
            + 1 / 6 * chain.form_matrix(rates[3], fission_yields))


def leqi_f1(chain, inputs, fission_yields):
    f1 = chain.form_matrix(inputs[0], fission_yields)
    f2 = chain.form_matrix(inputs[1], fission_yields)
    dt_l, dt = inputs[2], inputs[3]
    return -dt / (12 * dt_l) * f1 + (dt + 6 * dt_l) / (12 * dt_l) * f2


def leqi_f2(chain, inputs, fission_yields=None):
    f1 = chain.form_matrix(inputs[0], fission_yields)
    f2 = chain.form_matrix(inputs[1], fission_yields)
    dt_l, dt = inputs[2], inputs[3]
    return -5 * dt / (12 * dt_l) * f1 + (5 * dt + 6 * dt_l) / (12 * dt_l) * f2


def leqi_f3(chain, inputs, fission_yields=None):
    f1 = chain.form_matrix(inputs[0], fission_yields)
    f2 = chain.form_matrix(inputs[1], fission_yields)
    f3 = chain.form_matrix(inputs[2], fission_yields)
    dt_l, dt = inputs[3], inputs[4]
    return (-dt ** 2 / (12 * dt_l * (dt + dt_l)) * f1
            + (dt ** 2 + 6 * dt * dt_l + 5 * dt_l ** 2)
            / (12 * dt_l * (dt + dt_l)) * f2 + dt_l / (12 * (dt + dt_l)) * f3)


def leqi_f4(chain, inputs, fission_yields=None):
    f1 = chain.form_matrix(inputs[0], fission_yields)
    f2 = chain.form_matrix(inputs[1], fission_yields)
    f3 = chain.form_matrix(inputs[2], fission_yields)
    dt_l, dt = inputs[3], inputs[4]
    return (-dt ** 2 / (12 * dt_l * (dt + dt_l)) * f1
            + (dt ** 2 + 2 * dt * dt_l + dt_l ** 2)
            / (12 * dt_l * (dt + dt_l)) * f2
            + (4 * dt * dt_l + 5 * dt_l ** 2) / (12 * dt_l * (dt + dt_l)) * f3)
