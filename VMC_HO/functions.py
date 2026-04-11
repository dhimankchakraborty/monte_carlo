import numpy as np
import matplotlib.pyplot as plt




def trial_wave_function_t01(beta, x):
    return np.exp(-1 * beta * x * x)


def local_energy_t01(beta, x):
    return beta + (0.5 - 2 * beta * beta) * x * x


def probability_density_t01(beta, x):
    return (trial_wave_function_t01(beta, x) * trial_wave_function_t01(beta, x)) / np.power(2 * beta / np.pi, 0.25)


def acceptance_ratio_t01(beta, x_current, x_next):
    return (trial_wave_function_t01(beta, x_next) * trial_wave_function_t01(beta, x_next) / (trial_wave_function_t01(beta, x_current) * trial_wave_function_t01(beta, x_current)))
