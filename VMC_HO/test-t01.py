import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from functions import *




rng = np.random.default_rng()
rng_step = np.random.default_rng()
rng_acceptance = np.random.default_rng()

x_init = rng.random() * 2 - 1
step_size = 0.5
n_mc_steps = 10**7
n_therm_steps = int(n_mc_steps / 10)

beta = 0.51

x_current = x_init
for i in tqdm(range(n_therm_steps)):
    x_next = x_current + (rng_step.random() * 2 * step_size - step_size)
    acceptance_ratio = acceptance_ratio_t01(beta, x_current, x_next)

    # if (acceptance_ratio >= 1): # CHECK
    #     x_current = x_next

    if acceptance_ratio > rng_acceptance.random():
        x_current = x_next

total_energy = 0
total_energy_squared = 0

for i in tqdm(range(n_mc_steps)):
    x_next = x_current + (rng_step.random() * 2 * step_size - step_size)
    acceptance_ratio = acceptance_ratio_t01(beta, x_current, x_next)

    # if (acceptance_ratio >= 1): # CHECK
    #     x_current = x_next
    #     total_energy += local_energy_t01(beta, x_current) * probability_density_t01(beta, x_current)

    if acceptance_ratio > rng_acceptance.random():
        x_current = x_next
    
    local_energy = local_energy_t01(beta, x_current)
    total_energy += local_energy
    total_energy_squared += local_energy**2

energy = total_energy / n_mc_steps
varience = (total_energy_squared / n_mc_steps) - energy**2
standard_deviation = np.sqrt(varience)

print("Energy:", energy)
print("Varience:", varience)
print("Standard Deviation:", standard_deviation)