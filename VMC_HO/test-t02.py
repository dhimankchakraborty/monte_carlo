import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from functions import *






x_init = None
step_size = 0.5
n_mc_steps = 10**7
n_therm_steps = int(n_mc_steps / 10)
run_count = 101

beta_arr = np.linspace(0, 1, run_count)
energy_arr = np.zeros((run_count))

for j, beta in enumerate(beta_arr):
    
    rng_step = np.random.default_rng()
    rng_acceptance = np.random.default_rng()

    x_current = 0
    if x_init == None:
        rng = np.random.default_rng()
        x_current = rng.random() * 2 - 1
    
    for i in tqdm(range(n_therm_steps)):
        x_next = x_current + (rng_step.random() * 2 * step_size - step_size)
        acceptance_ratio = acceptance_ratio_t01(beta, x_current, x_next)

        # if (acceptance_ratio >= 1): # CHECK
        #     x_current = x_next

        if acceptance_ratio > rng_acceptance.random():
            x_current = x_next

    total_energy = 0

    for i in tqdm(range(n_mc_steps)):
        x_next = x_current + (rng_step.random() * 2 * step_size - step_size)
        acceptance_ratio = acceptance_ratio_t01(beta, x_current, x_next)

        # if (acceptance_ratio >= 1): # CHECK
        #     x_current = x_next
        #     total_energy += local_energy_t01(beta, x_current) * probability_density_t01(beta, x_current)

        if acceptance_ratio > rng_acceptance.random():
            x_current = x_next
        
        total_energy += local_energy_t01(beta, x_current)

    energy = total_energy / n_mc_steps
    print(rf"$\beta$ = {beta} ==> Energy = {energy}")

    energy_arr[j] = energy


plt.plot(energy_arr, marker='o')
plt.grid()
plt.savefig('result-2.png')

print("-----------------------------------------------------")
for i in range(run_count):
    print(rf"$\beta$ = {beta_arr[i]} ==> Energy = {energy_arr}")