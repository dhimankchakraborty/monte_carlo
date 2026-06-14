import numpy as np
from functions import *
from initialization import * 




# print(K)

measurement_phase = True

print("-"*80)
print("Equilibrium Phase")
print("-"*80)
for i_block in range(n_equilibrium_blocks):
    print(f"Block: {i_block + 1}")
    for j_step in tqdm(range(n_block_steps)):
        # print(f"        Step: {j_step}")
        Phi, weights, overlaps, E, W = step_walk(Phi, n_walkers, n_sites, weights, overlaps, E, W, K, K_half_projector, measurement_flag, Phi_t, n_elec_up, n_particles, U, fac_norm, auxiliary_field)

        if j_step % itv_modsvd == 0:
            Phi, overlaps = stabilize_walkers(Phi, n_walkers, overlaps, n_elec_up, n_particles)

        if j_step % itv_pc == 0:
            Phi, weights, overlaps = population_control(Phi, weights, overlaps, n_walkers, n_sites, n_particles)

measurement_phase = True
print("-"*80)
print("Measurement Phase")
print("-"*80)
for i_block in range(n_measurement_block):
    print(f"Block: {i_block + 1}")
    for j_step in tqdm(range(n_block_steps)):
        # print(f"        Step: {j_step}")
        if j_step % itv_Em == 0:
            measurement_flag = 1
        else:
            measurement_flag = 0

        Phi, weights, overlaps, E_block[i_block], W_block[i_block] = step_walk(Phi, n_walkers, n_sites, weights, overlaps, E_block[i_block], W_block[i_block], K, K_half_projector, measurement_flag, Phi_t, n_elec_up, n_particles, U, fac_norm, auxiliary_field)

        if j_step % itv_modsvd == 0:
            Phi, overlaps = stabilize_walkers(Phi, n_walkers, overlaps, n_elec_up, n_particles)

        if j_step % itv_pc == 0:
            Phi, weights, overlaps = population_control(Phi, weights, overlaps, n_walkers, n_sites, n_particles)

        if j_step % itv_Em == 0:
            fac_norm=(np.real(E_block[i_block] / W_block[i_block]) - 0.5 * U * n_particles) * delta_tau

    E_block[i_block] = E_block[i_block] / W_block[i_block]
    print(f"E[{i_block + 1}] = {np.real(E_block[i_block])}")

E = np.real(E_block);
E_average = np.mean(E)
E_error = np.std(E, ddof=1) / np.sqrt(n_measurement_block)

print(f"Average Energy                  : {E_average}")
print(f"Statistical Error in Energy     : {E_error}")



        
