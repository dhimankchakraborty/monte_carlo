import numpy as np
from scipy.linalg import expm
from functions import *
from tqdm import tqdm




lx = 2
ly = 2
n_sites = lx * ly
n_elec_up = 2
n_elec_dn = 2
n_particles = n_elec_up + n_elec_dn
t = 1.0
delta_tau = 0.01
U = 1.0
n_walkers = 500
# n_blocks = 10 # What are blocks?
measurement_flag = 0
E = 0
W = 0
n_equilibrium_blocks = 2
n_measurement_block = 5
n_block_steps = 500
itv_modsvd = 5
itv_pc = 10
itv_Em = 20
E_block = np.zeros(n_measurement_block)
W_block = np.zeros(n_measurement_block)


Phi_t = states(n_sites, n_elec_up, n_elec_dn)
# K = hubbard_1D_ham_kintic(n_sites, pbc=False)
K = hubbard_2D_ham_kinetic(lx, ly, t=1, pbc=False)
K_half_projector = expm(-0.5 * delta_tau * K); 

K_evals, K_evecs = np.linalg.eigh(K)
Phi_t = np.hstack((K_evecs[:, 0:n_elec_up], K_evecs[:, 0:n_elec_dn]))

E_K_t = np.sum(K_evals[0 : n_elec_up]) + np.sum(K_evals[0 : n_elec_dn])

n_r_up = np.diag(Phi_t[:, 0:n_elec_up] @ Phi_t[:, 0:n_elec_up].conj().T)
n_r_dn = np.diag(Phi_t[:, n_elec_up:n_particles] @ Phi_t[:, n_elec_up:n_particles].conjugate().transpose())
E_V_t = U * n_r_up.conj().T @ n_r_dn
E_t = E_K_t + E_V_t


Phi = np.zeros((n_walkers, n_sites, n_particles)) # Walker ensamble
for i in range(n_walkers):
    Phi[i, :, :] = Phi_t.copy()

weights = np.ones(n_walkers)
overlaps = np.ones(n_walkers)


fac_norm = (np.real(E_t) - 0.5 * U * n_particles) * delta_tau # ????
gamma = np.arccosh(np.exp(0.5 * delta_tau * U))
auxiliary_field = np.zeros((2, 2))
for i in range(2):
    for j in range(2):
        auxiliary_field[i, j] = np.exp(gamma * (-1)**(i+j))