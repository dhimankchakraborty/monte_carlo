import numpy as np
from debug_functions import *




def states(n_sites, n_elec_up, n_elec_dn):
    return np.zeros((n_sites, n_elec_up + n_elec_dn))


def hubbard_1D_ham_kintic(n_sites, t=1, pbc=False):
    K = np.zeros((n_sites, n_sites))

    for i in range(n_sites):
        if pbc:
            j = (i + 1) % n_sites
        else:
            if (i + 1) < n_sites:
                j = i + 1
            else:
                continue
        
        K[i, j] = -t
        K[j, i] = -t

    return K


def hubbard_2D_ham_kinetic(lx, ly, t=1, pbc=True):
    n_sites = lx * ly
    k = np.zeros((n_sites, n_sites))

    for y in range(ly):
        for x in range(lx):
            i = x + y * lx
            
            # Neighbors: right, left, up, down
            neighbors = []
            if pbc:
                neighbors.append(((x + 1) % lx, y))
                neighbors.append(((x - 1) % lx, y))
                neighbors.append((x, (y + 1) % ly))
                neighbors.append((x, (y - 1) % ly))
            else:
                if x + 1 < lx: neighbors.append((x + 1, y))
                if x - 1 >= 0: neighbors.append((x - 1, y))
                if y + 1 < ly: neighbors.append((x, y + 1))
                if y - 1 >= 0: neighbors.append((x, y - 1))
            
            for nx, ny in neighbors:
                j = nx + ny * lx
                k[i, j] = -t
                
    return k


def step_walk(Phi, n_walkers, n_sites, weights, overlaps, E, W, K, K_half_projector, measurement_flag, Phi_t, n_elec_up, n_particles, U, fac_norm, auxiliary_field):
    e = np.zeros(n_walkers)

    for walker_idx in range(n_walkers):
        phi = Phi[walker_idx, :, :]
        if weights[walker_idx] > 0:
            weights[walker_idx] = weights[walker_idx] * np.exp(fac_norm)

            phi, weights[walker_idx], overlaps[walker_idx], invO_matrix_up, invO_matrix_dn = apply_K_half_projector(phi, weights[walker_idx], overlaps[walker_idx], K_half_projector, Phi_t, n_elec_up, n_particles)

            if weights[walker_idx] > 0:
                for j_site in range(n_sites):
                    if weights[walker_idx] > 0:

                        phi[j_site,:], overlaps[walker_idx], weights[walker_idx], invO_matrix_up, invO_matrix_dn = apply_interaction_projector(phi[j_site,:], Phi_t[j_site,:], n_elec_up, n_particles, overlaps[walker_idx], weights[walker_idx], invO_matrix_up, invO_matrix_dn, auxiliary_field)
            
            if weights[walker_idx] > 0:
                phi, weights[walker_idx], overlaps[walker_idx], invO_matrix_up, invO_matrix_dn = apply_K_half_projector(phi, weights[walker_idx], overlaps[walker_idx], K_half_projector, Phi_t, n_elec_up, n_particles)

                if weights[walker_idx] > 0:
                    if measurement_flag == 1: 
                        e[walker_idx] = measure_energy(K, phi, Phi_t,  invO_matrix_up, invO_matrix_dn, n_elec_up, n_particles, U)
                
        Phi[walker_idx, :, :] = phi


    if measurement_flag == 1:
        for walker_idx in range(n_walkers):
            if weights[walker_idx] > 0:
                E = E + e[walker_idx] * weights[walker_idx]
                W = W + weights[walker_idx]

    return Phi, weights, overlaps, E, W




def apply_K_half_projector(phi, w, O, K_half_projector, Phi_t, n_elec_up, n_particles):
    phi = K_half_projector @ phi

    invO_matrix_up = np.linalg.inv(Phi_t[:, :n_elec_up].conj().T @ phi[:, :n_elec_up])
    invO_matrix_dn = np.linalg.inv(Phi_t[:, n_elec_up:n_particles].conj().T @ phi[:, n_elec_up:n_particles])

    O_new = 1 / (np.linalg.det(invO_matrix_up) * np.linalg.det(invO_matrix_dn))
    O_ratio = O_new / O

    if O_ratio > 0:
        O = O_new
        w = w * np.real(O_ratio)
    else:
        w = 0

    return phi, w, O, invO_matrix_up, invO_matrix_dn


def apply_interaction_projector(phi_row, Phi_t_row, n_elec_up, n_particles, O, w, invO_matrix_up, invO_matrix_dn, auxiliary_field):

    gii = np.zeros(2)
    rr = np.ones((2, 2))
    matone = rr.copy()

    temp1_up = phi_row[:n_elec_up] @ invO_matrix_up
    temp1_dn = phi_row[n_elec_up:n_particles] @ invO_matrix_dn

    temp2_up = invO_matrix_up @ Phi_t_row[:n_elec_up]
    temp2_dn = invO_matrix_dn @ Phi_t_row[n_elec_up:n_particles]

    gii[0] = np.sum(temp1_up * Phi_t_row[:n_elec_up]) 
    gii[1] = np.sum(temp1_dn * Phi_t_row[n_elec_up:n_particles])

    g_matrix = np.tile(gii, (2, 1)).T # ????
    rr = (auxiliary_field - matone) * g_matrix + matone # ????
    # whos()
    # print(temp1_up.shape)
    # print(temp1_dn.shape)
    # print(temp2_up.shape)
    # print(temp2_dn.shape)
    # print(invO_matrix_up.shape)
    # print(invO_matrix_dn.shape)
    # print(Phi_t_row.shape)
    # print(phi_row.shape)
    # print(gii.shape)
    # print(g_matrix.shape)
    exit()

    o_ratio_temp = rr[0, :] * rr[1, :]
    o_ratio_temp_real = np.maximum(np.real(o_ratio_temp), 0)
    sum_o_ratio_temp_real = np.sum(o_ratio_temp_real)

    if sum_o_ratio_temp_real <= 0:
        w = 0

    if w > 0:
        w = w * 0.5 * sum_o_ratio_temp_real

        if o_ratio_temp_real[0] / sum_o_ratio_temp_real >= np.random.rand():
            x_spin = 0
        else:
            x_spin = 1

        # print(temp1_up.shape)
        # print(temp2_up.shape)
        # print(np.outer(temp2_up, temp1_up).shape)
        # exit()
        phi_row[:n_elec_up] *= auxiliary_field[0, x_spin]
        phi_row[n_elec_up:n_particles] *= auxiliary_field[1, x_spin]

        O *= o_ratio_temp[x_spin]

        invO_matrix_up += ((1 - auxiliary_field[0, x_spin]) / rr[0, x_spin]) * np.outer(temp2_up, temp1_up)
        invO_matrix_dn += ((1 - auxiliary_field[1, x_spin]) / rr[1, x_spin]) * np.outer(temp2_dn, temp1_dn)

    return phi_row, O, w, invO_matrix_up, invO_matrix_dn


def measure_energy(K, phi, Phi_t,  invO_matrix_up, invO_matrix_dn, n_elec_up, n_particles, U):

    temp_up = phi[:, :n_elec_up] @ invO_matrix_up
    temp_dn = phi[:, n_elec_up:n_particles] @ invO_matrix_dn

    g_up = temp_up @ Phi_t[:, :n_elec_up].conj().T
    g_dn = temp_dn @ Phi_t[:, n_elec_up:n_particles].conj().T

    n_up = np.diag(g_up)
    n_dn = np.diag(g_dn)
    n_int = np.dot(n_up, n_dn)
    potential_energy = n_int * U

    kinetic_energy = np.sum(K.T * (g_up + g_dn))

    e = potential_energy + kinetic_energy

    # print(f"KE: {kinetic_energy}")
    # print(f"PE: {potential_energy}")

    return e


def stabilize_walkers(Phi, n_walkers, overlaps, n_elec_up, n_particles):
    for walker_idx in range(n_walkers):
        Phi[walker_idx, :, :n_elec_up], r_up = np.linalg.qr(Phi[walker_idx, :, :n_elec_up])
        Phi[walker_idx, :, n_elec_up:n_particles], r_dn = np.linalg.qr(Phi[walker_idx, :, n_elec_up:n_particles])

        _, logabsdet_up = np.linalg.slogdet(r_up)
        _, logabsdet_dn = np.linalg.slogdet(r_dn)
        overlaps[walker_idx] /= np.exp(logabsdet_up + logabsdet_dn)

        # overlaps[walker_idx] /= (np.linalg.det(r_up) * np.linalg.det(r_dn))

    return Phi, overlaps


def population_control(Phi, weights, overlaps, n_walkers, n_sites, n_particles): # ????
    new_Phi = np.zeros((n_walkers, n_sites, n_particles))
    new_overlaps = np.zeros(n_walkers)

    d = n_walkers / np.sum(weights)

    sum_w = -np.random.rand() # ????
    n_wlk_resampled = 0

    for walker_idx in range(n_walkers):
        sum_w += weights[walker_idx] * d
        n = int(np.ceil(sum_w))
        n = min(n, n_walkers)
        
        for j in range(n_wlk_resampled, n):
            new_Phi[j, :, :] = Phi[walker_idx, :, :]
            new_overlaps[j] = overlaps[walker_idx]
            
        n_wlk_resampled = n

    weights = np.ones(n_walkers)

    return new_Phi, weights, new_overlaps