import numpy as np




def site_index(x, y, length_x):
    return x + y * length_x


def non_interacting_fermi_hubbard_ed(length_x, length_y, t, use_pbc=False):
    H_dim = length_x * length_y
    H = np.zeros((H_dim, H_dim))

    for x in range(length_x):
        for y in range(length_y):

            i = site_index(x, y, length_x)

            x_plus = x + 1
            if x_plus < length_x:
                j = site_index(x_plus, y, length_x)
                H[i, j] = -t
                H[j, i] = -t
            elif use_pbc:
                j = site_index(0, y, length_x)
                H[i, j] = -t
                H[j, i] = -t

            y_plus = y + 1
            if y_plus < length_y:
                j = site_index(x, y_plus, length_x)
                H[i, j] = -t
                H[j, i] = -t
            elif use_pbc:
                # Wrap around
                j = site_index(x, 0, length_x)
                H[i, j] = -t
                H[j, i] = -t

    for i in range(H_dim):
        H[i, i] = 0
    
    eigenvalues, eigenvectors = np.linalg.eig(H)

    return H, eigenvalues, eigenvectors.transpose()


def hubbard_kinetic_2d_matrix(length_x, length_y, t, use_pbc=False):
    dim = length_x * length_y
    K = np.zeros((dim, dim))

    for x in range(length_x):
        for y in range(length_y):

            i = site_index(x, y, length_x)

            x_plus = x + 1
            if x_plus < length_x:
                j = site_index(x_plus, y, length_x)
                K[i, j] = -t
                K[j, i] = -t
            elif use_pbc:
                j = site_index(0, y, length_x)
                K[i, j] = -t
                K[j, i] = -t

            y_plus = y + 1
            if y_plus < length_y:
                j = site_index(x, y_plus, length_x)
                K[i, j] = -t
                K[j, i] = -t
            elif use_pbc:
                # Wrap around
                j = site_index(x, 0, length_x)
                K[i, j] = -t
                K[j, i] = -t

    for i in range(dim):
        K[i, i] = 0
    
    return K



