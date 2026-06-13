import sys
import os
from pprint import pprint
import numpy as np
from itertools import product
from scipy.linalg import expm

sys.path.append(os.path.dirname(os.getcwd()))
# from src.core import * # type: ignore
from src.hamiltonians import *
from src.lattices import *
from src.propagators import *
from src.states import *
from src.walkers import *


# ==========================================================
# PARAMETERS
# ==========================================================

L = 2

t = 1.0
U = 4.0

dtau = 0.05


# ==========================================================
# SYSTEM
# ==========================================================

lattice = Chain(
    n_sites=L,
    pbc=False
)

system = HubbardSystem(
    lattice=lattice,
    t=t,
    U=U
)

K = system.h_kin

prop = HubbardPropagator(
    K=K,
    U=U,
    dtau=dtau
)

print("gamma =", prop.gamma)

# ==========================================================
# EXACT INTERACTION PROPAGATOR
# ==========================================================

# local basis:
#
# |0>
# |up>
# |down>
# |updown>
#

local_V = np.diag([
    0.0,
    0.0,
    0.0,
    U
])

exact_local = expm(
    -dtau * local_V
)

print()
print("Exact local interaction propagator")
print(exact_local)

# ==========================================================
# HS AVERAGE
# ==========================================================

gamma = prop.gamma

hs_average = np.zeros((4, 4))

for x in [-1, +1]:

    # occupation numbers in local basis
    #
    # |0>         : (0,0)
    # |up>        : (1,0)
    # |down>      : (0,1)
    # |up down>   : (1,1)

    n_up = np.diag([0, 1, 0, 1])
    n_dn = np.diag([0, 0, 1, 1])

    hs_matrix = (
        expm(
            -0.5 * dtau * U *
            (n_up + n_dn)
        )
        @
        expm(
            gamma * x *
            (n_up - n_dn)
        )
    )

    hs_average += 0.5 * hs_matrix

print()
print("HS reconstructed propagator")
print(hs_average)

print()
print("Difference")
# print(exact_local)
# print(hs_average)

print(
    np.max(
        np.abs(
            exact_local - hs_average
        )
    )
)