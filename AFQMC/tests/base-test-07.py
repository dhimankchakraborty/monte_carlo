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



L = 2

Nup = 1
Ndown = 1

t = 1.0
U = 4.0

dtau = 0.05

n_steps = 200

n_walkers = 500


# ==================================================
# SYSTEM
# ==================================================

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

trial = SlaterDeterminantTwoSpinState(
    hamiltonian=system,
    n_electrons_up=Nup,
    n_electrons_down=Ndown
)

trial.initialize("non-interacting")

walkers = []

for _ in range(n_walkers):

    state = SlaterDeterminantTwoSpinState(
        hamiltonian=system,
        n_electrons_up=Nup,
        n_electrons_down=Ndown
    )

    state.initialize("random")

    walker = Walker(state)

    walker.compute_overlap(trial)

    walkers.append(walker)

# ==================================================
# PROPAGATION
# ==================================================

for step in range(n_steps):

    for walker in walkers:

        old_overlap = walker.overlap

        prop.propagate(walker)

        walker.orthogonalize()

        # walker.update_weight(trial)

        new_overlap = walker.compute_overlap(
            trial
        )

        if abs(old_overlap) > 1e-14:

            walker.weight *= abs(
                new_overlap / old_overlap
            )

        else:

            walker.weight = 0.0

    if step % 50 == 0:

        energies = [
            system.calculate_energy(w)
            for w in walkers
        ]

        weights = np.array([
            walker.weight.real
            for walker in walkers
        ])

        E = np.sum(weights * energies) / np.sum(weights)

        print(
            step,
            E
        )

# ==================================================
# FINAL RESULT
# ==================================================

energies = [
    system.calculate_energy(w)
    for w in walkers
]

weights = np.array([
    walker.weight.real
    for walker in walkers
])

E = np.sum(weights * energies) / np.sum(weights)


print()
print("Monte Carlo average:")
print(E)

print()
print("Exact Hubbard dimer:")
print(-0.828427124746)


weights = np.array(
    [walker.weight for walker in walkers]
)

print("max =", np.max(weights))
print("min =", np.min(weights))
print("mean =", np.mean(weights))

for i, walker in enumerate(walkers):
    print(f"Energy of walker {i} : {system.calculate_energy(walker)}")