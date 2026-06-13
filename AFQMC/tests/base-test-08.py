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

lattice = Chain(n_sites=L, pbc=False)
system = HubbardSystem(lattice=lattice, t=t, U=U)
prop = HubbardPropagator(K=system.h_kin, U=U, dtau=dtau)
trial = SlaterDeterminantTwoSpinState(hamiltonian=system, n_electrons_up=Nup, n_electrons_down=Ndown)
trial.initialize("non-interacting")


walkers = []
for i in range(n_walkers):
    state = SlaterDeterminantTwoSpinState(hamiltonian=system, n_electrons_up=Nup, n_electrons_down=Ndown)
    state.initialize("random")
    walker = Walker(state)
    walker.overlap = walker.compute_overlap(trial)
    walkers.append(walker)

    # print(i, "    ", system.calculate_variational_energy(walker))


for step in range(n_steps):
    for walker in walkers:
        prop.propagate(walker)
        walker.orthogonalize()
        walker.update_weight(trial)
        # print(walker.weight)

        # new_overlap = walker.compute_overlap(trial)

        # if abs(walker.overlap) < 1e-14: # type: ignore
        #     walker.weight = 0.0
        #     walker.overlap = new_overlap

        # walker.weight *= abs(new_overlap / walker.overlap)
        # walker.overlap = new_overlap
        # print(walker.weight)

    if step % 50 == 0:

        energies = [system.calculate_variational_energy(w) for w in walkers]
        weights = np.array([walker.weight.real for walker in walkers])
        E = np.sum(weights * energies) / np.sum(weights)
        print(step, E)


energies = [system.calculate_variational_energy(w) for w in walkers]
weights = np.array([walker.weight.real for walker in walkers])
E = np.sum(weights * energies) / np.sum(weights)

print()
print("Monte Carlo average:")
print(E)

print()
print("Exact Hubbard dimer:")
print(-0.828427124746)


for i in range(len(walkers)):
    print(f"{i}    {system.calculate_variational_energy(walkers[i])}    {walkers[i].weight}")

print("max =", np.max(weights))
print("min =", np.min(weights))
print("mean =", np.mean(weights))