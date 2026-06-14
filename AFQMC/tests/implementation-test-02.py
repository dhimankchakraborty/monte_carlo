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
# from src.propagators import *
from src.states import *
from src.walkers import *




class HubbardPropagator:

    def __init__(self, hamiltonian, dtau):
        self.ham = hamiltonian
        self.K = self.ham.h_kin
        self.U = self.ham.U
        self.dtau = dtau
        self.BK_half = expm(-0.5 * dtau * self.K)
        self.gamma = np.arccosh(np.exp(0.5 * dtau * self.U))
        self.prefector = 0.5 * np.exp(-0.5 * self.dtau * self.U)


    def propagate_walker_exact(self, walker):
        n_sites = walker.state.phi_up.shape[0]
        walker.state.phi_up = self.BK_half @ walker.state.phi_up
        walker.state.phi_down = self.BK_half @ walker.state.phi_down
        print('-'*80)
        print(walker.get_phi_up())
        print(walker.get_phi_down())

        # print(n_sites)
        
        for i in range(n_sites):
            walker.state.phi_up[i] = self.prefector * (np.exp(self.gamma) + np.exp(-1 * self.gamma)) * walker.state.phi_up[i]
            walker.state.phi_down[i] = self.prefector * (np.exp(-1 * self.gamma) + np.exp(self.gamma)) * walker.state.phi_down[i]        
        print('-'*80)
        print(walker.get_phi_up())
        print(walker.get_phi_down())

        print('-'*80)
        print(self.gamma)
        print(self.prefector)
        print((np.exp(self.gamma) + np.exp(-1 * self.gamma)))
        print(self.prefector * (np.exp(self.gamma) + np.exp(-1 * self.gamma)))

        


        walker.state.phi_up = self.BK_half @ walker.state.phi_up
        walker.state.phi_down = self.BK_half @ walker.state.phi_down
        



lattice = Chain(n_sites=2, pbc=False)
hamiltonian = HubbardSystem(lattice=lattice, t=1, U=1000)
trial_state = SlaterDeterminantTwoSpinState(hamiltonian, n_electrons_up=1, n_electrons_down=1)
trial_state.initialize("non-interacting")
state = SlaterDeterminantTwoSpinState(hamiltonian, n_electrons_up=1, n_electrons_down=1)
state.initialize("non-interacting")
walker = Walker(state)
print(walker.get_phi_up())
print(walker.get_phi_down())
    

propagator = HubbardPropagator(hamiltonian, dtau=0.01)

for _ in range(1):
    propagator.propagate_walker_exact(walker)
    walker.orthogonalize()

print('-'*80)
print(hamiltonian.calculate_variational_energy(walker))
print(hamiltonian.energy_mixed_estimator(trial_state, walker))
print('-'*80)
print(walker.get_phi_up())
print(walker.get_phi_down())