import numpy as np
from scipy.linalg import expm
# from ./walkers import *




class HubbardPropagator:

    def __init__(self, hamiltonian, dtau):

        self.K = hamiltonian.h_kin
        self.U = hamiltonian.U
        self.dtau = dtau
        self.BK_half = expm(-0.5 * dtau * self.K)
        self.gamma = np.arccosh(np.exp(0.5 * dtau * self.U))
        self.prefector = 0.5 * np.exp(-0.5 * self.dtau * self.U)


    # def propagate_walker_exact(self, walker):
    #     n_sites = walker.state.phi_up.shape[0]
    #     walker.state.phi_up = self.BK_half @ walker.state.phi_up
    #     walker.state.phi_down = self.BK_half @ walker.state.phi_down
    #     # print('-'*80)
    #     # print(walker.get_phi_up())
    #     # print(walker.get_phi_down())

    #     # print(n_sites)
        
    #     for i in range(n_sites):
    #         walker.state.phi_up[i] = self.prefector * (np.exp(self.gamma) + np.exp(-1 * self.gamma)) * walker.state.phi_up[i]
    #         walker.state.phi_down[i] = self.prefector * (np.exp(-1 * self.gamma) + np.exp(self.gamma)) * walker.state.phi_down[i]        
    #     # print('-'*80)
    #     # print(walker.get_phi_up())
    #     # print(walker.get_phi_down())

    #     # print('-'*80)
    #     # print(self.gamma)
    #     # print(self.prefector)
    #     # print((np.exp(self.gamma) + np.exp(-1 * self.gamma)))
    #     # print(self.prefector * (np.exp(self.gamma) + np.exp(-1 * self.gamma)))

    #     walker.state.phi_up = self.BK_half @ walker.state.phi_up
    #     walker.state.phi_down = self.BK_half @ walker.state.phi_down


    # def exact_propagator(self)


    # def sample_fields(self, nsites):
    #     """
    #     x_i = ±1
    #     """
    #     return np.random.choice([-1, 1], size=nsites)
    

    
    # def exact_interaction_propagator(self, x):
    #     prefactor = np.exp(-0.5 * self.dtau * self.U)
    #     Bup = np.diag(prefactor * np.exp(+self.gamma * x))
    #     Bdn = np.diag(prefactor * np.exp(-self.gamma * x))
    #     # Bup = np.diag(np.exp(+self.gamma * x))
    #     # Bdn = np.diag(np.exp(-self.gamma * x))

    #     Bup_x = 0
    #     Bdn_x = 0

    #     return Bup, Bdn
    

    # def exact_propagate(self, walker, x = None):
    #     if x.all() == None:
    #         x = self.sample_fields(walker.get_phi_up().shape[0])
    #     Bup, Bdn = self.exact_interaction_propagator(x)

    #     phi_up_new = (self.BK_half @ Bup @ self.BK_half @ walker.get_phi_up())
    #     walker.state.phi_up = phi_up_new
    #     phi_down_new = (self.BK_half @ Bdn @ self.BK_half @ walker.get_phi_down())
    #     walker.state.phi_down = phi_down_new

    #     return (phi_up_new, phi_down_new, x)

