import numpy as np
from scipy.linalg import expm
# from ./walkers import *




class HubbardPropagator:

    def __init__(self, K, U, dtau):

        self.K = K
        self.U = U
        self.dtau = dtau
        self.BK_half = expm(-0.5 * dtau * K)
        self.gamma = np.arccosh(np.exp(0.5 * dtau * U))


    def sample_fields(self, nsites):
        """
        x_i = ±1
        """
        return np.random.choice([-1, 1], size=nsites)
    

    
    def interaction_propagator(self, x):

        Bup = np.diag(np.exp(+self.gamma * x))
        Bdn = np.diag(np.exp(-self.gamma * x))

        return Bup, Bdn
    

    def propagate(self, walker):

        x = self.sample_fields(walker.get_phi_up().shape[0])

        Bup, Bdn = self.interaction_propagator(x)

        phi_up_new = (self.BK_half @ Bup @ self.BK_half @ walker.get_phi_up())
        walker.state.phi_up = phi_up_new

        phi_down_new = (self.BK_half @ Bdn @ self.BK_half @ walker.get_phi_down())
        walker.state.phi_down = phi_down_new

        return (phi_up_new, phi_down_new, x)

