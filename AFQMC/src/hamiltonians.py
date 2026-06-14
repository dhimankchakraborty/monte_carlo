import numpy as np




class System:
    def __init__(self, lattice):
        """
        Base class for physical systems.
        """
        self.lattice = lattice
        self.n_basis = lattice.n_sites
        self.h_kin = None
        self.evals = None
        self.evecs = None


    def _build_non_interacting_hamiltonian(self):
        raise NotImplementedError("Subclasses must implement this method.")


    def get_non_interacting_evals_evecs(self, return_evals=True, return_evecs=True):
        """Diagonalizes the kinetic energy matrix and returns eigenvectors."""
        if self.h_kin is None:
            self.h_kin = self._build_non_interacting_hamiltonian()
            
        if self.evals is None or self.evecs is None:
            self.evals, self.evecs = np.linalg.eigh(self.h_kin)
        
        if return_evals and return_evecs:
            return self.evals, self.evecs
        elif return_evals:
            return self.evals
        elif return_evecs:
            return self.evecs
        else:
            return None
    


class HubbardSystem(System):
    def __init__(self, lattice, t: float = 1.0, U: float = 10.0):
        """
        Initializes the Hubbard model parameters and builds the non-interacting 
        Hamiltonian (hopping) matrix based on the provided lattice.
        """
        super().__init__(lattice=lattice)
        self.t = t
        self.U = U
        
        self.h_kin = self._build_non_interacting_hamiltonian()


    def _build_non_interacting_hamiltonian(self):
        """Builds the N x N non-interacting kinetic energy matrix."""
        H0 = np.zeros((self.n_basis, self.n_basis))
        
        if self.lattice.dim == 1:
            for i in range(self.lattice.n_sites):
                idx1 = self.lattice.idx_dict[(i,)]
                
                if i + 1 < self.lattice.n_sites:
                    idx2 = self.lattice.idx_dict[(i + 1,)]
                    H0[idx1, idx2] = -self.t
                    H0[idx2, idx1] = -self.t
                elif self.lattice.pbc and self.lattice.n_sites > 2:
                    idx2 = self.lattice.idx_dict[(0,)]
                    H0[idx1, idx2] = -self.t
                    H0[idx2, idx1] = -self.t
                    
        elif self.lattice.dim == 2:
            for x in range(self.lattice.Lx):
                for y in range(self.lattice.Ly):
                    idx1 = self.lattice.idx_dict[(x, y)]
                    
                    nx = (x + 1) % self.lattice.Lx
                    if nx > x or (self.lattice.pbc and self.lattice.Lx > 2):
                        idx2 = self.lattice.idx_dict[(nx, y)]
                        H0[idx1, idx2] = -self.t
                        H0[idx2, idx1] = -self.t
                        
                    ny = (y + 1) % self.lattice.Ly
                    if ny > y or (self.lattice.pbc and self.lattice.Ly > 2):
                        idx3 = self.lattice.idx_dict[(x, ny)]
                        H0[idx1, idx3] = -self.t
                        H0[idx3, idx1] = -self.t
                        
        return H0
    

    def calculate_variational_energy(self, walker):

        K = self.h_kin
        U = self.U

        phi_up = walker.state.phi_up
        phi_dn = walker.state.phi_down

        G_up = phi_up @ phi_up.conj().T
        G_dn = phi_dn @ phi_dn.conj().T

        E_kin = np.trace(K @ G_up).real
        E_kin += np.trace(K @ G_dn).real

        n_up = np.diag(G_up).real # type: ignore
        n_dn = np.diag(G_dn).real # type: ignore

        E_int = U * np.sum(n_up * n_dn)

        return E_kin + E_int
    

    def energy_mixed_estimator(self, trial_state, walker):
        phi_up = walker.state.phi_up
        phi_dn = walker.state.phi_down
        phi_T_up = trial_state.phi_up
        phi_T_dn = trial_state.phi_down

        K = self.h_kin
        U = self.U

        G_up = phi_up @ phi_T_up.conj().T
        G_dn = phi_dn @ phi_T_dn.conj().T

        E_kin = np.trace(K @ G_up).real
        E_kin += np.trace(K @ G_dn).real

        n_up = np.diag(G_up).real # type: ignore
        n_dn = np.diag(G_dn).real # type: ignore

        E_int = U * np.sum(n_up * n_dn)

        return E_kin + E_int
    