import numpy as np




class SlaterDeterminantTwoSpinState:
    def __init__(self, hamiltonian, n_electrons_up: int | None = None,  n_electrons_down: int | None = None, initial_phi_up: np.ndarray | None = None, initial_phi_down: np.ndarray | None = None, dtype=np.complex128):
        """
        Initializes the Slater determinant with N x M_up and N x M_down matrices.
        N is the size of the basis, M is the number of occupied orbitals (electrons).
        """
        self.dtype = dtype
        if initial_phi_up is not None and initial_phi_down is not None:
            self.phi_up = np.array(initial_phi_up, dtype=dtype)
            self.phi_down = np.array(initial_phi_down, dtype=dtype)
            self.hamiltonian = hamiltonian
            self.n_basis = hamiltonian.n_basis
            _, self.n_electrons_up = self.phi_up.shape
            _, self.n_electrons_down = self.phi_down.shape
        
        elif n_electrons_up is not None and n_electrons_down is not None:
            self.hamiltonian = hamiltonian
            self.n_basis = self.hamiltonian.n_basis
            self.n_electrons_up = n_electrons_up
            self.n_electrons_down = n_electrons_down
            self.phi_up = np.zeros((self.n_basis, n_electrons_up), dtype=dtype)
            self.phi_down = np.zeros((self.n_basis, n_electrons_down), dtype=dtype)
        
        else:
            raise ValueError("Must provide either initial_phi_up and initial_phi_down or n_basis, n_electrons_up and n_electrons_down.")


    def __repr__(self):
        return f"SlaterDeterminant(n_basis={self.n_basis}, n_electrons_up={self.n_electrons_up}, n_electrons_down={self.n_electrons_down})"


    def norm(self):
        """
        Calculates the squared norm <self | self> of the Slater determinant.
        """
        return self.overlap_calculation_logdet(self)


    def initialize(self, method: str):
        if method == "random":
            self.phi_up = np.random.rand(self.n_basis, self.n_electrons_up)
            self.phi_down = np.random.rand(self.n_basis, self.n_electrons_down)
            self.orthogonalize_columns()

        elif method == "zero":
            self.phi_up = np.zeros((self.n_basis, self.n_electrons_up), dtype=self.dtype)
            self.phi_down = np.zeros((self.n_basis, self.n_electrons_down), dtype=self.dtype)

        elif method == "non-interacting":
            evecs = self.hamiltonian.get_non_interacting_evals_evecs(return_evals=False, return_evecs=True)
            self.phi_up = evecs[:, :self.n_electrons_up]
            self.phi_down = evecs[:, :self.n_electrons_down]

        else:
            raise ValueError("Invalid initialization method.")
        
    
    def orthogonalize_columns(self, return_R=False):
        Q_up, R_up = np.linalg.qr(self.phi_up)
        Q_down, R_down = np.linalg.qr(self.phi_down)
        self.phi_up = Q_up.copy()
        self.phi_down = Q_down.copy()

        if return_R:
            return R_up, R_down
    

    def overlap_calculation_raw(self, other):
        """
        Computes the overlap:
            <self | other>
        between two Slater determinants.
        """

        overlap_up = self.phi_up.conj().T @ other.phi_up
        overlap_down = self.phi_down.conj().T @ other.phi_down

        det_up = np.linalg.det(overlap_up)
        det_down = np.linalg.det(overlap_down)

        return det_up * det_down
    

    def overlap_calculation_logdet(self, other, return_log=False):
        """
        Computes overlap:

            <self | other>

        using numerically stable slogdet.
        """

        overlap_up = self.phi_up.conj().T @ other.phi_up
        overlap_down = self.phi_down.conj().T @ other.phi_down

        sign_up, logdet_up = np.linalg.slogdet(overlap_up)
        sign_down, logdet_down = np.linalg.slogdet(overlap_down)

        total_sign = sign_up * sign_down
        total_logdet = logdet_up + logdet_down

        if return_log:
            return total_sign, total_logdet

        return total_sign * np.exp(total_logdet)
        


class SlaterDeterminantSpinlessState:
    def __init__(self, n_basis: int | None = None, n_electrons: int | None = None, initial_phi: np.ndarray | None = None, dtype=np.complex128):
        """
        Initializes the Slater determinant with N x M matrices.
        N is the size of the basis, M is the number of occupied orbitals (electrons).
        """
        self.dtype = dtype
        if initial_phi is not None:
            self.phi = np.array(initial_phi, dtype=dtype)
            self.n_basis, self.n_electrons = self.phi.shape
        
        elif n_basis is not None and n_electrons is not None:
            self.n_basis = n_basis
            self.n_electrons = n_electrons
            self.phi = np.zeros((n_basis, n_electrons), dtype=dtype)
        
        else:
            raise ValueError("Must provide either initial_phi_up and initial_phi_down or n_basis, n_electrons_up and n_electrons_down.")


    def __repr__(self):
        return f"SlaterDeterminant(n_basis={self.n_basis}, n_electrons={self.n_electrons})"


    def norm(self):
        """
        Calculates the squared norm <self | self> of the Slater determinant.
        """
        return self.overlap_calculation_logdet(self)


    def initialize(self, system, method: str):
        if method == "random":
            self.phi = np.random.rand(self.n_basis, self.n_electrons)

        elif method == "zero":
            self.phi = np.zeros((self.n_basis, self.n_electrons), dtype=self.dtype)

        elif method == "non-interacting":
            evecs = system.get_non_interacting_evals_evecs(return_evals=False, return_evecs=True)
            self.phi = evecs[:, :self.n_electrons]

        else:
            raise ValueError("Invalid initialization method.")
        

    def orthogonalize_columns(self, return_R=False):
        Q, R = np.linalg.qr(self.phi)
        self.phi = Q.copy()

        if return_R:
            return R
        

    def overlap_calculation_raw(self, other):
        """
        Computes the overlap:
            <self | other>
        between two Slater determinants.
        """

        overlap = self.phi.conj().T @ other.phi
        det = np.linalg.det(overlap)

        return det
    

    def overlap_calculation_logdet(self, other, return_log=False):
        """
        Computes overlap:

            <self | other>

        using numerically stable slogdet.
        """

        overlap = self.phi.conj().T @ other.phi
        sign, logdet = np.linalg.slogdet(overlap)

        if return_log:
            return sign, logdet

        return sign * np.exp(logdet)