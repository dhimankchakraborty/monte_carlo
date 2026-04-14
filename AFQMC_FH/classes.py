import numpy as np




class slaterWaveFunction2Spin:
    def __init__(self, n_basis, n_particle_up, n_particle_down) -> None:
        self.n_basis = n_basis
        self.n_particle_up = n_particle_up
        self.n_particle_down = n_particle_down

        self.wave_function_up = np.zeros((self.n_basis, self.n_particle_up))
        self.wave_function_down = np.zeros((self.n_basis, self.n_particle_down))


