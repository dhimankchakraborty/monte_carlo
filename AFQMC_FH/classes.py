import numpy as np
from scipy.linalg import expm




class slaterWaveFunction2Spin:
    def __init__(self, n_basis, n_particle_up, n_particle_down) -> None:
        self.n_basis = n_basis
        self.n_particle_up = n_particle_up
        self.n_particle_down = n_particle_down

        self.wave_function_up = np.zeros((self.n_basis, self.n_particle_up))
        self.wave_function_down = np.zeros((self.n_basis, self.n_particle_down))

    def wave_function_initialization(self, up_sates, down_states):
        self.up_initializer = up_sates[:self.n_particle_up]
        self.down_initializer = down_states[:self.n_particle_down]

        for i in range(self.n_particle_up):
            self.wave_function_up[:, i] = self.up_initializer[i]
        
        for i in range(self.n_particle_down):
            self.wave_function_down[:, i] = self.down_initializer[i]



class imaginaryTimePropagator:
    def __init__(self, one_body_operator_matrix, delta_tau) -> None:
        self.one_body_operator_matrix = one_body_operator_matrix
        self.delta_tau = delta_tau
        self.propagator = expm(-self.delta_tau * self.one_body_operator_matrix)
        # self.propagator = np.exp(-self.delta_tau * self.one_body_operator_matrix) # Used exmp in code, if any problem arise, use np.exp to see if the result improves
    
    def apply_propagator(self, wave_function):
        self.n_basis = wave_function.n_basis
        self.n_particle_up = wave_function.n_particle_up
        self.n_particle_down = wave_function.n_particle_down

        modified_wave_function = slaterWaveFunction2Spin(self.n_basis, self.n_particle_up, self.n_particle_down)

        modified_wave_function.wave_function_up = self.propagator @ wave_function.wave_function_up

        modified_wave_function.wave_function_down = self.propagator @ wave_function.wave_function_down

        return modified_wave_function


