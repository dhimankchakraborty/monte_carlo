import numpy as np




class slaterWaveFunction:
    def __init__(self, single_particle_basis_number, particle_number) -> None:
        self.wave_function = np.zeros((single_particle_basis_number, particle_number))


wf = slaterWaveFunction(3, 2)