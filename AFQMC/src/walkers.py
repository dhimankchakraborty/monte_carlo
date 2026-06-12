import copy
import numpy as np


class Walker:
    """
    AFQMC walker.

    A walker contains:
        - A Slater determinant state
        - Monte Carlo weight
        - Overlap with trial state
        - Sign / phase information
    """

    def __init__(self, state, weight: complex = 1.0):

        self.state = state
        self.weight = np.complex128(weight)

        self.overlap = None
        self.local_energy = None
        self.phase = 1.0


    def copy(self):
        """
        Deep-copy the walker.
        Essential for branching/population control.
        """
        return copy.deepcopy(self)


    def compute_overlap(self, trial_state):
        self.overlap = trial_state.overlap_calculation_logdet(self.state)

        return self.overlap


    def orthogonalize(self):
        self.state.orthogonalize_columns()


    def get_phi_up(self):
        return self.state.phi_up


    def get_phi_down(self):
        return self.state.phi_down