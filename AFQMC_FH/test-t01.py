import numpy as np
from classes import *
from functions import *




H, evals, evecs = non_interacting_fermi_hubbard_ed(4, 1, 1, use_pbc=False)

print(H)
print(evals)
print(evecs)