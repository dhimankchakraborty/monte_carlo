import numpy as np




class Lattice:
    def __init__(self, dim):
        self.dim = dim



class Chain(Lattice):
    def __init__(self, n_sites: int, pbc: bool = True):
        self.n_sites = n_sites
        self.pbc = pbc
        self.dim = 1

        self.idx_dict = self.create_idx_dict()


    def create_idx_dict(self):
        return {(i,): i for i in range(self.n_sites)}
    


class Square(Lattice):
    def __init__(self, Lx: int, Ly: int, pbc: bool = True):

        self.Lx = Lx
        self.Ly = Ly
        self.n_sites = self.Lx * self.Ly
        self.pbc = pbc
        self.dim = 2

        self.idx_dict = self.create_idx_dict()
    

    def create_idx_dict(self):
        idx_dict = {}
        idx = 0
        for i in range(self.Lx):
            for j in range(self.Ly):
                idx_dict[(i, j)] = idx
                idx += 1
        return idx_dict