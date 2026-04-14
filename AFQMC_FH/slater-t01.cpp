#include <vector>
#include <iostream>


using namespace std;



class SlaterWaveFunction {
private:
    int n_basis;
    int n_particle;
    int n_species;
    vector<double> wave_function;

public:
    SlaterWaveFunction(int n_basis, int n_particle, int n_species)
        : n_basis(n_basis), n_particle(n_particle), n_species(n_species) {
        wave_function.resize(n_species * n_basis * n_particle, 0.0);
    }

    double& operator()(int s, int b, int p) {
        return wave_function[s * (n_basis * n_particle) + b * n_particle + p];
    }

    const double& operator()(int s, int b, int p) const {
        return wave_function[s * (n_basis * n_particle) + b * n_particle + p];
    }
};


int main() {
    int species = 2;   // e.g., Spin up and Spin down
    int basis = 4;   // Number of basis functions
    int particles = 3; // Number of electrons

    SlaterWaveFunction wf(basis, particles, species);

    for (size_t i = 0; i < species; i++)
    {
        for (size_t j = 0; j < basis; i++)
        {
            for (size_t k = 0; k < particles; i++)
            {
                cout << "wf(" << i << ',' << j << ',' << k << "):" << wf(i, j, k) << endl;
            }
            
        }
        
    }
    
    

    return 0;
}


// Error:  ./slater-t01 
// wf(0,0,0):0
// wf(1,0,0):0
// /usr/include/c++/15/bits/stl_vector.h:1263: std::vector<_Tp, _Alloc>::reference std::vector<_Tp, _Alloc>::operator[](size_type) [with _Tp = double; _Alloc = std::allocator<double>; reference = double&; size_type = long unsigned int]: Assertion '__n < this->size()' failed.
// Aborted (core dumped)