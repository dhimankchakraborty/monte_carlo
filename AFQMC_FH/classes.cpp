#include "classes.hpp"
#include <iostream>

SlaterWaveFunction::SlaterWaveFunction(int n_basis, int n_particle, int n_species)
    : n_basis(n_basis), n_particle(n_particle), n_species(n_species),
      wave_function(n_species, n_basis, n_particle) {
    wave_function.setZero();
}

double& SlaterWaveFunction::operator()(int s, int i, int j) {
    return wave_function(s, i, j);
}

const double& SlaterWaveFunction::operator()(int s, int i, int j) const {
    return wave_function(s, i, j);
}

void SlaterWaveFunction::print_dimensions() const {
    std::cout << "Species: " << wave_function.dimension(0) 
              << ", Basis: " << wave_function.dimension(1) 
              << ", Particles: " << wave_function.dimension(2) << std::endl;
}