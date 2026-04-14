#ifndef CLASSES_HPP
#define CLASSES_HPP

#include <unsupported/Eigen/CXX11/Tensor>

class SlaterWaveFunction {
private:
    int n_basis;
    int n_particle;
    int n_species;
    Eigen::Tensor<double, 3> wave_function;

public:
    SlaterWaveFunction(int n_basis, int n_particle, int n_species);

    double& operator()(int s, int i, int j);
    const double& operator()(int s, int i, int j) const;
    
    void print_dimensions() const;
};

#endif