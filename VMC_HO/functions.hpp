#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <iostream>
#include <random>
#include <cmath>
#include <iomanip>
#include <cstdint>

using namespace std;


double psi_squared(double x, double beta);

double local_energy(double x, double beta);

double metropolis_step(double x_current, double step_size, double beta,  uniform_real_distribution<double>& uniform, mt19937& rng);

#endif