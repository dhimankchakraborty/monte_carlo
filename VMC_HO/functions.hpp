#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <iostream>
#include <random>
#include <cmath>
#include <iomanip>
#include <cstdint>
#include <vector>
#include <stdexcept>

using namespace std;


double trial_wave_function_t01(double beta, double x);

double local_energy_t01(double beta, double x);

double probability_density_t01(double beta, double x);

double acceptance_ratio_t01(double beta, double x_current, double x_next);

vector<double> linspace(double start, double end, size_t num);

#endif