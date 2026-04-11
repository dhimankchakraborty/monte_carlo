#include "functions.hpp"

#include <random>
#include <omp.h>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;




double trial_wave_function_t01(double beta, double x)
{
    return exp(-1 * beta * x * x);
}


double local_energy_t01(double beta, double x)
{
    return beta + (0.5 - 2 * beta * beta) * x * x;
}


double probability_density_t01(double beta, double x)
{
    return (trial_wave_function_t01(beta, x) * trial_wave_function_t01(beta, x));
}


double acceptance_ratio_t01(double beta, double x_current, double x_next)
{
    return (trial_wave_function_t01(beta, x_next) * trial_wave_function_t01(beta, x_next) / (trial_wave_function_t01(beta, x_current) * trial_wave_function_t01(beta, x_current)));
}


// double metropolis_step(double x_current, double step_size, double beta, uniform_real_distribution<double> uniform, mt19937 rng)
// {
//     double delta = step_size * (uniform(rng) - 0.5);
//     double x_trial = x_current + delta;


//     double ratio = psi_squared(x_trial, beta) / psi_squared(x_current, beta);


//     if (ratio >= 1.0)
//     {
//         return x_trial;
//     }
//     else
//     {
//         double r = uniform(rng);
//         if (r < ratio)
//         {
//             return x_trial;
//         }
//         else
//         {
//             return x_current;
//         }
//     }
// }