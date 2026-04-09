#include "functions.hpp"

#include <random>
#include <omp.h>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;




double psi_squared(double x, double beta) 
{
    return exp(-2.0 * beta * x * x);
}


double local_energy(double x, double beta)
{
    return beta + (0.5 - 2.0 * beta * beta) * x * x;
}


double metropolis_step(double x_current, double step_size, double beta, uniform_real_distribution<double> uniform, mt19937 rng)
{
    double delta = step_size * (uniform(rng) - 0.5);
    double x_trial = x_current + delta;


    double ratio = psi_squared(x_trial, beta) / psi_squared(x_current, beta);


    if (ratio >= 1.0)
    {
        return x_trial;
    }
    else
    {
        double r = uniform(rng);
        if (r < ratio)
        {
            return x_trial;
        }
        else
        {
            return x_current;
        }
    }
}