#include "functions.hpp"

#include <random>
#include <omp.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <stdexcept>


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


vector<double> linspace(double start, double end, size_t num)
{
    if (num <= 0)
    {
        throw invalid_argument("Number of points must be > 0");
    }
    else if (num == 1)
    {
        return {start};
    }

    vector<double> result(num);
    double step = (end - start) / (num - 1);

    for (size_t i = 0; i < num; ++i)
    {
        result[i] = start + i * step;
    }

    return result;
}