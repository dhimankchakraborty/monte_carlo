// Compilation: g++ -O3 -march=native -ffast-math -std=c++17 test-t01.cpp functions.cpp -o test-t01
// Run: ./test-t01
#include <iostream>
#include <random>
#include <cmath>
#include <iomanip>
#include <cstdint>

#include "functions.hpp"



// const int N_TOTAL = 1e6;
// const int N_THERM = 1e4;
// const double STEP_SIZE = 1.0;
// const double BETA = 0.5;






int main(int argc, char* argv[])
{
    if (argc < 5) 
    {
        std::cerr << "Usage: " << argv[0] << " <total_number_of_steps> <thermalization_steps> <step_size> <beta>\n";
        return 1;
    }

    uint64_t n_total = 0;
    uint64_t n_therm = 0;
    double step_size = 0.0;
    double beta = 0.0;

    try 
    {
        n_total = stoull(argv[1]);
        n_therm = stoull(argv[2]);
        step_size = stod(argv[3]);
        beta = stod(argv[4]);
    } 
    catch (...) 
    {
        cerr << "Error: Incorrect input datatype.\n";
        return 1;
    }

    double x = 0.0;

    double energy_sum = 0.0;
    double energy_sq_sum = 0.0;

    int accepted = 0;

    random_device rd;
    mt19937 rng(rd());
    uniform_real_distribution<double> uniform(0.0, 1.0);

    for (int i = 0; i < n_total; ++i) {

        double x_new = metropolis_step(x, step_size, beta, uniform, rng);

        if (x_new != x) accepted++;

        x = x_new;

        if (i >= n_therm) {

            double e = local_energy(x, beta);

            energy_sum += e;
            energy_sq_sum += e * e;
        }
    }

    int n_samples = n_total - n_therm;

    double E_avg = energy_sum / n_samples;
    double E2_avg = energy_sq_sum / n_samples;

    double variance = E2_avg - E_avg * E_avg;
    double std_error = sqrt(variance / n_samples);

    double acceptance_ratio = (double)accepted / n_total;

    // ------------------------------
    // Output results
    // ------------------------------
    cout << std::fixed << std::setprecision(6);

    cout << "====================================\n";
    cout << "VMC Results (1D Harmonic Oscillator)\n";
    cout << "====================================\n";

    cout << "Beta                = " << beta << "\n";
    cout << "Energy              = " << E_avg << "\n";
    cout << "Variance            = " << variance << "\n";
    cout << "Std Error           = " << std_error << "\n";
    cout << "Acceptance Ratio    = " << acceptance_ratio << "\n";

    cout << "====================================\n";

    return 0;
}