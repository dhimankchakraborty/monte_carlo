// Compilation: g++ -O3 -march=native -ffast-math -std=c++17 test-t02.cpp functions.cpp -o test-t02
// Run: ./test-t02

#include <iostream>
#include <random>
#include <cmath>
#include <iomanip>
#include <cstdint>
#include <vector>
#include <stdexcept>
#include <fstream>


#include "functions.hpp"




int main()
{
    double x_init = 0.0;
    double step_size = 0.5;
    uint64_t n_mc_steps = 100'000'000;
    uint64_t n_therm_steps = n_mc_steps / 10;

    double beta_start = 0.0;
    double beta_end = 1.0;
    size_t run_count = 101;
    const string filename = "result-01.csv";

    vector<double> beta_arr = linspace(beta_start, beta_end, run_count);

    double x_current, x_next, acceptance_ratio, total_energy, total_energy_squared, local_energy, beta;
    vector<double> energy_arr(run_count), varience_arr(run_count), standard_deviation_arr(run_count);

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(0.0, 1.0);

    random_device rd_step;
    mt19937 gen_step(rd_step());
    uniform_real_distribution<double> dist_step(-1*step_size, step_size);

    random_device rd_acceptance;
    mt19937 gen_acceptance(rd_acceptance());
    uniform_real_distribution<double> dist_acceptance(0.0, 1.0);


    for (size_t i = 0; i < run_count; i++)
    {
        beta = beta_arr[i];

        if (x_init == 0.0)
        {
            x_current = dist(gen);
        }
        else
        {
            x_current = x_init;
        }


        for (size_t i = 0; i < n_therm_steps; i++)
        {
            x_next = x_current + dist_step(gen);
            acceptance_ratio = acceptance_ratio_t01(beta, x_current, x_next);

            if (acceptance_ratio > dist_acceptance(gen))
            {
                x_current = x_next;
            }
        }

        total_energy = 0;
        total_energy_squared = 0;
        for (size_t i = 0; i < n_mc_steps; i++)
        {
            x_next = x_current + dist_step(gen);
            acceptance_ratio = acceptance_ratio_t01(beta, x_current, x_next);

            if (acceptance_ratio > dist_acceptance(gen))
            {
                x_current = x_next;
            }

            local_energy = local_energy_t01(beta, x_current);
            total_energy += local_energy;
            total_energy_squared += local_energy*local_energy;
        }

        double energy = total_energy / n_mc_steps;
        double varience = (total_energy_squared / n_mc_steps) - (energy*energy);
        double standard_deviation = sqrt(varience);

        cout << "Beta                  :" << beta << endl;
        cout << "Energy                :" << energy << endl;
        cout << "Varience              :" << varience << endl;
        cout << "Standard Deviation    :" << standard_deviation << endl;
        cout << "----------------------------------------------------------------------" << endl;

        energy_arr[i] = energy;
        varience_arr[i] = varience;
        standard_deviation_arr[i] = standard_deviation;
    }

    ofstream file(filename);
    file << "beta,energy,varience,standard_deviation";
    for (std::size_t i = 0; i < run_count; ++i) {
        file << "\n";
        file << beta_arr[i] << "," << energy_arr[i] << "," << varience_arr[i] << "," << standard_deviation_arr[i];
    }

    file.close();
    cout << "Result written to " << filename << "\n";
    
    return 0;
}