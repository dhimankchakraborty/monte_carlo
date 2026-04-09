// Compilation: g++ -O3 -march=native -fopenmp -std=c++17 test-t03.cpp functions.cpp -o test-t03
// Run: ./test-t03

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "functions.hpp"

using namespace std;

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        cerr << "Usage: " << argv[0] << " <max_points>\n";
        return 1;
    }

    uint64_t max_points = stoull(argv[1]);

    // Define step sizes (log-scale is better for convergence)
    vector<uint64_t> steps;

    for (uint64_t n = 1e3; n <= max_points; n *= 2)
    {
        steps.push_back(n);
    }
    steps.push_back(max_points);

    // Open CSV file
    ofstream file("results.csv");

    // CSV header
    file << "Points,Estimated_Pi,Error(%)\n";

    cout << setprecision(18);

    // ================================
    // Loop over step sizes (SEQUENTIAL)
    // ================================
    for (auto N : steps)
    {
        long double pi_est = monte_carlo_pi(N);

        long double error =
            ((pi_est - M_PI) / M_PI) * 100.0;

        // Print to terminal
        cout << "N = " << N
             << " | Pi = " << pi_est
             << " | Error (%) = " << fixed << setprecision(5)
             << error << endl;

        // Write to CSV
        file << N << ","
             << pi_est << ","
             << error << "\n";
    }

    file.close();

    cout << "\nResults saved to results.csv\n";

    return 0;
}