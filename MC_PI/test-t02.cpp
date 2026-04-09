// Compilation: g++ -O3 -march=native -fopenmp -std=c++17 test-t02.cpp -o test-t02
// Run: .\test-02

#include <iostream>
#include <random>
#include <cstdint>
#include <cmath>
#include <iomanip>
#include <omp.h>


using namespace std;



int main(int argc, char* argv[])
{

    if (argc < 2) 
    {
        cerr << "Usage: " << argv[0] << " <number_of_points>\n";
        return 1;
    }

    uint64_t total_points_count = stoull(argv[1]);

    uint64_t points_inside_circle_count = 0;


    #pragma omp parallel
    {
        random_device rd;
        mt19937_64 gen(rd() + omp_get_thread_num());

        uniform_real_distribution<double> dist(0.0, 1.0);


        uint64_t points_inside_circle_count_local = 0;

        #pragma omp for
        for (uint64_t i = 0; i < total_points_count; ++i)
        {
            double x = dist(gen);
            double y = dist(gen);

            if (x*x + y*y <= 1.0)
            {
                points_inside_circle_count_local++;
            }
        }

        #pragma omp atomic
        points_inside_circle_count += points_inside_circle_count_local;
    }

    long double estimated_pi = 4.0L * points_inside_circle_count / total_points_count;

    cout << setprecision(18);
    cout << "Estimated Pi: " << estimated_pi << "\n";

    long double error =
        ((estimated_pi - M_PI) / M_PI) * 100.0;

    cout << "Error (%): "
        << fixed << setprecision(5)
        << error << "\n";

    return 0;
}