#include "functions.hpp"

#include <random>
#include <omp.h>

using namespace std;


long double monte_carlo_pi(std::uint64_t total_points_count)
{
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

    return 4.0L * points_inside_circle_count / total_points_count;
}