#include <iostream>
#include <random>
#include <cstdint>
#include <cmath>
#include <iomanip>
#include <chrono>


using namespace std;



int main(int argc, char* argv[])
{

    if (argc < 2) 
    {
        std::cerr << "Usage: " << argv[0] << " <number_of_points>\n";
        return 1;
    }

    uint64_t total_points_count = 0;
    try 
    {
        total_points_count = stoull(argv[1]);
    } 
    catch (...) 
    {
        cerr << "Error: Input must be a positive integer.\n";
        return 1;
    }


    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> dist(0.0, 1.0);


    uint64_t points_inside_circle_count = 0;

    auto start = chrono::high_resolution_clock::now();

    for (uint64_t i = 0; i < total_points_count; ++i)
    {
        double x = dist(gen);
        double y = dist(gen);

        if (x*x + y*y <= 1.0)
        {
            points_inside_circle_count++;
        }
    }

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;


    long double estimated_pi = 4.0L * points_inside_circle_count / total_points_count;


    cout << setprecision(18);

    cout << "Estimated value of Pi: " << estimated_pi << "\n";

    long double error_percent =
        ((estimated_pi - M_PI) / M_PI) * 100.0;

    cout << "Error in estimated value of Pi: "
              << fixed << setprecision(5)
              << error_percent << " %\n";

    cout << "Execution time: " << elapsed.count() << " seconds\n";

    return 0;
}