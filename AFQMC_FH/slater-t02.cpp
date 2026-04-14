// Compilation: g++ -O3 -I /home/dhiman/packages/eigen-5.0.0 classes.cpp slater-t02.cpp -o slater-t02
// Run: ./slater-t02

#include <iostream>
#include "classes.hpp"




int main() {
    int s = 2, b = 4, p = 3;
    SlaterWaveFunction swf(b, p, s);

    
    swf.print_dimensions();
    for (size_t i = 0; i < s; i++)
    {
        for (size_t j = 0; j < b; j++)
        {
            for (size_t k = 0; k < p; k++)
            {
                std::cout << "swf(" << i << ',' << j << ',' << k << "):" << swf(i, j, k) << std::endl;
            }
            
        }
        
    }

    return 0;
}