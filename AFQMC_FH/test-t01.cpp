// Compilation: g++ -O3 -std=c++17 -I /home/dhiman/packages/eigen-5.0.0 functions.cpp classes.cpp test-t01.cpp -o test-t01
// Run: ./test-t01

#include <iostream>
#include "classes.hpp"
#include "functions.hpp"


using namespace std;




int main() {

    int Lx = 4;
    int Ly = 1;
    double t = 1.0;

    auto result = non_interacting_fermi_hubbard_ed(Lx, Ly, t, false);

    cout << "Hamiltonian:\n" << result.H << endl << endl;

    cout << "Eigenvalues:\n" << result.eigenvalues << endl << endl;

    cout << "Eigenvectors-1:\n" << result.eigenvectors.col(0) << endl;

    cout << "Eigenvectors-2:\n" << result.eigenvectors.col(1) << endl;

    cout << "Eigenvectors-3:\n" << result.eigenvectors.col(2) << endl;

    cout << "Eigenvectors-4:\n" << result.eigenvectors.col(3) << endl;

    return 0;
}