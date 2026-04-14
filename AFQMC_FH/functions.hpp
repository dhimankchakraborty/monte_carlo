#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <iostream>
#include <vector>
#include <Eigen/Dense>


using namespace Eigen;
using namespace std;




int site_index(int x, int y, int length_x);

struct HubbardResult {
    MatrixXd H;
    VectorXd eigenvalues;
    MatrixXd eigenvectors;
};

HubbardResult non_interacting_fermi_hubbard_ed(
    int length_x,
    int length_y,
    double t,
    bool use_pbc = false
);




#endif