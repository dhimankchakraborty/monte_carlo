#include <iostream>
#include <vector>
#include <Eigen/Dense>


using namespace Eigen;
using namespace std;




int site_index(int x, int y, int length_x)
{
    return x + y * length_x;
}


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
) 
{
    int H_dim = length_x * length_y;

    MatrixXd H = MatrixXd::Zero(H_dim, H_dim);

    for (int x = 0; x < length_x; ++x)
    {
        for (int y = 0; y < length_y; ++y)
        {

            int i = site_index(x, y, length_x);

            int x_plus = x + 1;

            if (x_plus < length_x)
            {
                int j = site_index(x_plus, y, length_x);
                H(i, j) = -t;
                H(j, i) = -t;
            }
            else if (use_pbc)
            {
                int j = site_index(0, y, length_x);
                H(i, j) = -t;
                H(j, i) = -t;
            }


            int y_plus = y + 1;

            if (y_plus < length_y)
            {
                int j = site_index(x, y_plus, length_x);
                H(i, j) = -t;
                H(j, i) = -t;
            }
            else if (use_pbc)
            {
                int j = site_index(x, 0, length_x);
                H(i, j) = -t;
                H(j, i) = -t;
            }
        }
    }


    for (int i = 0; i < H_dim; ++i)
    {
        H(i, i) = 0.0;
    }


    SelfAdjointEigenSolver<MatrixXd> solver(H);

    if (solver.info() != Success)
    {
        throw runtime_error("Eigenvalue decomposition failed!");
    }

    HubbardResult result;
    result.H = H;
    result.eigenvalues = solver.eigenvalues();
    result.eigenvectors = solver.eigenvectors(); // columns = eigenvectors

    return result;
}
