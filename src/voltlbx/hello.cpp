#include "hello.h"

#include <Eigen/Dense>

namespace voltlbx
{
    
    int add(int i, int j)
    {        
        return i + j;
    }

    using Matrix = Eigen::MatrixXd;
    using Vector = Eigen::VectorXd;
    void test_eigen()
    {
        Matrix m1 = Matrix::Zero(2, 2);
        m1(0, 0) = 3.0;
        m1(1, 1) = 50.0;

        Vector v(2);
        v(0) = 1.0;
        v(1) = 1.0;

        auto r = m1 * v;       
    }

}
