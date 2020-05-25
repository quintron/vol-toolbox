#include <gtest/gtest.h>
#include <voltlbx/optim.h>

using namespace voltlbx;

TEST(MinPack, Optim)
{

    double a = 1.0;
    double b = 10.0;

    auto rosenbrock = [&](std::size_t n, const double *x, std::size_t m, double *f)
    {
        ASSERT_EQ(n, 2);
        ASSERT_EQ(m, 2);
        f[0] = a - x[0];
        f[1] = b * (x[1] - x[0] * x[0]);
    };

    std::vector<double> x = { 0.0, 0.0 };
    std::vector<double> f(2);
    double tol = 1.0e-7;

    levmarq(rosenbrock, x, f, tol);
    EXPECT_NEAR(x[0], 1.0, tol);
    EXPECT_NEAR(x[1], 1.0, tol);
    EXPECT_NEAR(f[0], 0.0, tol);
    EXPECT_NEAR(f[1], 0.0, tol);
}