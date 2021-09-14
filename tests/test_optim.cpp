#include <gtest/gtest.h>
#include <voltlbx/optim.h>

using namespace voltlbx;


TEST(OptimTest, Brent)

{
    auto f = [](double x) {
        return (x - 3) * (x - 3) + 2;
    };
    OptimBracket b;
    b.a = 0;
    b.b = 1;
    b.c = 9;
    b.fa = f(b.a);
    b.fb = f(b.b);
    b.fc = f(b.c);

    static const auto tol = std::sqrt(std::numeric_limits<double>::epsilon());
    double fmin;
    double xmin = optimize_brent(f, b, tol, &fmin);
    EXPECT_NEAR(xmin, 3, tol);
    EXPECT_NEAR(fmin, 2, tol * tol);
}


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