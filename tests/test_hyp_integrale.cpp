#include <gtest/gtest.h>
#include <voltlbx/hyp_integrale.h>
#include <voltlbx/utils.h>

using namespace voltlbx;

void check_hypint_derivative(const HyperbolicIntegral hyp_int, std::pair<double, double> domain)
{
    auto[inf, sup] = domain;
    ASSERT_TRUE(inf < sup);
    auto xs = linspace(inf, sup, 101);
    
    auto fd_deriv = [&](double x)
    {
        constexpr double EPS = 1.0e-6;
        return (hyp_int(x + EPS) - hyp_int(x - EPS)) / (2.0 * EPS);
    };

    for (auto x : xs)
    {
        const double deriv = fd_deriv(x);
        const double ref_deriv = 1.0 / std::sqrt(hyp_int.a * x * x + hyp_int.b * x + hyp_int.c);
        ASSERT_NEAR(deriv, ref_deriv, 1.0e-9);
    }
}

void check_hypint(const HyperbolicIntegral hyp_int)
{
    auto[s, e] = hyp_int.domain();
    ASSERT_LT(0.0, e);
    ASSERT_LT(s, 0.0);

    ASSERT_DOUBLE_EQ(hyp_int(0.0), 0.0);

    constexpr double EPS_BOUND = 0.1;
    check_hypint_derivative(hyp_int, { std::max(-10.0, s + EPS_BOUND), std::min(10.0, e - EPS_BOUND) });
}

TEST(HyperbolicIntegral, CheckDerivative)
{
    check_hypint(HyperbolicIntegral(2.50, 2.10, 0.51));
    check_hypint(HyperbolicIntegral(1.90, 2.10, 0.5));
    check_hypint(HyperbolicIntegral(1.90, -2.10, 0.5));
}


void check_value_solve(double b, double c, double x)
{
    const double max_value = HyperbolicIntegral::max_value(b, c, x);
    ASSERT_GT(x * max_value, 0.0);

    auto vs = x > 0.0 ? linspace(max_value * 0.01, max_value * 0.99, 101)
                      : linspace(max_value * 0.99, max_value * 0.01, 101);

    for (auto value : vs)
    {
        auto hyp_int = HyperbolicIntegral::solve_a_from_value(b, c, x, value);
        const double solved_value = hyp_int(x);
        ASSERT_NEAR(solved_value, value, 1.0e-10);
    }
}

TEST(HyperbolicIntegral, CheckValueSolve)
{
    check_value_solve(-1.0, 1.0, 0.5);
    check_value_solve(-1.0, 1.0, 1.0);
    check_value_solve(-1.0, 1.0, -1.0);
    check_value_solve(-1.0, 1.0, -0.5);
}

