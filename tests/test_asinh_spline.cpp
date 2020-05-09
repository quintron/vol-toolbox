#include <gtest/gtest.h>
#include <voltlbx/asinh_spline.cpp>

#include <voltlbx/utils.h>

using namespace voltlbx;


void check_spline_solving(double second_deriv0)
{
    const double max_v1 = max_value1_from_left_second_deriv(second_deriv0);
    auto v1s = linspace(max_v1 * 0.01, max_v1 * (1.0 - 1.0e-5), 50);
    for (auto v1 : v1s)
    {
        auto spline_elem = solve_from_left_second_deriv_and_right_value(second_deriv0, v1);
        ASSERT_NEAR(spline_elem(1.0), v1, 1.0e-10);
        
        const double EPS = 1.0e-6;
        const double deriv = (spline_elem(EPS) - spline_elem(-EPS)) / (2.0 * EPS);
        ASSERT_NEAR(deriv, 1.0, 1.0e-5);
    }
}


TEST(AsinhSpline, TestSolve)
{   
    check_spline_solving(3.90);
    check_spline_solving(1.0);
    check_spline_solving(-5.0);
    check_spline_solving(-1.0);
    check_spline_solving(-0.5);            

}
