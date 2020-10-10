#include <gtest/gtest.h>
#include <voltlbx/roots.h>
#include <cmath>

using namespace voltlbx;

void check_bracket(const std::function<double(double)>& f, double x1, double x2)
{
    auto[xa, xb] = bracket_root(f, x1, x2);
    ASSERT_TRUE(f(xa) * f(xb) <= 0.0);
}

TEST(Bracket, Check)
{
    auto f = [](double x) {return  x * x - 9.0; };
    check_bracket(f, 0.0, 1.0);
    check_bracket(f, 1.0, 0.0);
    check_bracket(f, 5.0, 4.0);
    check_bracket(f, 4.0, 5.0);
    check_bracket(f, 2.0, 4.0);
    check_bracket(f, 3.0, 4.0);

    auto f2 = [](double x) {return  9.0 - x * x; };
    check_bracket(f2, 0.0, 1.0);
    check_bracket(f2, 1.0, 0.0);
    check_bracket(f2, 5.0, 4.0);
    check_bracket(f2, 4.0, 5.0);
    check_bracket(f2, 2.0, 4.0);
    check_bracket(f2, 3.0, 4.0);
}


TEST(Brent, Check)
{
    auto f = [](double x)
    {
        return std::sqrt(1.0 + x * x) - 3.0;
    };
    
    const auto[x1, x2] = bracket_root(f, 0.0, 1.0);
    const double root = brenth(f, x1, x2, 1.0e-8);       
    const double ref_root = std::sqrt(8.0);
    
    ASSERT_TRUE(std::abs(root - ref_root) < 1e-8);

}