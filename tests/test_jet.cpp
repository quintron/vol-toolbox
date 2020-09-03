#include <gtest/gtest.h>
#include <voltlbx/jet.h>

using namespace voltlbx;


void check_jet_entire_pow(const Jet& j, int n)
{
    assert(n > 0);
    auto pow_j = pow(j, 1.0 / n);
    auto j_tilde = pow_j;
    for (std::size_t i = 0; i + 1 < n; ++i)
        j_tilde = j_tilde * pow_j;

    ASSERT_DOUBLE_EQ(j.y, j_tilde.y);
    ASSERT_DOUBLE_EQ(j.dy_dx, j_tilde.dy_dx);
    ASSERT_DOUBLE_EQ(j.d2y_d2x, j_tilde.d2y_d2x);
}


TEST(Jet, Pow)
{
    check_jet_entire_pow(Jet{ 0.5, 1.5, 0.5 }, 1);
    check_jet_entire_pow(Jet{ 0.5, 1.5, 0.5 }, 2);
    check_jet_entire_pow(Jet{ 0.5, 1.5, 0.5 }, 3);
}