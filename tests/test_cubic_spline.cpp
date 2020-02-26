#include <gtest/gtest.h>
#include <voltlbx/cubic_spline.h>

using namespace voltlbx;

void check_natural_spline_continuity(const std::vector<double>& xs,
                                     const std::vector<double>& ys)
{
    ASSERT_EQ(xs.size(), ys.size());
    CubicSpline s(xs, ys);
    for (long i = 0; i < xs.size(); ++i)
    {
        double x = xs[i];
        ASSERT_DOUBLE_EQ(s(x), ys[i]);
        
        auto[x_l, x_p] = [x]()
        {
            if (std::abs(x) > 0.0)
            {
                return std::tuple { x * (1.0 - std::numeric_limits<double>::epsilon()),
                                    x * (1.0 + std::numeric_limits<double>::epsilon()) };
            }                
            else
            {
                return std::tuple { -std::numeric_limits<double>::min(), 
                                    std::numeric_limits<double>::min() };
            }
        }();

        ASSERT_DOUBLE_EQ(s(x_l), ys[i]);
        ASSERT_DOUBLE_EQ(s(x_p), ys[i]);
    }
}

TEST(CubicSpline, CheckContinuity)
{
    check_natural_spline_continuity({ -2.0, -1.0, 0.0, 1.0, 2.0 },
                                    { 1.6, 1.1, 1.0, 0.9, 0.7 });

    check_natural_spline_continuity({ -1.0, 1.0 },
                                    { 1.0, 0.5 });

    check_natural_spline_continuity({ 0.0, 0.5, 0.99, 2.5 },
                                    { 5.0, 5.4, 3.1, 1.0 });
}