#include <gtest/gtest.h>
#include <voltlbx/utils.h>

using namespace voltlbx;

void check_locate_index(const std::vector<double>& xs)
{
    for (long i = 0; i < xs.size(); ++i)
    {
        double x = xs[i];
        if (std::abs(x) > 0.0)
        {
            const double x_l = x * (1.0 - std::numeric_limits<double>::epsilon());
            ASSERT_EQ(locate_left_index(xs, x_l), i - 1);

            const double x_p = x * (1.0 + std::numeric_limits<double>::epsilon());
            ASSERT_EQ(locate_left_index(xs, x_p), i);
        }
        else
        {
            const double x_l = -std::numeric_limits<double>::min();
            ASSERT_EQ(locate_left_index(xs, x_l), i - 1);

            const double x_p = std::numeric_limits<double>::min();
            ASSERT_EQ(locate_left_index(xs, x_p), i);
        }
    }
}

TEST(Utils, CheckLocate)
{
    check_locate_index({ 0.0 });
    check_locate_index({ 3.14 });
    check_locate_index({ 0.0, 1.0, 2.0 });
}

TEST(Utils, LongestIncreasingSeq)
{
    std::vector<double> a = { 10, 22, 9, 33, 21, 50, 41, 60 };
    auto seq = longest_increasing_subsequence(a);
    std::vector<std::size_t> ref_seq = { 0, 1, 3, 6, 7 };
    ASSERT_EQ(seq.size(), ref_seq.size());
    for (std::size_t i = 0; i < seq.size(); ++i)
    {
        ASSERT_EQ(seq[i], ref_seq[i]);
    }
}