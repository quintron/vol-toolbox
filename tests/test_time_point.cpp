#include <gtest/gtest.h>

#include "voltlbx/time_utils.h"

using namespace voltlbx;

TEST(TimePoint, Create)
{
    auto check_datetime = [](unsigned year,
                             unsigned month,
                             unsigned day,
                             unsigned hour,
                             unsigned minute,
                             unsigned sec,
                             int ref_time_since_epoch)
    {
        auto t = utc_datetime(year, month, day, hour, minute, sec);
        ASSERT_EQ(t.time_since_epoch().count(), ref_time_since_epoch);

        auto [y, m, d] = get_year_month_day(t);
        ASSERT_EQ(y, year);
        ASSERT_EQ(m, month);
        ASSERT_EQ(d, day);

        auto [h, min, s] = get_hour_minute_second(t);
        ASSERT_EQ(h, hour);
        ASSERT_EQ(min, minute);
        ASSERT_EQ(s, sec);
    };

    check_datetime(2020, 12, 31, 15, 0, 0, 1609426800);
    check_datetime(2020, 2, 29, 15, 45, 10, 1582991110);
    check_datetime(2022, 12, 31, 23, 59, 59, 1672531199);

}
