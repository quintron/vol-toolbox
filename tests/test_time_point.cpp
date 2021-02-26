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

auto test_calendar()
{
    std::vector<chrono::Date> holidays;
    for (int i = 0; i < 50; ++i)
    {
        int year = 2020 + i;
        holidays.emplace_back(create_date(year, 1, 1));
        holidays.emplace_back(create_date(year, 7, 4));
        holidays.emplace_back(create_date(year, 12, 25));
    }
    return Calendar({ date::Saturday , date::Sunday }, holidays);
}


TEST(Calendar, IsClosed)
{
    auto cal = test_calendar();

    auto check_closed = [&](chrono::Date d, bool expected_is_closed)
    {
        auto is_closed = cal.is_closed(d);
        ASSERT_EQ(is_closed, expected_is_closed);
    };
    
    check_closed(create_date(2024, 1, 5), false); //friday
    check_closed(create_date(2024, 1, 6), true); //saturday
    check_closed(create_date(2024, 1, 7), true); //sunday
    check_closed(create_date(2024, 1, 8), false); //monday

    check_closed(create_date(2020, 1, 1), true);
    check_closed(create_date(2021, 1, 1), true);
    check_closed(create_date(2022, 1, 1), true);
    check_closed(create_date(2023, 1, 1), true);
}

TEST(Calendar, CountOpenDays)
{
    auto cal = test_calendar();

    auto check_open_days = [&](chrono::Date s, chrono::Date e, int expected_nb_opens)
    {
        auto nb_opens = cal.count_open_days(s, e);
        ASSERT_EQ(nb_opens, expected_nb_opens);
    };
    
    check_open_days(create_date(2021, 2, 11),
                    create_date(2021, 2, 16), 
                    3);

    check_open_days(create_date(2020, 12, 30),
                    create_date(2021, 1, 4),
                    2);
    
    check_open_days(create_date(2020, 12, 30),
                    create_date(2021, 1, 5),
                    3);

    check_open_days(create_date(2023, 12, 23),
                    create_date(2024, 1, 2),
                    4);

    check_open_days(create_date(2020, 2, 25),
                    create_date(2021, 2, 25),
                    260);

    check_open_days(create_date(2021, 1, 2),
                    create_date(2022, 1, 2),
                    260);

}
