#pragma once
#include <chrono>
#include <tuple>

namespace voltlbx
{


    using time_point = std::chrono::time_point<std::chrono::system_clock, std::chrono::seconds>;


    time_point utc_datetime(unsigned year,
                            unsigned month,
                            unsigned day,
                            unsigned hour,
                            unsigned minute,
                            unsigned sec);


    std::tuple<int, unsigned, unsigned> get_year_month_day(const time_point& t);


    std::tuple<int, int, int> get_hour_minute_second(const time_point& t);


}
