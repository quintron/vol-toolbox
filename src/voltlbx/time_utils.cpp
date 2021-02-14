#include "time_utils.h"

#include <date/date.h>

namespace voltlbx
{


    time_point utc_datetime(unsigned year,
                            unsigned month,
                            unsigned day,
                            unsigned hour,
                            unsigned minute,
                            unsigned sec)
    {
        auto ymd = date::month{ month } / day / year;
        auto t_loc = date::local_days{ ymd }
            + std::chrono::hours(hour)
            + std::chrono::minutes(minute)
            + std::chrono::seconds(sec);
        return time_point(t_loc.time_since_epoch());
    }


    std::tuple<int, unsigned, unsigned>
    get_year_month_day(const time_point& t)
    {
        auto t_day = date::floor<date::days>(t);
        const auto ymd = date::year_month_day{ t_day };
        date::year y = ymd.year();
        date::month m = ymd.month();
        date::day d = ymd.day();

        return std::tuple<int, unsigned, unsigned>{ y, m, d };
    }


    std::tuple<int, int, int>
    get_hour_minute_second(const time_point& t)
    {
        auto t_day = date::floor<date::days>(t);
        const auto ymd = date::year_month_day{ t_day };
        const auto time = date::make_time(t - t_day);
        auto h = time.hours().count();
        auto m = time.minutes().count();
        auto s = time.seconds().count();        
        return std::tuple<int, int, int>{ h, m, s };
    }


}
