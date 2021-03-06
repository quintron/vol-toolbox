#pragma once

#include <chrono>
#include <tuple>
#include <vector>
#include <date/date.h>

#include "pimpl.h"

namespace chrono
{
    using DateTime = date::local_seconds;
    using Date = date::local_days;
}


namespace voltlbx
{
    chrono::Date to_date(const chrono::DateTime& t);

    chrono::Date create_date(unsigned year,
                             unsigned month,
                             unsigned day);

    chrono::DateTime utc_datetime(unsigned year,
                                  unsigned month,
                                  unsigned day,
                                  unsigned hour=0,
                                  unsigned minute=0,
                                  unsigned sec=0);

    std::tuple<int, unsigned, unsigned> get_year_month_day(const chrono::Date& d);

    inline std::tuple<int, unsigned, unsigned> get_year_month_day(const chrono::DateTime& t)
    {
        return get_year_month_day(to_date(t));
    }

    std::tuple<int, int, int> get_hour_minute_second(const chrono::DateTime& t);


    class Calendar: public Pimpl<Calendar>
    {
    public:
        Calendar(const std::vector<date::weekday>& close_week_days,
                 const std::vector<chrono::Date>& holidays);

        Calendar(const std::vector<chrono::Date>& holidays);

        int count_open_days(const chrono::Date& start, const chrono::Date& end) const;

        bool is_closed(const chrono::Date d) const;
    };

}
