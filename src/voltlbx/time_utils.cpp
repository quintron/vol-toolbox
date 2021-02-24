#include "time_utils.h"

#include <array>
#include <set>

namespace voltlbx
{
    using namespace chrono;

    Date to_date(const DateTime& t) 
    {
        return date::floor<date::days>(t);
    }

    Date create_date(unsigned year,
                     unsigned month,
                     unsigned day)
    {
        auto ymd = date::month{ month } / day / year;
        return Date{ ymd };
    }


    DateTime utc_datetime(unsigned year,
                            unsigned month,
                            unsigned day,
                            unsigned hour,
                            unsigned minute,
                            unsigned sec)
    {        
        auto t_loc = create_date(year, month, day)
                    + std::chrono::hours(hour)
                    + std::chrono::minutes(minute)
                    + std::chrono::seconds(sec);
        return DateTime(t_loc.time_since_epoch());
    }


    std::tuple<int, unsigned, unsigned>
    get_year_month_day(const Date& t)
    {
        const auto ymd = date::year_month_day{ t };
        date::year y = ymd.year();
        date::month m = ymd.month();
        date::day d = ymd.day();
        return std::tuple<int, unsigned, unsigned>{ y, m, d };
    }


    std::tuple<int, int, int>
    get_hour_minute_second(const DateTime& t)
    {
        const auto time = date::make_time(t - to_date(t));
        auto h = time.hours().count();
        auto m = time.minutes().count();
        auto s = time.seconds().count();        
        return std::tuple<int, int, int>{ h, m, s };
    }


    template<>
    struct Pimpl<Calendar>::Implementation
    {
        explicit Implementation(const std::vector<date::weekday>& close_week_days, 
                                const std::vector<Date>& holidays)
        {
            for (std::size_t i = 0; i < 7; ++i)
            {
                closed_weekdays[i] = false;
            }
            for (auto w : close_week_days)
            {
                assert(w.ok());
                closed_weekdays[w.c_encoding()] = true;
            }
            
            nb_open_weekdays = 0;
            for (auto is_closed : closed_weekdays)
            {
                nb_open_weekdays += !is_closed;
            }

            for (int i = 0; i < 7; i++)
            {
                int k = 0;
                for (int j = i; j < i + 7; j++)
                {
                    weekday_adjustments[i][j % 7] = k;

                    if (!closed_weekdays[j % 7])
                        k++;
                }
            }

            for (const auto& h : holidays)
            {
                // filter open week days
                auto wd = date::year_month_weekday{ h }.weekday();
                if (!closed_weekdays[wd.c_encoding()])
                {
                    week_holidays.insert(h);
                }
            }

        }

        unsigned count_open_days_no_holidays(const Date& start, 
                                             const Date& end) const
        {
            if (nb_open_weekdays == 0)
                return 0;

            auto ymd1 = date::year_month_weekday{ start };
            auto ymd2 = date::year_month_weekday{ end };

            auto d1 = ymd1.weekday().c_encoding();
            auto d2 = ymd2.weekday().c_encoding();

            int weeks = (end - start).count() / 7;
            int adjustment = weekday_adjustments[d1][d2];

            return nb_open_weekdays * weeks + adjustment;
        }

        int count_holidays(const Date& start,
                                const Date& end) const
        {
            auto i0 = week_holidays.lower_bound(start);
            auto i1 = week_holidays.lower_bound(end);
            return std::distance(i0, i1);
        }

        int open_days_count(const Date& start, const Date& end) const
        {
            if(start > end)
                throw std::out_of_range("start should be previous to end");

            if (start == end)
                return 0;

            int open_days = count_open_days_no_holidays(start, end);
            open_days -= count_holidays(start, end);

            return open_days;
        }

        std::array<bool, 7> closed_weekdays;
        int nb_open_weekdays;
        std::array<std::array<int, 7>, 7> weekday_adjustments;
        std::set<Date> week_holidays;
    };


    Calendar::Calendar(const std::vector<date::weekday>& close_week_days, 
                       const std::vector<Date>& holidays)
        : Pimpl<Calendar>(close_week_days, holidays)
    {
    }

    int Calendar::open_days_count(const Date& start, const Date& end) const
    {
        return impl->open_days_count(start, end);
    }

}
