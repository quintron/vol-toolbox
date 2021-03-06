#include "time.h"
using namespace std::chrono;


namespace voltlbx
{
    using namespace chrono;

    template<>
    struct Pimpl<BusinessTimeMeasure>::Implementation
    {
        Implementation(std::shared_ptr<Calendar> calendar,
            double close_weight,
            double yearly_nb_open)
            : calendar(std::move(calendar)),
            close_weight(close_weight),
            yearly_nb_open(yearly_nb_open)
        {}

        double open_to_open(const Date& d0, const Date& d1) const
        {
            const int open_days = calendar->count_open_days(d0, d1);
            const int total_days = std::chrono::duration_cast<date::days>(d1 - d0).count();
            const int close_days = total_days - open_days;
            return open_days + close_days * close_weight;
        }

        double distance_to_open(const DateTime& t) const
        {
            using namespace std::chrono;

            const auto d = to_date(t);
            const double weight = calendar->is_closed(d) ? close_weight : 1.0;
            return weight * duration_cast<seconds>(t - d).count() / (24 * 3600.0);
        }

        double norm() const
        {
            return yearly_nb_open + close_weight * (365.0 - yearly_nb_open);
        }

        double distance(const DateTime& t0, const DateTime& t1) const
        {
            const auto d0 = to_date(t0);
            const auto d1 = to_date(t1);
            double dist = open_to_open(d0, d1)
                            - distance_to_open(t0)
                            + distance_to_open(t1);
            return dist / norm();
        }

        const std::shared_ptr<Calendar> calendar;
        const double close_weight;
        const double yearly_nb_open;
    };

    BusinessTimeMeasure::BusinessTimeMeasure
        (std::shared_ptr<Calendar> calendar,
         double close_weight,
         double yearly_nb_open)
        :Pimpl<BusinessTimeMeasure>(calendar, close_weight, yearly_nb_open)
    {
    }

    double BusinessTimeMeasure::distance(const DateTime& t0, const DateTime& t1) const
    {
        return impl->distance(t0, t1);
    }

}