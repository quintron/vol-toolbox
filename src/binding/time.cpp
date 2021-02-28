#include <pybind11/pybind11.h>
#include <pybind11/chrono.h>

#include <voltlbx/time_utils.h>


namespace pybind11 {
    namespace detail {

        template <> class type_caster<chrono::DateTime> {
        public:

            bool load(handle src, bool) {
                using namespace std::chrono;

                // Lazy initialise the PyDateTime import
                if (!PyDateTimeAPI) { PyDateTime_IMPORT; }

                if (!src) return false;

                unsigned year, month, day, hour, minute, sec;

                if (PyDateTime_Check(src.ptr())) {
                    sec = PyDateTime_DATE_GET_SECOND(src.ptr());
                    minute = PyDateTime_DATE_GET_MINUTE(src.ptr());
                    hour = PyDateTime_DATE_GET_HOUR(src.ptr());
                    day = PyDateTime_GET_DAY(src.ptr());
                    month = PyDateTime_GET_MONTH(src.ptr());
                    year = PyDateTime_GET_YEAR(src.ptr());
                }
                else if (PyDate_Check(src.ptr())) {
                    sec = 0;
                    minute = 0;
                    hour = 0;
                    day = PyDateTime_GET_DAY(src.ptr());
                    month = PyDateTime_GET_MONTH(src.ptr());
                    year = PyDateTime_GET_YEAR(src.ptr());
                }
                else return false;

                value = voltlbx::utc_datetime(year, month, day, hour, minute, sec);
                return true;
            }

            static handle cast(const chrono::DateTime& src, return_value_policy /* policy */, handle /* parent */) {
                using namespace std::chrono;

                // Lazy initialise the PyDateTime import
                if (!PyDateTimeAPI) { PyDateTime_IMPORT; }

                auto [year, month, day] = voltlbx::get_year_month_day(src);
                auto[hour, minute, sec] = voltlbx::get_hour_minute_second(src);

                return PyDateTime_FromDateAndTime(year, month, day, hour, minute, sec, 0);
            }
            PYBIND11_TYPE_CASTER(chrono::DateTime, _("datetime.datetime"));
        };

        template <> class type_caster<chrono::Date> {
        public:

            bool load(handle src, bool) {
                using namespace std::chrono;

                // Lazy initialise the PyDateTime import
                if (!PyDateTimeAPI) { PyDateTime_IMPORT; }

                if (!src) return false;

                unsigned year, month, day;

                if (PyDate_Check(src.ptr())) {
                    day = PyDateTime_GET_DAY(src.ptr());
                    month = PyDateTime_GET_MONTH(src.ptr());
                    year = PyDateTime_GET_YEAR(src.ptr());
                }
                else return false;

                value = voltlbx::create_date(year, month, day);
                return true;
            }

            static handle cast(const chrono::Date& src, return_value_policy /* policy */, handle /* parent */) {
                using namespace std::chrono;

                // Lazy initialise the PyDateTime import
                if (!PyDateTimeAPI) { PyDateTime_IMPORT; }

                auto [year, month, day] = voltlbx::get_year_month_day(src);
                return PyDate_FromDate(year, month, day);
            }
            PYBIND11_TYPE_CASTER(chrono::Date, _("datetime.date"));
        };

    }
}
