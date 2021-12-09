#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <voltlbx/utils.h>
#include <voltlbx/time.h>
#include <voltlbx/smile.h>
#include <voltlbx/black_scholes.h>
#include <voltlbx/smile_filter.h>
#include "time.cpp"

#include <vector>

namespace py = pybind11;
using namespace voltlbx;

namespace voltlbx_bind
{
    auto check_datetime(const chrono::DateTime& d)
    {
        return d;
    }
    auto check_date(const chrono::Date& d)
    {
        return d;
    }

    auto longest_increasing_subsequence(const std::vector<double>& xs)
    {
        auto indexes = voltlbx::longest_increasing_subsequence(xs);
        return voltlbx::map(indexes, [](std::size_t i)
        {
            return static_cast<int>(i);
        });
    }
}


PYBIND11_MODULE(_voltoolbox, m) {

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif

    m.def("check_datetime", &voltlbx_bind::check_datetime, "A datetime test function");
    m.def("check_date", &voltlbx_bind::check_date, "A date test function");

    m.def("longest_increasing_subsequence",
          &voltlbx_bind::longest_increasing_subsequence,
          "Return indexes of longest increasing subsequence");

    py::class_<Calendar, std::shared_ptr<Calendar>>(m, "Calendar")
        .def(py::init<std::vector<chrono::Date>>(), py::arg("holidays"))
        .def("is_closed", &Calendar::is_closed, py::arg("date"))
        .def("count_open_days",
            [](const Calendar& cal, const chrono::Date& start, const chrono::Date& end)
            {
                return cal.count_open_days(start, end);
            });


    py::class_<BusinessTimeMeasure>(m, "BusinessTimeMeasure")
        .def(py::init<std::shared_ptr<Calendar>, double, double>(),
             py::arg("calendar"), py::arg("close_weight"), py::arg("yearly_nb_open"))
        .def("distance", &BusinessTimeMeasure::distance, py::arg("t0"), py::arg("t1"));


    py::class_<SmileCurve>(m, "Smile")
        .def("vol", &SmileCurve::vol, py::arg("x"))
        .def("vol_jet", 
            [](const SmileCurve& s, double x)
            {
                auto j = s.vol_jet(x);
                return std::make_tuple(j.y, j.dy_dx, j.d2y_d2x);
            },
            py::arg("x"));


    py::class_<SmileSplineCurve, SmileCurve>(m, "SmileSplineCurve")
        .def(py::init<std::vector<double>, std::vector<double>>(),
             py::arg("xs"), py::arg("vols"));


    py::class_<SplineCurve, SmileCurve>(m, "SplineCurve")
        .def(py::init<std::vector<double>, std::vector<double>>(),
             py::arg("xs"), py::arg("vols"));


    py::class_<AverageSplineCurve, SmileCurve>(m, "AverageSplineCurve")
        .def(py::init<std::vector<double>, std::vector<double>>(),
             py::arg("xs"), py::arg("vs"));


    m.def("bs_implied_volatility", &voltlbx::bs_implied_volatility,  
          py::arg("forward"), py::arg("strike"), py::arg("price"), py::arg("time"), py::arg("option_type"),
          "Compute Black-Scholes implied volatility");

    py::class_<SmileVariationFilter>(m, "SmileVariationFilter")
        .def(py::init(&SmileVariationFilter::create),
             py::arg("zs"), py::arg("dvols"), py::arg("error_devs"), 
             py::arg("atm_dev"), py::arg("atm_skew_dev"), py::arg("z_ref"))
        .def("dvol_with_error", &SmileVariationFilter::dvol_with_error, py::arg("z"));

}
