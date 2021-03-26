#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <voltlbx/smile.h>
#include <voltlbx/time.h>
#include "time.cpp"

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
}


PYBIND11_MODULE(_voltoolbox, m) {

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif

    m.def("check_datetime", &voltlbx_bind::check_datetime, "A datetime test function");
    m.def("check_date", &voltlbx_bind::check_date, "A date test function");


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


    py::class_<SplineSmileCurve, SmileCurve>(m, "SplineSmileCurve")
        .def(py::init<std::vector<double>, std::vector<double>>(),
             py::arg("xs"), py::arg("vols"));

}


