#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/chrono.h>
#include <voltlbx/smile.h>

namespace py = pybind11;
using namespace voltlbx;

PYBIND11_MODULE(_voltoolbox, m) {

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif

    py::class_<Smile>(m, "Smile")
        .def("vol", &Smile::vol, py::arg("x"))
        .def("density_ratio", &Smile::density_ratio, py::arg("x"))
        .def("pseudo_local_vol", &Smile::pseudo_local_vol, py::arg("x"));

    py::class_<CubicSplineSmile, Smile>(m, "CubicSplineSmile")
        .def(py::init<double, double, std::vector<double>, std::vector<double>>(),
             py::arg("time_to_maturity"), py::arg("atf_vol"), py::arg("zs"), py::arg("vol_ratios"));

}