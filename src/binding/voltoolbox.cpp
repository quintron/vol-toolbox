#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <voltlbx/smile.h>

namespace py = pybind11;
using namespace voltlbx;

PYBIND11_MODULE(voltoolbox, m) {

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif

    py::class_<Smile>(m, "Smile")
        .def("vol", &Smile::vol);

    py::class_<CubicSplineSmile, Smile>(m, "CubicSplineSmile")
        .def(py::init<double, double, std::vector<double>, std::vector<double>>());

}