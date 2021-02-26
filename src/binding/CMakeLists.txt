pybind11_add_module(_voltoolbox
                    time.cpp
                    voltoolbox.cpp)
target_link_libraries (_voltoolbox PUBLIC VolTlbx Eigen3::Eigen3 date::date)
set_target_properties (_voltoolbox PROPERTIES FOLDER binding)
