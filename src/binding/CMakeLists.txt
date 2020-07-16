pybind11_add_module(voltoolbox 
                    voltoolbox.cpp)
target_link_libraries (voltoolbox PUBLIC VolTlbx Eigen3::Eigen3)
set_target_properties (voltoolbox PROPERTIES FOLDER binding)
