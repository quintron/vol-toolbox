add_library (VolTlbx
             cubic_spline.h
             cubic_spline.cpp
             hyp_integrale.h
             hyp_integrale.cpp
             roots.h
             roots.cpp
             utils.h
             utils.cpp
             hello.h
             hello.cpp)

target_include_directories (VolTlbx
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>) 

target_compile_features(VolTlbx PRIVATE cxx_std_17)
                       
target_link_libraries (VolTlbx
                       PRIVATE Eigen3::Eigen3)
