add_library (VolTlbx
             black_scholes.h
             black_scholes.cpp
             cubic_spline.h
             cubic_spline.cpp
             hyp_integral.h
             hyp_integral.cpp
             jet.h             
             linear_filter.h
             linear_filter.cpp
             optim.h
             optim.cpp
             pimpl.h
             roots.h
             roots.cpp
             smile.h
             smile.cpp
             time_utils.h
             time_utils.cpp
             utils.h
             utils.cpp
             hello.h
             hello.cpp)

target_include_directories (VolTlbx
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>) 

target_compile_features(VolTlbx PRIVATE cxx_std_17)

target_link_libraries (VolTlbx
                       PUBLIC Eigen3::Eigen3 date::date
                       PRIVATE cminpack LetsBeRational )
 