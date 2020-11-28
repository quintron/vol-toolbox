add_library (LetsBeRational
             erf_cody.cpp
             importexport.h
             lets_be_rational.cpp
             lets_be_rational.h             
             normaldistribution.cpp
             normaldistribution.h
             rationalcubic.cpp
             rationalcubic.h
             version.h
             )

target_include_directories (LetsBeRational
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>) 