package_add_test (voltlbx_test
                  test_asinh_spline.cpp
                  test_cubic_spline.cpp 
                  test_roots.cpp
                  test_utils.cpp				  
                  test_add.cpp)
target_compile_features(voltlbx_test PRIVATE cxx_std_17)
target_link_libraries (voltlbx_test VolTlbx)
