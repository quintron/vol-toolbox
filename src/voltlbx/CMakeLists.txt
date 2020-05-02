add_library (voltlbx
             cubic_spline.h
             cubic_spline.cpp
             roots.h
             roots.cpp
             utils.h
             utils.cpp
             hello.h
             hello.cpp)
set_target_properties (voltlbx PROPERTIES 
                       LINKER_LANGUAGE CXX
                       FOLDER src)
target_include_directories (voltlbx PUBLIC
                            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>)                        
target_link_libraries (voltlbx
                       PRIVATE eigen3)
add_dependencies (voltlbx eigen3)	