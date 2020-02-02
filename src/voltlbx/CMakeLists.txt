add_library (voltlbx
             hello.h
             hello.cpp
             )
set_target_properties (voltlbx PROPERTIES 
                       LINKER_LANGUAGE CXX
                       FOLDER src)
target_include_directories (voltlbx PUBLIC
                            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
                           )
                          
