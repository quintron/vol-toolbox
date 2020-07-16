cmake_minimum_required(VERSION 3.8)
project (vol-toolbox VERSION 1.0 LANGUAGES CXX)

#CMake tools
list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
include (UpdateGitSubmodule)
include (CompilerOptions)
include (AddTest)

#Lib
find_package (Eigen3 REQUIRED)
find_package (CMinPack REQUIRED)
add_subdirectory (lib/pybind11)
add_subdirectory (src/voltlbx)
add_subdirectory (src/binding)

#Tests
option (BUILD_TESTING "Test building option" ON)
if (BUILD_TESTING)
    #Google test
    set (gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    add_subdirectory (lib/googletest)
    set_target_properties (gtest PROPERTIES FOLDER external)
    set_target_properties (gtest_main PROPERTIES FOLDER external)
    set_target_properties (gmock PROPERTIES FOLDER external)
    set_target_properties (gmock_main PROPERTIES FOLDER external)    
    
    #Project tests
    enable_testing ()
    add_subdirectory (tests)
endif()

#Install
set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/install")
install(TARGETS voltoolbox DESTINATION "${CMAKE_CURRENT_LIST_DIR}/dist")
