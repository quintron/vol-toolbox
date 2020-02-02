cmake_minimum_required(VERSION 3.8)
project (vol-toolbox VERSION 1.0 LANGUAGES CXX)

#Update repo submodules
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

#VS project options
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

#Global compile options
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#Include
include(cmake/AddTest.cmake)

#Lib
add_subdirectory (src/voltlbx)

#Tests
option(BUILD_TESTING "Test building option" ON)
if(BUILD_TESTING)
    #Google test
    set (gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    add_subdirectory (lib/googletest)
    set_target_properties (gtest PROPERTIES FOLDER external)
    set_target_properties (gtest_main PROPERTIES FOLDER external)
    set_target_properties (gmock PROPERTIES FOLDER external)
    set_target_properties (gmock_main PROPERTIES FOLDER external)
    
    enable_testing ()
    add_subdirectory (tests)
endif()

