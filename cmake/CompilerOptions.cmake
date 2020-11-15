set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (MSVC)
    # All warnings
    add_compile_options ("/W4")
    
    # Parallel build
    add_compile_options ("/MP")
else()
    # Need for pybind11  
    add_compile_options (-D_hypot=hypot -fPIC)
    
    # Need to compile with mingw in debug mode 
    if (WIN32)
         add_compile_options (-Wa,-mbig-obj)
    endif()

endif ()

#Useful defines
if (WIN32)
    add_definitions (-DSTRICT)
    add_definitions (-DNOMINMAX)
endif ()

if (MSVC)
    add_definitions (-D_SCL_SECURE_NO_WARNINGS)
    add_definitions (-D_CRT_SECURE_NO_WARNINGS)
endif ()

# Use folders in Visual Studio solutions
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
