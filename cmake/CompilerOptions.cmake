set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (MSVC)
    # All warnings
    add_compile_options ("/W4")
    
    # Parallel build
    add_compile_options ("/MP")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_hypot=hypot ")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj ")    
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
