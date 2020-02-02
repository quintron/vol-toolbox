if(NOT TARGET Eigen3::Eigen3)
    find_path(EIGEN3_INCLUDE_DIRS
              NAMES "Eigen/Dense"
              HINTS "lib/eigen")
    add_library(eigen3 INTERFACE)
    target_include_directories( eigen3 INTERFACE "${EIGEN3_INCLUDE_DIRS}")
    #target_compile_definitions(eigen3 INTERFACE "EIGEN_USE_MKL_ALL")
    target_compile_options(eigen3 INTERFACE "/wd4127" "/wd4146")

    # Starting with msvc 2017 15.8, a byte alignment issue has been resolved.
    # It prevents vectorization of Eigen instructions and can lead to crash.
    if (MSVC AND ${MSVC_VERSION} GREATER_EQUAL 1915)
        target_compile_definitions(eigen3 INTERFACE "_ENABLE_EXTENDED_ALIGNED_STORAGE")
    else()
        target_compile_definitions(eigen3 INTERFACE "_DISABLE_EXTENDED_ALIGNED_STORAGE")
    endif()

    target_sources(eigen3 INTERFACE "${EIGEN3_INCLUDE_DIRS}/debug/msvc/eigen.natvis")
    add_library(Eigen3::Eigen3 ALIAS eigen3)
endif()
