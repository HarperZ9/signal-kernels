# =============================================================================
# signal-kernels: header-only signal-processing, entropy, change-point,
# information-theory, curvature, and forecast algorithms.
#
# Header-only (signal-kernels is an INTERFACE target — no compiled sources).
# All algorithm headers live in include/algorithms/.
# No external runtime dependencies (C++ standard library only).
# =============================================================================

cmake_minimum_required(VERSION 3.24)

project(signal-kernels
    VERSION     1.0.0
    DESCRIPTION "Signal processing, entropy, change-point, and forecast algorithms"
    LANGUAGES   CXX
)

if(NOT WIN32)
    message(FATAL_ERROR "signal-kernels targets Windows x64 only.")
endif()

# Header-only interface library — consumers inherit include path
add_library(signal-kernels INTERFACE)
add_library(signal-kernels::signal-kernels ALIAS signal-kernels)
target_compile_features(signal-kernels INTERFACE cxx_std_23)

target_include_directories(signal-kernels INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# MSVC warning guard — consumers that link signal-kernels get correct /W4 level
if(MSVC)
    target_compile_options(signal-kernels INTERFACE
        /W4 /WX /permissive- /Zc:__cplusplus
    )
endif()

# -----------------------------------------------------------------------------
# Tests — built and run via CTest when this is the top-level project.
# Toggle with -DSIGNAL_KERNELS_BUILD_TESTS=ON/OFF.
# -----------------------------------------------------------------------------
option(SIGNAL_KERNELS_BUILD_TESTS
    "Build the signal-kernels unit tests"
    ${PROJECT_IS_TOP_LEVEL}
)

if(SIGNAL_KERNELS_BUILD_TESTS)
    enable_testing()

    # Fetch doctest (single-header test framework). If a vendored copy exists at
    # tests/third_party/doctest/doctest.h it is used instead of downloading.
    set(_vendored_doctest
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/third_party/doctest/doctest.h")

    if(EXISTS "${_vendored_doctest}")
        add_library(doctest INTERFACE)
        target_include_directories(doctest INTERFACE
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/third_party")
        add_library(doctest::doctest ALIAS doctest)
        message(STATUS "signal-kernels: using vendored doctest at ${_vendored_doctest}")
    else()
        include(FetchContent)
        FetchContent_Declare(doctest
            GIT_REPOSITORY https://github.com/doctest/doctest.git
            GIT_TAG        v2.4.11
        )
        FetchContent_MakeAvailable(doctest)
    endif()

    set(_sk_test_sources
        tests/test_algorithms_numeric.cpp
        tests/test_algorithms_fft.cpp
        tests/test_algorithms_entropy.cpp
        tests/test_algorithms_information.cpp
        tests/test_algorithms_causal.cpp
        tests/test_algorithms_changepoint.cpp
        tests/test_algorithms_curvature.cpp
        tests/test_algorithms_forecast.cpp
    )

    add_executable(signal-kernels-tests
        tests/doctest_main.cpp
        ${_sk_test_sources}
    )
    target_link_libraries(signal-kernels-tests
        PRIVATE signal-kernels doctest::doctest
    )

    add_test(NAME signal-kernels-tests COMMAND signal-kernels-tests)
endif()

