# =============================================================================
# algorithms: In-house signal processing + change-point + forecast library
#
# Header-only (we-algorithms is an INTERFACE target — no compiled sources).
# All algorithm headers live in include/algorithms/.
# Tier 51 — no external dependencies.
# =============================================================================

cmake_minimum_required(VERSION 3.24)

project(we-algorithms
    VERSION     1.0.0
    DESCRIPTION "Signal processing, entropy, change-point, and forecast algorithms"
    LANGUAGES   CXX
)

if(NOT WIN32)
    message(FATAL_ERROR "we-algorithms targets Windows x64 only.")
endif()

# Header-only interface library — consumers inherit include path
add_library(we-algorithms INTERFACE)
target_compile_features(we-algorithms INTERFACE cxx_std_23)

target_include_directories(we-algorithms INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# MSVC warning guard — consumers that link we-algorithms get correct /W4 level
if(MSVC)
    target_compile_options(we-algorithms INTERFACE
        /W4 /WX /permissive- /Zc:__cplusplus
    )
endif()
