# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/rrahn/lib_template/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

# This file provides functionality common to the different test modules used by
# lib-template. To build tests, run cmake on one of the sub-folders in this directory
# which contain a CMakeLists.txt.

cmake_minimum_required (VERSION 3.7)

# require lib_template package
find_package (lib_template REQUIRED
              HINTS ${CMAKE_CURRENT_LIST_DIR}/../../build_system)

include (CheckCXXSourceCompiles)
include (FindPackageHandleStandardArgs)
include (FindPackageMessage)

option (LIB_TEMPLATE_TEST_BUILD_OFFLINE "Skip the update step of external projects." OFF)

# Force alignment of benchmarked loops so that numbers are reliable.
# For large loops and erratic seeming bench results the value might
# have to be adapted or the option deactivated.
option (LIB_TEMPLATE_BENCHMARK_ALIGN_LOOPS "Pass -falign-loops=32 to the benchmark builds." ON)

# ----------------------------------------------------------------------------
# Paths to folders.
# ----------------------------------------------------------------------------

find_path (LIB_TEMPLATE_TEST_CMAKE_MODULE_DIR NAMES seqan3_test_component.cmake HINTS "${CMAKE_CURRENT_LIST_DIR}/../../lib/seqan3/test/cmake/")
list(APPEND CMAKE_MODULE_PATH "${LIB_TEMPLATE_TEST_CMAKE_MODULE_DIR}")

set (SEQAN3_BENCHMARK_CLONE_DIR "${PROJECT_BINARY_DIR}/vendor/benchmark")
set (SEQAN3_TEST_CLONE_DIR "${PROJECT_BINARY_DIR}/vendor/googletest")

# needed for add_library (lib_template::test::* INTERFACE IMPORTED)
# see cmake bug https://gitlab.kitware.com/cmake/cmake/issues/15052
file(MAKE_DIRECTORY ${SEQAN3_BENCHMARK_CLONE_DIR}/include/)
file(MAKE_DIRECTORY ${SEQAN3_TEST_CLONE_DIR}/googletest/include/)

# ----------------------------------------------------------------------------
# Interface targets for the different test modules in lib_template.
# ----------------------------------------------------------------------------

# seqan::lib_template::test exposes a base set of required flags, includes, definitions and
# libraries which are in common for **all** lib_template tests
add_library (lib_template_test INTERFACE)
target_compile_options (lib_template_test INTERFACE "-pedantic"  "-Wall" "-Wextra" "-Werror")
target_link_libraries (lib_template_test INTERFACE "seqan::lib_template" "pthread")
target_include_directories (lib_template_test INTERFACE "${LIB_TEMPLATE_TEST_INCLUDE_DIR}")
add_library (seqan::lib_template::test ALIAS lib_template_test)

# lib_template::test::performance specifies required flags, includes and libraries
# needed for performance test cases in lib_template/test/performance
add_library (lib_template_test_performance INTERFACE)
target_link_libraries (lib_template_test_performance INTERFACE "seqan::lib_template::test" "gbenchmark")

if (LIB_TEMPLATE_BENCHMARK_ALIGN_LOOPS)
    target_compile_options (lib_template_test_performance INTERFACE "-falign-loops=32")
endif ()

target_include_directories (lib_template_test_performance INTERFACE "${SEQAN3_BENCHMARK_CLONE_DIR}/include/")
add_library (seqan::lib_template::test::performance ALIAS lib_template_test_performance)

# lib_template::test::unit specifies required flags, includes and libraries
# needed for unit test cases in lib_template/test/unit
add_library (lib_template_test_unit INTERFACE)
target_link_libraries (lib_template_test_unit INTERFACE "seqan::lib_template::test" "gtest_main" "gtest")
target_include_directories (lib_template_test_unit INTERFACE "${SEQAN3_TEST_CLONE_DIR}/googletest/include/")
add_library (seqan::lib_template::test::unit ALIAS lib_template_test_unit)

# lib_template::test::coverage specifies required flags, includes and libraries
# needed for coverage test cases in lib_template/test/coverage
add_library (lib_template_test_coverage INTERFACE)
target_compile_options (lib_template_test_coverage INTERFACE "--coverage" "-fprofile-arcs" "-ftest-coverage")
target_link_libraries (lib_template_test_coverage INTERFACE "seqan::lib_template::test::unit" "gcov")
add_library (seqan::lib_template::test::coverage ALIAS lib_template_test_coverage)

# lib_template::test::header specifies required flags, includes and libraries
# needed for header test cases in lib_template/test/header
add_library (lib_template_test_header INTERFACE)
target_link_libraries (lib_template_test_header INTERFACE "seqan::lib_template::test::unit")
target_link_libraries (lib_template_test_header INTERFACE "seqan::lib_template::test::performance")
target_compile_definitions (lib_template_test_header INTERFACE -DLIB_TEMPLATE_DISABLE_DEPRECATED_WARNINGS)
target_compile_definitions (lib_template_test_header INTERFACE -DLIB_TEMPLATE_HEADER_TEST)
add_library (seqan::lib_template::test::header ALIAS lib_template_test_header)

# ----------------------------------------------------------------------------
# Commonly shared options for external projects.
# ----------------------------------------------------------------------------

set (SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}")

# ----------------------------------------------------------------------------
# Commonly used macros for the different test modules in lib_template.
# ----------------------------------------------------------------------------

include (seqan3_test_component)
include (seqan3_test_files)
include (seqan3_require_ccache)
include (seqan3_require_benchmark)
include (seqan3_require_test)
include (add_subdirectories)
