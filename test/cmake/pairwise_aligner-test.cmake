# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

# This file provides functionality common to the different test modules used by
# lib-template. To build tests, run cmake on one of the sub-folders in this directory
# which contain a CMakeLists.txt.

cmake_minimum_required (VERSION 3.7)

# require pairwise_aligner package
find_package (pairwise_aligner REQUIRED
              HINTS ${CMAKE_CURRENT_LIST_DIR}/../../build_system)

include (CheckCXXSourceCompiles)
include (FindPackageHandleStandardArgs)
include (FindPackageMessage)

option (PAIRWISE_ALIGNER_TEST_BUILD_OFFLINE "Skip the update step of external projects." OFF)

# Force alignment of benchmarked loops so that numbers are reliable.
# For large loops and erratic seeming bench results the value might
# have to be adapted or the option deactivated.
option (PAIRWISE_ALIGNER_BENCHMARK_ALIGN_LOOPS "Pass -falign-loops=32 to the benchmark builds." ON)

# ----------------------------------------------------------------------------
# Paths to folders.
# ----------------------------------------------------------------------------

find_path (SEQAN3_TEST_INCLUDE_DIR NAMES seqan3/test/tmp_filename.hpp HINTS "${CMAKE_CURRENT_LIST_DIR}/../../lib/seqan3/test/include/")
find_path (PAIRWISE_ALIGNER_TEST_CMAKE_MODULE_DIR NAMES seqan3_test_component.cmake HINTS "${CMAKE_CURRENT_LIST_DIR}/../../lib/seqan3/test/cmake/")
list(APPEND CMAKE_MODULE_PATH "${PAIRWISE_ALIGNER_TEST_CMAKE_MODULE_DIR}" "${CMAKE_CURRENT_LIST_DIR}")

message (STATUS "[DEBUG] CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}")

# set (SEQAN3_BENCHMARK_CLONE_DIR "${PROJECT_BINARY_DIR}/vendor/benchmark")
# set (SEQAN3_TEST_CLONE_DIR "${PROJECT_BINARY_DIR}/vendor/googletest")
# set (LIBUNIFEX_CLONE_DIR "${PROJECT_BINARY_DIR}/vendor/libunifex")

# needed for add_library (pairwise_aligner::test::* INTERFACE IMPORTED)
# see cmake bug https://gitlab.kitware.com/cmake/cmake/issues/15052
# file(MAKE_DIRECTORY ${SEQAN3_BENCHMARK_CLONE_DIR}/include/)
# file(MAKE_DIRECTORY ${SEQAN3_TEST_CLONE_DIR}/googletest/include/)
# file(MAKE_DIRECTORY ${LIBUNIFEX_CLONE_DIR}/include/)

# ----------------------------------------------------------------------------
# Interface targets for the different test modules in pairwise_aligner.
# ----------------------------------------------------------------------------

# seqan::pairwise_aligner::test exposes a base set of required flags, includes, definitions and
# libraries which are in common for **all** pairwise_aligner tests
if (NOT TARGET seqan::pairwise_aligner::test)
    add_library (pairwise_aligner_test INTERFACE)
    target_compile_options (pairwise_aligner_test INTERFACE "-pedantic"  "-Wall" "-Wextra" "-Werror")
    target_link_libraries (pairwise_aligner_test INTERFACE "seqan::pairwise_aligner" "pthread" "unifex")
    target_include_directories (pairwise_aligner_test INTERFACE "${SEQAN3_TEST_INCLUDE_DIR}")
    add_library (seqan::pairwise_aligner::test ALIAS pairwise_aligner_test)
endif ()

# seqan::pairwise_aligner::test::performance specifies required flags, includes and libraries
# needed for performance test cases in pairwise_aligner/test/performance
if (NOT TARGET seqan::pairwise_aligner::test::performance)
    add_library (pairwise_aligner_test_performance INTERFACE)
    target_link_libraries (pairwise_aligner_test_performance INTERFACE "seqan::pairwise_aligner::test" "gbenchmark")

    if (PAIRWISE_ALIGNER_BENCHMARK_ALIGN_LOOPS)
        target_compile_options (pairwise_aligner_test_performance INTERFACE "-falign-loops=32")
    endif ()

    add_library (seqan::pairwise_aligner::test::performance ALIAS pairwise_aligner_test_performance)
endif ()

# seqan::pairwise_aligner::test::unit specifies required flags, includes and libraries
# needed for unit test cases in pairwise_aligner/test/unit
if (NOT TARGET seqan::pairwise_aligner::test::unit)
    add_library (pairwise_aligner_test_unit INTERFACE)
    target_link_libraries (pairwise_aligner_test_unit INTERFACE "seqan::pairwise_aligner::test" "gtest_main" "gtest")
    add_library (seqan::pairwise_aligner::test::unit ALIAS pairwise_aligner_test_unit)
endif ()

# seqan::pairwise_aligner::test::coverage specifies required flags, includes and libraries
# needed for coverage test cases in pairwise_aligner/test/coverage
if (NOT TARGET seqan::pairwise_aligner::test::coverage)
    add_library (pairwise_aligner_test_coverage INTERFACE)
    target_compile_options (pairwise_aligner_test_coverage INTERFACE "--coverage" "-fprofile-arcs" "-ftest-coverage")
    # -fprofile-abs-path requires at least gcc8, it forces gcov to report absolute instead of relative paths.
    # gcovr has trouble detecting the headers otherwise.
    # ccache is not aware of this option, so it needs to be skipped with `--ccache-skip`.
    find_program (CCACHE_PROGRAM ccache)
    if (CCACHE_PROGRAM)
        target_compile_options (pairwise_aligner_test_coverage INTERFACE "--ccache-skip" "-fprofile-abs-path")
    else ()
        target_compile_options (pairwise_aligner_test_coverage INTERFACE "-fprofile-abs-path")
    endif ()
    target_link_libraries (pairwise_aligner_test_coverage INTERFACE "seqan::pairwise_aligner::test::unit" "gcov")
    add_library (seqan::pairwise_aligner::test::coverage ALIAS pairwise_aligner_test_coverage)
endif ()

# seqan::pairwise_aligner::test::header specifies required flags, includes and libraries
# needed for header test cases in pairwise_aligner/test/header
if (NOT TARGET seqan::pairwise_aligner::test::header)
    add_library (pairwise_aligner_test_header INTERFACE)
    target_link_libraries (pairwise_aligner_test_header INTERFACE "seqan::pairwise_aligner::test::unit")
    target_link_libraries (pairwise_aligner_test_header INTERFACE "seqan::pairwise_aligner::test::performance")
    target_compile_definitions (pairwise_aligner_test_header INTERFACE -DPAIRWISE_ALIGNER_DISABLE_DEPRECATED_WARNINGS)
    target_compile_definitions (pairwise_aligner_test_header INTERFACE -DPAIRWISE_ALIGNER_HEADER_TEST)
    add_library (seqan::pairwise_aligner::test::header ALIAS pairwise_aligner_test_header)
endif ()

# ----------------------------------------------------------------------------
# Commonly shared options for external projects.
# ----------------------------------------------------------------------------

set (SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}")

# ----------------------------------------------------------------------------
# Commonly used macros for the different test modules in pairwise_aligner.
# ----------------------------------------------------------------------------

# libunifex is not called
include (configure_libunifex)
include (seqan3_test_component)
include (seqan3_test_files)
include (seqan3_require_ccache)
include (seqan3_require_benchmark)
include (seqan3_require_test)
include (add_subdirectories)
