# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.7)
project (pairwise_aligner_test_unit CXX)

include (../cmake/pairwise_aligner-test.cmake)
include (GoogleTest OPTIONAL)

include (diagnostics/list_missing_unit_tests)
include (diagnostics/list_unused_unit_tests)
include (include_dependencies/add_include_dependencies)

option (PAIRWISE_ALIGNER_VERBOSE_TESTS "Run each test case individually" OFF)
option (PAIRWISE_ALIGNER_USE_INCLUDE_DEPENDENCIES "Build tests in an hierarchical order (by an include graph, i.e. tests with less dependencies are build first)" OFF)

if (PAIRWISE_ALIGNER_USE_INCLUDE_DEPENDENCIES)
    cmake_minimum_required (VERSION 3.8)
endif ()

# Adds the test target to the test suite.
# If the variable PAIRWISE_ALIGNER_SDE_TESTS is set, the test is run with the given SDE tests unless the SDE64_BIN is
# not found. In this case no tests will be emitted to the cmake test suite.
# If the variable PAIRWISE_ALIGNER_SDE_TESTS is not set, the test is run natively with out emulation.
macro (add_test_target test_name target)
    find_program (SDE64_BIN sde64)
    if ("${PAIRWISE_ALIGNER_SDE_TESTS}" STREQUAL "")
        message (STATUS "Adding native tests for ${test_name}")
        add_test (NAME "${test_name}" COMMAND ${target})
    else ()
        if ("${SDE64_BIN}" STREQUAL "SDE64_BIN-NOTFOUND")
            message (NOTICE "No sde64 found in PATH. No ctest for ${test_name} generated.")
        else ()
            message (STATUS "Adding emulated tests for ${test_name}")
            set(TARGET_WRAPPER "${SDE64_BIN}" "-${PAIRWISE_ALIGNER_SDE_TESTS}" "--" "./${target}")
            add_test (NAME "${test_name}" COMMAND ${TARGET_WRAPPER})
        endif ()
    endif ()
endmacro ()

macro (pairwise_aligner_test unit_test_cpp)
    cmake_parse_arguments(PAIRWISE_ALIGNER_TEST "" "" "CYCLIC_DEPENDING_INCLUDES" ${ARGN})

    file (RELATIVE_PATH unit_test "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_LIST_DIR}/${unit_test_cpp}")
    seqan3_test_component (target "${unit_test}" TARGET_NAME)
    seqan3_test_component (test_name "${unit_test}" TEST_NAME)

    add_executable (${target} ${unit_test_cpp})
    target_link_libraries (${target} seqan::pairwise_aligner::test::unit)
    add_include_dependencies (${target} "${PAIRWISE_ALIGNER_TEST_CYCLIC_DEPENDING_INCLUDES}")
    collect_used_test (${target})
    if (PAIRWISE_ALIGNER_VERBOSE_TESTS AND NOT CMAKE_VERSION VERSION_LESS 3.10 AND "${PAIRWISE_ALIGNER_SDE_TESTS}" STREQUAL "") # cmake >= 3.10
        gtest_discover_tests(${target} TEST_PREFIX "${test_name}::" PROPERTIES TIMEOUT "30")
    else ()
        add_test_target(${test_name} ${target})
    endif ()

    unset (unit_test)
    unset (target)
    unset (test_name)
endmacro ()

seqan3_require_ccache ()
seqan3_require_test ()

add_subdirectories ()
