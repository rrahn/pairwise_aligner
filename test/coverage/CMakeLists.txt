# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.10)
project (pairwise_aligner_test_coverage CXX)

include (../cmake/pairwise_aligner-test.cmake)

if (CMAKE_BUILD_TYPE AND NOT ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug"))
    message(WARNING "Coverage test must be build in debug mode [build type = ${CMAKE_BUILD_TYPE}]")
endif ()

find_program (GCOVR_COMMAND NAMES gcovr)

if (NOT GCOVR_COMMAND)
    message(FATAL_ERROR "gcovr not found! Aborting...")
endif ()

# Holds all target's defined by pairwise_aligner_test
set_property (GLOBAL PROPERTY GLOBAL_TEST_COVERAGE_ALL_TESTS "")

set (PAIRWISE_ALIGNER_COVERAGE_PARALLEL_LEVEL "1" CACHE STRING "Number of threads to use for coverage report generation.")

# Arguments for gcovr.
# Each line is an argument, ';' will end up being a whitespace due to expansion by CMake.
# Special characters like '"', '\', ')' need to be escaped by a number of '\', depending how often they are evaluated.
# gcovr uses python regex expressions.
# See https://gcovr.com/en/5.0/guide.html#the-gcovr-command for an overview of parameters.
# See https://gcovr.com/en/5.0/faq.html#why-does-c-code-have-so-many-uncovered-branches for an explanation on branches.
set (PAIRWISE_ALIGNER_GCOVR_ARGUMENTS "")
# The directory of the CMakeLists.txt used for invoking cmake.
list (APPEND PAIRWISE_ALIGNER_GCOVR_ARGUMENTS "--root;${CMAKE_CURRENT_LIST_DIR}")
# The build directory.
list (APPEND PAIRWISE_ALIGNER_GCOVR_ARGUMENTS "${PROJECT_BINARY_DIR}")
# Include all files whose path match '${PAIRWISE_ALIGNER_CLONE_DIR}/include/pairwise_aligner/.*'.
list (APPEND PAIRWISE_ALIGNER_GCOVR_ARGUMENTS "--filter;\'${PAIRWISE_ALIGNER_CLONE_DIR}/include/pairwise_aligner/.*\'")
# Include all files whose path match '${PAIRWISE_ALIGNER_CLONE_DIR}/test/include/pairwise_aligner/test/.*'.
list (APPEND PAIRWISE_ALIGNER_GCOVR_ARGUMENTS "--filter;\'${PAIRWISE_ALIGNER_CLONE_DIR}/test/include/pairwise_aligner/test/.*\'")
# Include all files whose path match '${PAIRWISE_ALIGNER_CLONE_DIR}/lib/seqan3/include/seqan3/.*'.
list (APPEND SEQAN3_GCOVR_ARGUMENTS "--filter;\'${PAIRWISE_ALIGNER_CLONE_DIR}/lib/seqan3/include/seqan3/.*\'")
# Include all files whose path match '${PAIRWISE_ALIGNER_CLONE_DIR}/lib/seqan3/test/include/seqan3/test/.*'.
list (APPEND SEQAN3_GCOVR_ARGUMENTS "--filter;\'${PAIRWISE_ALIGNER_CLONE_DIR}/lib/seqan3/test/include/seqan3/test/.*\'")
# Remove all files whose path match '${PAIRWISE_ALIGNER_CLONE_DIR}/lib/seqan3/include/seqan3/contrib/.*'.
list (APPEND SEQAN3_GCOVR_ARGUMENTS "--exclude;\'${PAIRWISE_ALIGNER_CLONE_DIR}/lib/seqan3/include/seqan3/contrib/.*\'")
# Remove all files whose path match '${PAIRWISE_ALIGNER_CLONE_DIR}/lib/seqan3/include/seqan3/std/.*'.
list (APPEND SEQAN3_GCOVR_ARGUMENTS "--exclude;\'${PAIRWISE_ALIGNER_CLONE_DIR}/lib/seqan3/include/seqan3/std/.*\'")
# Remove line coverage for all lines that match '^\s*$', i.e. empty lines.
list (APPEND PAIRWISE_ALIGNER_GCOVR_ARGUMENTS "--exclude-lines-by-pattern;\\\'^\\\\s*$$\\\'")
# Remove line coverage for all lines that match "^\s*[\{\}]\s*\;*\s*(//.*)?$", i.e. lines with a single '{' or '}' which
# might be followed by a semicolon or a comment '//'.
list (APPEND PAIRWISE_ALIGNER_GCOVR_ARGUMENTS "--exclude-lines-by-pattern;\\\'^\\\\s*[\\{\\}]\\\\s*\\\\\\;*\\\\s*\\(//.*\\)?$$\\\'")
# Will exclude branches that are unreachable.
list (APPEND PAIRWISE_ALIGNER_GCOVR_ARGUMENTS "--exclude-unreachable-branches")
# Will exclude branches that are only generated for exception handling.
list (APPEND PAIRWISE_ALIGNER_GCOVR_ARGUMENTS "--exclude-throw-branches")
# Run up to this many gcov instances in parallel.
list (APPEND PAIRWISE_ALIGNER_GCOVR_ARGUMENTS "-j;${PAIRWISE_ALIGNER_COVERAGE_PARALLEL_LEVEL}")

add_custom_command (
    OUTPUT ${PROJECT_BINARY_DIR}/pairwise_aligner_coverage.xml

    # Run tests.
    COMMAND ${CMAKE_CTEST_COMMAND} -j ${PAIRWISE_ALIGNER_COVERAGE_PARALLEL_LEVEL} --output-on-failure

    # Run gcovr and create XML report.
    COMMAND ${GCOVR_COMMAND} ${PAIRWISE_ALIGNER_GCOVR_ARGUMENTS} --xml --output ${PROJECT_BINARY_DIR}/pairwise_aligner_coverage.xml

    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    COMMENT "Processing code coverage counters and generating XML report."
    COMMAND_EXPAND_LISTS
)

add_custom_target (
    coverage ALL
    DEPENDS
    ${PROJECT_BINARY_DIR}/pairwise_aligner_coverage.xml
)

add_custom_command (
    OUTPUT ${PROJECT_BINARY_DIR}/html/index.html

    # Run tests.
    COMMAND ${CMAKE_CTEST_COMMAND} -j ${PAIRWISE_ALIGNER_COVERAGE_PARALLEL_LEVEL} --output-on-failure

    # Create output directory.
    COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/html/

    # Run gcovr and create HTML report.
    COMMAND ${GCOVR_COMMAND} ${PAIRWISE_ALIGNER_GCOVR_ARGUMENTS} --html-details --output ${PROJECT_BINARY_DIR}/html/index.html

    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    COMMENT "Processing code coverage counters and generating HTML report."
    COMMAND_EXPAND_LISTS
)

add_custom_target (
    coverage_html
    DEPENDS
    ${PROJECT_BINARY_DIR}/html/index.html
    COMMENT "Generate coverage report."
)

add_custom_command(
    TARGET coverage_html POST_BUILD
    COMMAND ;
    COMMENT "Open ${PROJECT_BINARY_DIR}/html/index.html in your browser to view the coverage report."
)

macro (pairwise_aligner_test unit_test_cpp)
    file (RELATIVE_PATH unit_test "${CMAKE_SOURCE_DIR}/../unit" "${CMAKE_CURRENT_LIST_DIR}/${unit_test_cpp}")
    seqan3_test_component (target "${unit_test}" TARGET_NAME)
    seqan3_test_component (test_name "${unit_test}" TEST_NAME)

    add_executable (${target} ${unit_test_cpp})
    target_link_libraries (${target} seqan::pairwise_aligner::test::coverage)
    add_test (NAME "${test_name}" COMMAND ${target})

    # any change of a target will invalidate the coverage result;
    # NOTE that this is a GLOBAL variable, because a normal
    # `set(GLOBAL_TEST_COVERAGE_ALL_TESTS)` would not propagate the result when
    # CMakeLists.txt goes out of scope due to a `add_subdirectory`
    set_property(GLOBAL APPEND PROPERTY GLOBAL_TEST_COVERAGE_ALL_TESTS ${target})

    unset (unit_test)
    unset (target)
    unset (test_name)
endmacro ()

seqan3_require_ccache ()
seqan3_require_test ()

# add all unit tests
add_subdirectories_of ("${CMAKE_CURRENT_SOURCE_DIR}/../unit")

# add collected test cases as dependency
get_property(TEST_COVERAGE_ALL_TESTS GLOBAL PROPERTY GLOBAL_TEST_COVERAGE_ALL_TESTS)
add_custom_command (
    OUTPUT ${PROJECT_BINARY_DIR}/pairwise_aligner_coverage.xml
    DEPENDS ${TEST_COVERAGE_ALL_TESTS}
    APPEND
)
