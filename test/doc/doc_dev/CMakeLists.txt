# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

message (STATUS "Configuring devel doc.")

set (PAIRWISE_ALIGNER_DOXYGEN_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}")
set (PAIRWISE_ALIGNER_DOXYGEN_EXCLUDE_SYMBOLS "")
set (PAIRWISE_ALIGNER_DOXYGEN_PREDEFINED_NDEBUG "")
set (PAIRWISE_ALIGNER_DOXYGEN_ENABLED_SECTIONS "DEV")
set (PAIRWISE_ALIGNER_DOXYGEN_EXTRACT_PRIVATE "YES")

configure_file (${PAIRWISE_ALIGNER_DOXYFILE_IN} ${PAIRWISE_ALIGNER_DOXYGEN_OUTPUT_DIR}/Doxyfile)

add_custom_target(doc_dev ALL
                  COMMAND ${DOXYGEN_EXECUTABLE}
                  WORKING_DIRECTORY ${PAIRWISE_ALIGNER_DOXYGEN_OUTPUT_DIR}
                  BYPRODUCTS html/
                  DEPENDS download-cppreference-doxygen-web-tag
                  COMMENT "Generating developer API documentation with Doxygen"
                  VERBATIM)

# Install doc_dev documentation in ./install_doc_dev folder
# cmake --install . --prefix install_doc_dev --component doc_dev
install (
    DIRECTORY "${PAIRWISE_ALIGNER_DOXYGEN_OUTPUT_DIR}/html"
    DESTINATION "${CMAKE_INSTALL_DOCDIR}"
    COMPONENT doc-dev
    EXCLUDE_FROM_ALL
    PATTERN "*.md5" EXCLUDE
    PATTERN "*.map" EXCLUDE
    PATTERN "formula.repository" EXCLUDE
)

### Enable testing

enable_testing()

message (STATUS "Add dev doc test.")
add_test(NAME doc_dev_no_warnings_test
         COMMAND bash -c "${PAIRWISE_ALIGNER_TEST_DOXYGEN_FAIL_ON_WARNINGS}"
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME doc_dev_no_uncond_requires_test
         COMMAND bash -c "${PAIRWISE_ALIGNER_TEST_DOXYGEN_FAIL_ON_UNCOND_REQUIRES}"
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
