# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

# Minimum cmake version
cmake_minimum_required(VERSION 3.7)

### Find pairwise_aligner
include (../../build_system/pairwise_aligner-config-version.cmake)
set (PAIRWISE_ALIGNER_VERSION "${PACKAGE_VERSION}")

project (pairwise_aligner LANGUAGES NONE VERSION "${PAIRWISE_ALIGNER_PROJECT_VERSION}")

if (NOT EXISTS "${PAIRWISE_ALIGNER_INCLUDE_DIR}/pairwise_aligner/version.hpp")
    message (FATAL_ERROR "Could not find pairwise_aligner. Not building documentation.")
endif ()
message (STATUS "Found: ${PAIRWISE_ALIGNER_INCLUDE_DIR} (${PAIRWISE_ALIGNER_VERSION})")

set (PAIRWISE_ALIGNER_DOXYGEN_INPUT_DIR "${CMAKE_SOURCE_DIR}")

include (pairwise_aligner-doxygen.cmake)

enable_testing ()

if (PAIRWISE_ALIGNER_USER_DOC)
    add_subdirectory(doc_usr)
endif ()

if (PAIRWISE_ALIGNER_DEV_DOC)
    add_subdirectory(doc_dev)
endif ()
