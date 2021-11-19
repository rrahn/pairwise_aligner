# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/rrahn/lib_template/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

# Minimum cmake version
cmake_minimum_required(VERSION 3.7)

### Find lib_template
include (../../build_system/lib_template-config-version.cmake)
set (LIB_TEMPLATE_VERSION "${PACKAGE_VERSION}")

project (lib_template LANGUAGES NONE VERSION "${LIB_TEMPLATE_PROJECT_VERSION}")

if (NOT EXISTS "${LIB_TEMPLATE_INCLUDE_DIR}/lib_template/version.hpp")
    message (FATAL_ERROR "Could not find lib_template. Not building documentation.")
endif ()
message (STATUS "Found: ${LIB_TEMPLATE_INCLUDE_DIR} (${LIB_TEMPLATE_VERSION})")

set (LIB_TEMPLATE_DOXYGEN_INPUT_DIR "${CMAKE_SOURCE_DIR}")

include (lib_template-doxygen.cmake)

enable_testing ()

if (LIB_TEMPLATE_USER_DOC)
    add_subdirectory(doc_usr)
endif ()

if (LIB_TEMPLATE_DEV_DOC)
    add_subdirectory(doc_dev)
endif ()
