# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------
#
# This CMake module will try to find pairwise_aligner and its dependencies.  You can use
# it the same way you would use any other CMake module.
#
#   find_package (pairwise_aligner [REQUIRED] ...)
#
# pairwise_aligner has the following platform requirements:
# <platform requirements>
#
# pairwise_aligner requires the following libraries:
#
# <mandatory library dependencies>
#
# pairwise_aligner has the following optional dependencies:
#
# <optional library dependencies>
#
# Once the search has been performed, the following variables will be set.
#
#   PAIRWISE_ALIGNER_FOUND            -- Indicate whether pairwise_aligner was found and requirements met.
#
#   PAIRWISE_ALIGNER_VERSION          -- The version as string, e.g. "3.0.0"
#   PAIRWISE_ALIGNER_VERSION_MAJOR    -- e.g. 3
#   PAIRWISE_ALIGNER_VERSION_MINOR    -- e.g. 0
#   PAIRWISE_ALIGNER_VERSION_PATCH    -- e.g. 0
#
#   PAIRWISE_ALIGNER_INCLUDE_DIRS     -- to be passed to include_directories ()
#   PAIRWISE_ALIGNER_LIBRARIES        -- to be passed to target_link_libraries ()
#   PAIRWISE_ALIGNER_DEFINITIONS      -- to be passed to add_definitions ()
#   PAIRWISE_ALIGNER_CXX_FLAGS        -- to be added to CMAKE_CXX_FLAGS
#
# Additionally, the following [IMPORTED][IMPORTED] targets are defined:
#
#   seqan::pairwise_aligner          -- interface target where
#                                  target_link_libraries(target seqan::pairwise_aligner)
#                              automatically sets
#                                  target_include_directories(target $PAIRWISE_ALIGNER_INCLUDE_DIRS),
#                                  target_link_libraries(target $PAIRWISE_ALIGNER_LIBRARIES),
#                                  target_compile_definitions(target $PAIRWISE_ALIGNER_DEFINITIONS) and
#                                  target_compile_options(target $PAIRWISE_ALIGNER_CXX_FLAGS)
#                              for a target.
#
#   [IMPORTED]: https://cmake.org/cmake/help/v3.10/prop_tgt/IMPORTED.html#prop_tgt:IMPORTED
#
# ============================================================================

cmake_minimum_required (VERSION 3.4...3.12)

# ----------------------------------------------------------------------------
# Set initial variables
# ----------------------------------------------------------------------------

# make output globally quiet if required by find_package, this effects cmake functions like `check_*`
set(CMAKE_REQUIRED_QUIET_SAVE ${CMAKE_REQUIRED_QUIET})
set(CMAKE_REQUIRED_QUIET ${${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY})

# ----------------------------------------------------------------------------
# Greeter
# ----------------------------------------------------------------------------

string (ASCII 27 Esc)
set (ColourBold "${Esc}[1m")
set (ColourReset "${Esc}[m")

if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
    message (STATUS "${ColourBold}Finding pairwise_aligner and checking requirements:${ColourReset}")
endif ()

# ----------------------------------------------------------------------------
# Includes
# ----------------------------------------------------------------------------

include (CheckIncludeFileCXX)
include (CheckCXXSourceCompiles)
include (FindPackageHandleStandardArgs)

# ----------------------------------------------------------------------------
# Pretty printing and error handling
# ----------------------------------------------------------------------------

macro (config_print text)
    if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
        message (STATUS "  ${text}")
    endif ()
endmacro ()

macro (config_error text)
    if (${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED)
        message (FATAL_ERROR ${text})
    else ()
        if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
            message (WARNING ${text})
        endif ()
        return ()
    endif ()
endmacro ()

# ----------------------------------------------------------------------------
# Find pairwise_aligner include path
# ----------------------------------------------------------------------------

# Note that pairwise_aligner-config.cmake can be standalone and thus PAIRWISE_ALIGNER_CLONE_DIR might be empty.
# * `PAIRWISE_ALIGNER_CLONE_DIR` was already found in pairwise_aligner-config-version.cmake
# * `PAIRWISE_ALIGNER_INCLUDE_DIR` was already found in pairwise_aligner-config-version.cmake
find_path (PAIRWISE_ALIGNER_SUBMODULES_DIR NAMES lib/seqan3 HINTS "${PAIRWISE_ALIGNER_CLONE_DIR}" "${PAIRWISE_ALIGNER_INCLUDE_DIR}/pairwise_aligner")

if (PAIRWISE_ALIGNER_INCLUDE_DIR)
    config_print ("pairwise_aligner include dir found:   ${PAIRWISE_ALIGNER_INCLUDE_DIR}")
else ()
    config_error ("pairwise_aligner include directory could not be found (PAIRWISE_ALIGNER_INCLUDE_DIR: '${PAIRWISE_ALIGNER_INCLUDE_DIR}')")
endif ()

# ----------------------------------------------------------------------------
# Detect if we are a clone of repository and if yes auto-add submodules
# ----------------------------------------------------------------------------

if (PAIRWISE_ALIGNER_CLONE_DIR)
    config_print ("Detected as running from a repository checkout…")
endif ()

if (PAIRWISE_ALIGNER_SUBMODULES_DIR)
    file (GLOB submodules ${PAIRWISE_ALIGNER_SUBMODULES_DIR}/lib/*/include)
    foreach (submodule ${submodules})
        if (IS_DIRECTORY ${submodule})
            config_print ("  …adding submodule include:  ${submodule}")
            set (PAIRWISE_ALIGNER_DEPENDENCY_INCLUDE_DIRS ${submodule} ${PAIRWISE_ALIGNER_DEPENDENCY_INCLUDE_DIRS})
        endif ()
    endforeach ()
endif ()

# ----------------------------------------------------------------------------
# Options for CheckCXXSourceCompiles
# ----------------------------------------------------------------------------

# deactivate messages in check_*
set (CMAKE_REQUIRED_QUIET       1)
# use global variables in Check* calls
set (CMAKE_REQUIRED_INCLUDES    ${CMAKE_INCLUDE_PATH} ${PAIRWISE_ALIGNER_INCLUDE_DIR} ${PAIRWISE_ALIGNER_DEPENDENCY_INCLUDE_DIRS})
set (CMAKE_REQUIRED_FLAGS       ${CMAKE_CXX_FLAGS})

# ----------------------------------------------------------------------------
# thread support (pthread, windows threads)
# ----------------------------------------------------------------------------

find_package (SeqAn3 REQUIRED HINTS ${CMAKE_CURRENT_LIST_DIR}/../lib/seqan3/build_system)

if (SeqAn3_FOUND)
    set (PAIRWISE_ALIGNER_LIBRARIES ${PAIRWISE_ALIGNER_LIBRARIES} seqan3::seqan3)
    config_print ("Required dependency:                 SeqAn3 found.")
else ()
    config_error ("The required SeqAn3 library was marked as required, but wasn't found.")
endif ()

# ----------------------------------------------------------------------------
# Export targets
# ----------------------------------------------------------------------------

separate_arguments (PAIRWISE_ALIGNER_CXX_FLAGS_LIST UNIX_COMMAND "${PAIRWISE_ALIGNER_CXX_FLAGS}")

add_library (pairwise_aligner INTERFACE)
target_compile_definitions (pairwise_aligner INTERFACE ${PAIRWISE_ALIGNER_DEFINITIONS})
target_compile_options (pairwise_aligner INTERFACE ${PAIRWISE_ALIGNER_CXX_FLAGS_LIST})
target_link_libraries (pairwise_aligner INTERFACE "${PAIRWISE_ALIGNER_LIBRARIES}")
# include seqan3/include/ as -I, because seqan3 should never produce warnings.
target_include_directories (pairwise_aligner INTERFACE "${PAIRWISE_ALIGNER_INCLUDE_DIR}")
# include everything except seqan3/include/ as -isystem, i.e.
# a system header which suppresses warnings of external libraries.
target_include_directories (pairwise_aligner SYSTEM INTERFACE "${PAIRWISE_ALIGNER_DEPENDENCY_INCLUDE_DIRS}")
add_library (seqan::pairwise_aligner ALIAS pairwise_aligner)

# propagate PAIRWISE_ALIGNER_INCLUDE_DIR into PAIRWISE_ALIGNER_INCLUDE_DIRS
set (PAIRWISE_ALIGNER_INCLUDE_DIRS ${PAIRWISE_ALIGNER_INCLUDE_DIR} ${PAIRWISE_ALIGNER_DEPENDENCY_INCLUDE_DIRS})

# ----------------------------------------------------------------------------
# Finish find_package call
# ----------------------------------------------------------------------------

find_package_handle_standard_args (${CMAKE_FIND_PACKAGE_NAME} REQUIRED_VARS PAIRWISE_ALIGNER_INCLUDE_DIR)

# Set PAIRWISE_ALIGNER_* variables with the content of ${CMAKE_FIND_PACKAGE_NAME}_(FOUND|...|VERSION)
# This needs to be done, because `find_package(pairwise_aligner)` might be called in any case-sensitive way and we want to
# guarantee that PAIRWISE_ALIGNER_* are always set.
foreach (package_var FOUND DIR ROOT CONFIG VERSION VERSION_MAJOR VERSION_MINOR VERSION_PATCH VERSION_TWEAK VERSION_COUNT)
    set (PAIRWISE_ALIGNER_${package_var} "${${CMAKE_FIND_PACKAGE_NAME}_${package_var}}")
endforeach ()

set (CMAKE_REQUIRED_QUIET ${CMAKE_REQUIRED_QUIET_SAVE})

if (PAIRWISE_ALIGNER_FIND_DEBUG)
  message ("Result for ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt")
  message ("")
  message ("  CMAKE_BUILD_TYPE            ${CMAKE_BUILD_TYPE}")
  message ("  CMAKE_SOURCE_DIR            ${CMAKE_SOURCE_DIR}")
  message ("  CMAKE_INCLUDE_PATH          ${CMAKE_INCLUDE_PATH}")
  message ("  PAIRWISE_ALIGNER_INCLUDE_DIR          ${PAIRWISE_ALIGNER_INCLUDE_DIR}")
  message ("")
  message ("  ${CMAKE_FIND_PACKAGE_NAME}_FOUND                ${${CMAKE_FIND_PACKAGE_NAME}_FOUND}")
  message ("  PAIRWISE_ALIGNER_HAS_ZLIB             ${ZLIB_FOUND}")
  message ("  PAIRWISE_ALIGNER_HAS_BZIP2            ${BZIP2_FOUND}")
  message ("")
  message ("  PAIRWISE_ALIGNER_INCLUDE_DIRS         ${PAIRWISE_ALIGNER_INCLUDE_DIRS}")
  message ("  PAIRWISE_ALIGNER_LIBRARIES            ${PAIRWISE_ALIGNER_LIBRARIES}")
  message ("  PAIRWISE_ALIGNER_DEFINITIONS          ${PAIRWISE_ALIGNER_DEFINITIONS}")
  message ("  PAIRWISE_ALIGNER_CXX_FLAGS            ${PAIRWISE_ALIGNER_CXX_FLAGS}")
  message ("")
  message ("  PAIRWISE_ALIGNER_VERSION              ${PAIRWISE_ALIGNER_VERSION}")
  message ("  PAIRWISE_ALIGNER_VERSION_MAJOR        ${PAIRWISE_ALIGNER_VERSION_MAJOR}")
  message ("  PAIRWISE_ALIGNER_VERSION_MINORG       ${PAIRWISE_ALIGNER_VERSION_MINOR}")
  message ("  PAIRWISE_ALIGNER_VERSION_PATCH        ${PAIRWISE_ALIGNER_VERSION_PATCH}")
endif ()
