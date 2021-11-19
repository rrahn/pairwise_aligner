# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/rrahn/lib_template/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------
#
# This CMake module will try to find lib_template and its dependencies.  You can use
# it the same way you would use any other CMake module.
#
#   find_package (lib_template [REQUIRED] ...)
#
# lib_template has the following platform requirements:
# <platform requirements>
#
# lib_template requires the following libraries:
#
# <mandatory library dependencies>
#
# lib_template has the following optional dependencies:
#
# <optional library dependencies>
#
# Once the search has been performed, the following variables will be set.
#
#   LIB_TEMPLATE_FOUND            -- Indicate whether lib_template was found and requirements met.
#
#   LIB_TEMPLATE_VERSION          -- The version as string, e.g. "3.0.0"
#   LIB_TEMPLATE_VERSION_MAJOR    -- e.g. 3
#   LIB_TEMPLATE_VERSION_MINOR    -- e.g. 0
#   LIB_TEMPLATE_VERSION_PATCH    -- e.g. 0
#
#   LIB_TEMPLATE_INCLUDE_DIRS     -- to be passed to include_directories ()
#   LIB_TEMPLATE_LIBRARIES        -- to be passed to target_link_libraries ()
#   LIB_TEMPLATE_DEFINITIONS      -- to be passed to add_definitions ()
#   LIB_TEMPLATE_CXX_FLAGS        -- to be added to CMAKE_CXX_FLAGS
#
# Additionally, the following [IMPORTED][IMPORTED] targets are defined:
#
#   seqan::lib_template          -- interface target where
#                                  target_link_libraries(target seqan::lib_template)
#                              automatically sets
#                                  target_include_directories(target $LIB_TEMPLATE_INCLUDE_DIRS),
#                                  target_link_libraries(target $LIB_TEMPLATE_LIBRARIES),
#                                  target_compile_definitions(target $LIB_TEMPLATE_DEFINITIONS) and
#                                  target_compile_options(target $LIB_TEMPLATE_CXX_FLAGS)
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
    message (STATUS "${ColourBold}Finding lib_template and checking requirements:${ColourReset}")
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
# Find lib_template include path
# ----------------------------------------------------------------------------

# Note that lib_template-config.cmake can be standalone and thus LIB_TEMPLATE_CLONE_DIR might be empty.
# * `LIB_TEMPLATE_CLONE_DIR` was already found in lib_template-config-version.cmake
# * `LIB_TEMPLATE_INCLUDE_DIR` was already found in lib_template-config-version.cmake
find_path (LIB_TEMPLATE_SUBMODULES_DIR NAMES lib/seqan3 HINTS "${LIB_TEMPLATE_CLONE_DIR}" "${LIB_TEMPLATE_INCLUDE_DIR}/lib_template")

if (LIB_TEMPLATE_INCLUDE_DIR)
    config_print ("lib_template include dir found:   ${LIB_TEMPLATE_INCLUDE_DIR}")
else ()
    config_error ("lib_template include directory could not be found (LIB_TEMPLATE_INCLUDE_DIR: '${LIB_TEMPLATE_INCLUDE_DIR}')")
endif ()

# ----------------------------------------------------------------------------
# Detect if we are a clone of repository and if yes auto-add submodules
# ----------------------------------------------------------------------------

if (LIB_TEMPLATE_CLONE_DIR)
    config_print ("Detected as running from a repository checkout…")
endif ()

if (LIB_TEMPLATE_SUBMODULES_DIR)
    file (GLOB submodules ${LIB_TEMPLATE_SUBMODULES_DIR}/lib/*/include)
    foreach (submodule ${submodules})
        if (IS_DIRECTORY ${submodule})
            config_print ("  …adding submodule include:  ${submodule}")
            set (LIB_TEMPLATE_DEPENDENCY_INCLUDE_DIRS ${submodule} ${LIB_TEMPLATE_DEPENDENCY_INCLUDE_DIRS})
        endif ()
    endforeach ()
endif ()

# ----------------------------------------------------------------------------
# Options for CheckCXXSourceCompiles
# ----------------------------------------------------------------------------

# deactivate messages in check_*
set (CMAKE_REQUIRED_QUIET       1)
# use global variables in Check* calls
set (CMAKE_REQUIRED_INCLUDES    ${CMAKE_INCLUDE_PATH} ${LIB_TEMPLATE_INCLUDE_DIR} ${LIB_TEMPLATE_DEPENDENCY_INCLUDE_DIRS})
set (CMAKE_REQUIRED_FLAGS       ${CMAKE_CXX_FLAGS})

# ----------------------------------------------------------------------------
# thread support (pthread, windows threads)
# ----------------------------------------------------------------------------

find_package (SeqAn3 REQUIRED HINTS ${CMAKE_CURRENT_LIST_DIR}/../lib/seqan3/build_system)

if (SeqAn3_FOUND)
    set (LIB_TEMPLATE_LIBRARIES ${LIB_TEMPLATE_LIBRARIES} seqan3::seqan3)
    config_print ("Required dependency:                 SeqAn3 found.")
else ()
    config_error ("The required SeqAn3 library was marked as required, but wasn't found.")
endif ()

# ----------------------------------------------------------------------------
# Export targets
# ----------------------------------------------------------------------------

separate_arguments (LIB_TEMPLATE_CXX_FLAGS_LIST UNIX_COMMAND "${LIB_TEMPLATE_CXX_FLAGS}")

add_library (lib_template INTERFACE)
target_compile_definitions (lib_template INTERFACE ${LIB_TEMPLATE_DEFINITIONS})
target_compile_options (lib_template INTERFACE ${LIB_TEMPLATE_CXX_FLAGS_LIST})
target_link_libraries (lib_template INTERFACE "${LIB_TEMPLATE_LIBRARIES}")
# include seqan3/include/ as -I, because seqan3 should never produce warnings.
target_include_directories (lib_template INTERFACE "${LIB_TEMPLATE_INCLUDE_DIR}")
# include everything except seqan3/include/ as -isystem, i.e.
# a system header which suppresses warnings of external libraries.
target_include_directories (lib_template SYSTEM INTERFACE "${LIB_TEMPLATE_DEPENDENCY_INCLUDE_DIRS}")
add_library (seqan::lib_template ALIAS lib_template)

# propagate LIB_TEMPLATE_INCLUDE_DIR into LIB_TEMPLATE_INCLUDE_DIRS
set (LIB_TEMPLATE_INCLUDE_DIRS ${LIB_TEMPLATE_INCLUDE_DIR} ${LIB_TEMPLATE_DEPENDENCY_INCLUDE_DIRS})

# ----------------------------------------------------------------------------
# Finish find_package call
# ----------------------------------------------------------------------------

find_package_handle_standard_args (${CMAKE_FIND_PACKAGE_NAME} REQUIRED_VARS LIB_TEMPLATE_INCLUDE_DIR)

# Set LIB_TEMPLATE_* variables with the content of ${CMAKE_FIND_PACKAGE_NAME}_(FOUND|...|VERSION)
# This needs to be done, because `find_package(lib_template)` might be called in any case-sensitive way and we want to
# guarantee that LIB_TEMPLATE_* are always set.
foreach (package_var FOUND DIR ROOT CONFIG VERSION VERSION_MAJOR VERSION_MINOR VERSION_PATCH VERSION_TWEAK VERSION_COUNT)
    set (LIB_TEMPLATE_${package_var} "${${CMAKE_FIND_PACKAGE_NAME}_${package_var}}")
endforeach ()

set (CMAKE_REQUIRED_QUIET ${CMAKE_REQUIRED_QUIET_SAVE})

if (LIB_TEMPLATE_FIND_DEBUG)
  message ("Result for ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt")
  message ("")
  message ("  CMAKE_BUILD_TYPE            ${CMAKE_BUILD_TYPE}")
  message ("  CMAKE_SOURCE_DIR            ${CMAKE_SOURCE_DIR}")
  message ("  CMAKE_INCLUDE_PATH          ${CMAKE_INCLUDE_PATH}")
  message ("  LIB_TEMPLATE_INCLUDE_DIR          ${LIB_TEMPLATE_INCLUDE_DIR}")
  message ("")
  message ("  ${CMAKE_FIND_PACKAGE_NAME}_FOUND                ${${CMAKE_FIND_PACKAGE_NAME}_FOUND}")
  message ("  LIB_TEMPLATE_HAS_ZLIB             ${ZLIB_FOUND}")
  message ("  LIB_TEMPLATE_HAS_BZIP2            ${BZIP2_FOUND}")
  message ("")
  message ("  LIB_TEMPLATE_INCLUDE_DIRS         ${LIB_TEMPLATE_INCLUDE_DIRS}")
  message ("  LIB_TEMPLATE_LIBRARIES            ${LIB_TEMPLATE_LIBRARIES}")
  message ("  LIB_TEMPLATE_DEFINITIONS          ${LIB_TEMPLATE_DEFINITIONS}")
  message ("  LIB_TEMPLATE_CXX_FLAGS            ${LIB_TEMPLATE_CXX_FLAGS}")
  message ("")
  message ("  LIB_TEMPLATE_VERSION              ${LIB_TEMPLATE_VERSION}")
  message ("  LIB_TEMPLATE_VERSION_MAJOR        ${LIB_TEMPLATE_VERSION_MAJOR}")
  message ("  LIB_TEMPLATE_VERSION_MINORG       ${LIB_TEMPLATE_VERSION_MINOR}")
  message ("  LIB_TEMPLATE_VERSION_PATCH        ${LIB_TEMPLATE_VERSION_PATCH}")
endif ()
