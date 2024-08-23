cmake_minimum_required (VERSION 3.7)

# ----------------------------------------------------------------------------
# Initialize cache variable for setting the SIMD ISA to be used for testing
# ----------------------------------------------------------------------------

set (SIMD_ISA_LIST "none" "native" "SSE4" "AVX2" "AVX512BW" "AVX512VBMI")
set (PAIRWISE_ALIGNER_SIMD_ISA "native" CACHE STRING "The SIMD-ISA to select for the pairwise aligner tests. Default is <native>.")
set_property(CACHE PAIRWISE_ALIGNER_SIMD_ISA PROPERTY STRINGS ${SIMD_ISA_LIST})

# ----------------------------------------------------------------------------
# Options for CheckCXXSourceCompiles
# ----------------------------------------------------------------------------

include (CheckCXXSourceCompiles)

# deactivate messages in check_*
set (CMAKE_REQUIRED_QUIET 1)
# use global variables in check_* calls
set (CMAKE_REQUIRED_INCLUDES ${CMAKE_INCLUDE_PATH})
set (CMAKE_REQUIRED_FLAGS ${CMAKE_CXX_FLAGS})

# ----------------------------------------------------------------------------
# Select SIMD CPU Flags for the respective SIMD ISAs
# ----------------------------------------------------------------------------

set (PAIRWISE_ALIGNER_SIMD_CXX_FLAGS "")
set (PAIRWISE_ALIGNER_SDE_TESTS "")

# Check if the selected SIMD ISA is supported natively by the current processor architecture.
# Sets the result_var to TRUE if the selected SIMD ISA is supported natively, otherwise FALSE.
function (check_native_simd_support simd_isa_idx result_var)
    set (CXX_CPU_MACROS "defined(__SSE4_1__) && defined(__SSE4_2__)"
                        "defined(__AVX2__)"
                        "defined(__AVX512BW__)"
                        "defined(__AVX512VBMI__)")
    list (GET CXX_CPU_MACROS ${simd_isa_idx} CXX_TEST_MACRO)
    set (CXX_MACRO_TEST_SOURCE
        "#if !(${CXX_TEST_MACRO})
        #error NO ${PAIRWISE_ALIGNER_SIMD_ISA}
        #endif
        int main() {}")
    set (CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -march=native")
    check_cxx_source_compiles("${CXX_MACRO_TEST_SOURCE}" NATIVE_SIMD_ISA_SUPPORTED)
    set (${result_var} ${NATIVE_SIMD_ISA_SUPPORTED} PARENT_SCOPE)
endfunction ()

# Do nothing if no SIMD ISA is selected
if (NOT ${PAIRWISE_ALIGNER_SIMD_ISA} STREQUAL "none")

    # Use native architecture flags if the selected SIMD ISA is `native`
    if (${PAIRWISE_ALIGNER_SIMD_ISA} STREQUAL "native")
        message (STATUS "Using native cpu compile flags.")
        set (PAIRWISE_ALIGNER_SIMD_CXX_FLAGS "-march=native")
    else()
        # Choose explicit SIMD falgs.
        list (SUBLIST SIMD_ISA_LIST 2 -1 SIMD_ISA_LIST_SUBSET)
        list(FIND SIMD_ISA_LIST_SUBSET ${PAIRWISE_ALIGNER_SIMD_ISA} ITEM_INDEX)
        if(${ITEM_INDEX} EQUAL -1)
            message(STATUS "Item '${PAIRWISE_ALIGNER_SIMD_ISA}' not found in the list. Fall back to use no SIMD flags.")
        else()
            message (STATUS "Using explicit SIMD compile flags.")
            set (EXPLICIT_SIMD_CXX_FLAGS "-msse4"
                                         "-mavx2"
                                         "-mavx512f -mavx512dq -mavx512cd -mavx512vl -mavx512bw"
                                         "-mavx512f -mavx512dq -mavx512cd -mavx512vl -mavx512bw -mavx512vbmi")
            list (GET EXPLICIT_SIMD_CXX_FLAGS ${ITEM_INDEX} EXPLICIT_SIMD_CXX_FLAG)
            set (PAIRWISE_ALIGNER_SIMD_CXX_FLAGS ${EXPLICIT_SIMD_CXX_FLAG})

            check_native_simd_support(${ITEM_INDEX} HAS_NATIVE_SIMD_SUPPORT)

            # Selecting SDE cpu type for emulating respective SIMD ISA if the native platform does not support it
            if (NOT HAS_NATIVE_SIMD_SUPPORT)
                set (SDE_CPU_TYPES "snb"
                                   "hsw"
                                   "skx"
                                   "icx")
                list (GET SDE_CPU_TYPES ${ITEM_INDEX} SDE_CPU_TYPE)
                set (PAIRWISE_ALIGNER_SDE_TESTS ${SDE_CPU_TYPE})
                message (STATUS "SIMD flag not natively supported. Using SDE cpu-type ${SDE_CPU_TYPE} for emulation.")
            endif()
        endif ()
    endif ()
endif ()
