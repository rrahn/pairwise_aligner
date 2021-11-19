// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <cstdint>

/*!\file
 * \brief Provides pairwise_aligner version macros and global variables.
 * \author name <name [AT] domain>
 */

//!\brief The major version as MACRO.
#define PAIRWISE_ALIGNER_VERSION_MAJOR 0
//!\brief The minor version as MACRO.
#define PAIRWISE_ALIGNER_VERSION_MINOR 0
//!\brief The patch version as MACRO.
#define PAIRWISE_ALIGNER_VERSION_PATCH 1

//!\brief The full version as MACRO (number).
#define PAIRWISE_ALIGNER_VERSION (PAIRWISE_ALIGNER_VERSION_MAJOR * 10000 \
                     + PAIRWISE_ALIGNER_VERSION_MINOR * 100 \
                     + PAIRWISE_ALIGNER_VERSION_PATCH)

/*!\brief Converts a number to a string. Preprocessor needs this indirection to
 * properly expand the values to strings.
 */
#define PAIRWISE_ALIGNER_VERSION_CSTRING_HELPER_STR(str) #str

//!\brief Converts version numbers to string.
#define PAIRWISE_ALIGNER_VERSION_CSTRING_HELPER_FUNC(MAJOR, MINOR, PATCH) \
        PAIRWISE_ALIGNER_VERSION_CSTRING_HELPER_STR(MAJOR) "."\
        PAIRWISE_ALIGNER_VERSION_CSTRING_HELPER_STR(MINOR) "."\
        PAIRWISE_ALIGNER_VERSION_CSTRING_HELPER_STR(PATCH)

//!\brief The full version as null terminated string.
#define PAIRWISE_ALIGNER_VERSION_CSTRING \
    PAIRWISE_ALIGNER_VERSION_CSTRING_HELPER_FUNC(PAIRWISE_ALIGNER_VERSION_MAJOR, \
                                             PAIRWISE_ALIGNER_VERSION_MINOR, \
                                             PAIRWISE_ALIGNER_VERSION_PATCH)

namespace pairwise_aligner
{
inline namespace v1
{
//!\brief The major version.
constexpr uint8_t pairwise_aligner_version_major = PAIRWISE_ALIGNER_VERSION_MAJOR;
//!\brief The minor version.
constexpr uint8_t pairwise_aligner_version_minor = PAIRWISE_ALIGNER_VERSION_MINOR;
//!\brief The patch version.
constexpr uint8_t pairwise_aligner_version_patch = PAIRWISE_ALIGNER_VERSION_PATCH;

//!\brief The full version as `std::size_t`.
constexpr std::size_t pairwise_aligner_version = PAIRWISE_ALIGNER_VERSION;

//!\brief The full version as null terminated string.
constexpr char const* pairwise_aligner_version_cstring = PAIRWISE_ALIGNER_VERSION_CSTRING;
} // inline namespace v1
} // namespace lib-template

#undef PAIRWISE_ALIGNER_VERSION_CSTRING_HELPER_STR
#undef PAIRWISE_ALIGNER_VERSION_CSTRING_HELPER_FUNC
