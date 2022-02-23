// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::simd::detail::host_simd_tag.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <cinttypes>

#include <seqan3/utility/simd/concept.hpp>
#include <seqan3/utility/simd/simd_traits.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace simd::detail {

template <std::size_t _size, std::size_t _width>
struct host_simd_tag
{
    static constexpr std::size_t size = _size;
    static constexpr std::size_t width = _width;
};

template <seqan3::simd::simd_concept host_simd_t>
using host_simd_tag_t = host_simd_tag<seqan3::simd::simd_traits<host_simd_t>::length,
                                      seqan3::simd::simd_traits<host_simd_t>::max_length>;

} // namespace simd::detail
} // inline namespace v1
} // namespace seqan::pairwise_aligner
