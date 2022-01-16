// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise::simd_mask.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <seqan3/std/concepts>
#include <cstdint>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

namespace detail
{
#if defined(__AVX512F__)
inline constexpr size_t max_simd_size = 64;
#elif defined(__AVX2__)
inline constexpr size_t max_simd_size = 32;
#else
inline constexpr size_t max_simd_size = 16;
#endif

} // namespace detail

template <std::integral score_t, size_t simd_size = detail::max_simd_size / sizeof(score_t)>
class alignas(detail::max_simd_size) simd_score;

template <std::unsigned_integral score_t, size_t simd_size = detail::max_simd_size / sizeof(score_t)>
class alignas(detail::max_simd_size) simd_mask;

} // v1
} // namespace seqan::pairwise_aligner
