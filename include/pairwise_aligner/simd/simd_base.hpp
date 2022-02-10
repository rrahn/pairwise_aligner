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

namespace detail {

template <typename simd_score_t>
struct make_unsigned_impl;

template <std::integral simd_score_t>
struct make_unsigned_impl<simd_score_t>
{
    using type = std::make_unsigned_t<simd_score_t>;
};

template <std::integral score_t, size_t simd_size_v>
struct make_unsigned_impl<simd_score<score_t, simd_size_v>>
{
    using type = simd_score<std::make_unsigned_t<score_t>, simd_size_v>;
};

template <typename simd_score_t>
using make_unsigned_t = typename make_unsigned_impl<simd_score_t>::type;

} // namespace detail

} // v1
} // namespace seqan::pairwise_aligner
