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

#include <pairwise_aligner/simd/concept.hpp>
#include <pairwise_aligner/simd/simd_saturated_score.hpp>

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

template <typename simd_offset_t, typename selector_tag_t>
struct simd_selector;

template <size_t operand_count, size_t operand_bit_width, size_t simd_bit_width>
struct selector_tag;

template <typename score_t, size_t simd_size, template <typename > typename ...policies_t>
struct simd_score_base;

} // namespace detail

template <std::integral score_t, size_t simd_size = detail::max_simd_size / sizeof(score_t)>
using simd_score = detail::simd_score_base<score_t, simd_size>;

template <std::integral score_t, size_t simd_size = detail::max_simd_size / sizeof(score_t)>
using simd_score_saturated = detail::simd_score_base<score_t, simd_size, saturated_score>;

template <std::unsigned_integral score_t, size_t simd_size = detail::max_simd_size / sizeof(score_t)>
class alignas(detail::max_simd_size) simd_mask;

namespace detail {

template <typename score_t>
struct make_unsigned_impl;

template <std::integral score_t>
struct make_unsigned_impl<score_t> : public std::make_unsigned<score_t>
{};

template <simd::simd_type simd_score_t>
struct make_unsigned_impl<simd_score_t>
{
    using type = simd_score<std::make_unsigned_t<typename simd_score_t::value_type>, simd_score_t::size>;
};

template <simd::saturated_simd_type simd_score_t>
struct make_unsigned_impl<simd_score_t>
{
    using type = simd_score_saturated<std::make_unsigned_t<typename simd_score_t::value_type>, simd_score_t::size>;
};

template <typename simd_score_t>
using make_unsigned_t = typename make_unsigned_impl<simd_score_t>::type;

template <typename simd_score_t>
struct make_signed_impl;

template <std::integral simd_score_t>
struct make_signed_impl<simd_score_t>
{
    using type = std::make_signed_t<simd_score_t>;
};

template <std::integral score_t, size_t simd_size_v>
struct make_signed_impl<simd_score<score_t, simd_size_v>>
{
    using type = simd_score<std::make_signed_t<score_t>, simd_size_v>;
};

template <typename simd_score_t>
using make_signed_t = typename make_signed_impl<simd_score_t>::type;

} // namespace detail

} // v1
} // namespace seqan::pairwise_aligner
