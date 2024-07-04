// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides concept defintions for simd_type.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <concepts>

#include <pairwise_aligner/simd/math.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace simd {

// ----------------------------------------------------------------------------
// Concept simd_type
// ----------------------------------------------------------------------------

template <typename simd_t>
concept simd_type = requires (simd_t const & a, simd_t const & b, typename simd_t::value_type c) {
    typename simd_t::simd_type;
    typename simd_t::mask_type;

    requires std::integral<decltype(simd_t::size)>;

    { a + b } -> std::same_as<simd_t>;
    { a + c } -> std::same_as<simd_t>;
    { a * b } -> std::same_as<simd_t>;
    { a * c } -> std::same_as<simd_t>;
    { a - b } -> std::same_as<simd_t>;
    { a.le(b) } -> std::same_as<typename simd_t::mask_type>;
    { a.lt(b) } -> std::same_as<typename simd_t::mask_type>;
};

// ----------------------------------------------------------------------------
// Concept saturated_simd_type
// ----------------------------------------------------------------------------

template <typename simd_t>
concept saturated_simd_type = simd_type<simd_t> && requires (simd_t const & a, simd_t const & b) {
    { simd::add_saturated(a, b) } -> std::same_as<simd_t>;
    { simd::subtract_saturated(a, b) } -> std::same_as<simd_t>;
};

} // namespace simd
} // inline namespace v1
} // namespace seqan::pairwise_aligner
