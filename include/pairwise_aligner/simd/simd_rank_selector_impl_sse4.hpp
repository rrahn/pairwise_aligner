// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::simd_rank_selector_impl_sse4.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <immintrin.h>
#include <array>
#include <bit>
#include <ranges>
#include <utility>

#include <seqan3/utility/detail/bits_of.hpp>
#include <seqan3/utility/container/aligned_allocator.hpp>
#include <seqan3/utility/type_traits/basic.hpp>

#include <pairwise_aligner/simd/simd_base.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename key_t>
    requires (key_t::size_v == detail::max_simd_size)
struct _simd_rank_selector_impl_sse4
{
    struct type;
};

template <typename key_t>
    requires (key_t::size_v == detail::max_simd_size)
using simd_rank_selector_impl_sse4 = typename _simd_rank_selector_impl_sse4<key_t>::type;

template <typename key_t>
    requires (key_t::size_v == detail::max_simd_size)
struct _simd_rank_selector_impl_sse4<key_t>::type
{
protected:
    using scalar_t = typename key_t::value_type;
    using native_key_t = typename key_t::native_simd_type;
    using rank_map_t = std::vector<key_t, seqan3::aligned_allocator<key_t, alignof(key_t)>>;

    template <std::ranges::random_access_range ranks_t>
        requires (std::ranges::range_value_t<ranks_t>::size_v == key_t::size_v)
    static rank_map_t initialise_rank_map(ranks_t const & ranks) noexcept
    {
        rank_map_t tmp{std::ranges::begin(ranks), std::ranges::end(ranks)};
        return tmp;
    }

    static key_t select_rank_for(rank_map_t const & rank_map, key_t const & key) noexcept
    {
        auto [offset, local_key] = to_offset(key);

        key_t tmp{};
        for (scalar_t i = 0; i < std::ranges::ssize(rank_map); ++i) {
            tmp = tmp | blend(offset.eq(key_t{i}), select_rank_for_impl(rank_map[i], local_key), key_t{});
        }
        return tmp;
    }

private:
    using offset_type = std::pair<key_t, key_t>;

    static offset_type to_offset(key_t const & index) noexcept
    {
        constexpr scalar_t bit_index = std::bit_width(detail::max_simd_size) - 1;
        constexpr key_t modulo_mask{static_cast<scalar_t>((1ull << bit_index) - 1)};

        return std::pair{index >> bit_index, index & modulo_mask};
    }

    static key_t select_rank_for_impl(key_t const & ranks, key_t const & local_key)
        noexcept
    {
        return to_packed(_mm_shuffle_epi8(to_native(ranks), to_native(local_key)));
    }

    static __m128i const & to_native(key_t const & packed) noexcept
    {
        return reinterpret_cast<__m128i const &>(packed);
    }

    static key_t to_packed(__m128i const & native) noexcept
    {
        return reinterpret_cast<key_t const &>(native);
    }
};
} // inline namespace v1
} // namespace seqan::pairwise_aligner
