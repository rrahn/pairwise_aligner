// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::simd_rank_selector_impl_avx2.
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
    requires (key_t::size == detail::max_simd_size)
struct _simd_rank_selector_impl_avx2
{
    struct type;
};

template <typename key_t>
    requires (key_t::size == detail::max_simd_size)
using simd_rank_selector_impl_avx2 = typename _simd_rank_selector_impl_avx2<key_t>::type;

template <typename key_t>
    requires (key_t::size == detail::max_simd_size)
struct _simd_rank_selector_impl_avx2<key_t>::type
{
protected:
    using scalar_t = typename key_t::value_type;
    using native_key_t = typename key_t::native_simd_type;
    using split_key_t = std::pair<key_t, key_t>;
    using rank_map_t = std::vector<split_key_t, seqan3::aligned_allocator<split_key_t, alignof(key_t)>>;

    template <std::ranges::random_access_range ranks_t>
        requires (std::same_as<std::ranges::range_value_t<ranks_t>, key_t>)
    static rank_map_t initialise_rank_map(ranks_t && ranks) noexcept
    {
        rank_map_t tmp{};
        tmp.resize(std::ranges::size(ranks));

        for (int32_t i = 0; i < std::ranges::ssize(ranks); ++i) {
            // Unpack ranks into low and high ranks.
            tmp[i] = {to_packed(_mm256_permute2x128_si256(to_native(ranks[i]), to_native(ranks[i]), 0b0000'0000)),
                      to_packed(_mm256_permute2x128_si256(to_native(ranks[i]), to_native(ranks[i]), 0b0001'0001))};
        }

        return tmp;
    }

    static key_t select_rank_for(rank_map_t const & rank_map, key_t const & key) noexcept
    {
        auto [offset, key_low, key_high] = to_offset(key);

        key_t tmp{};
        for (scalar_t i = 0; i < std::ranges::ssize(rank_map); ++i) {
            tmp = tmp | blend(offset.eq(key_t{i}), select_rank_for_impl(rank_map[i], key_low, key_high), key_t{});
        }
        return tmp;
    }

private:
    using offset_type = std::tuple<key_t, key_t, key_t>;

    static offset_type to_offset(key_t const & index) noexcept
    {
        constexpr scalar_t bit_index = std::bit_width(detail::max_simd_size) - 1;
        constexpr key_t modulo_mask{static_cast<scalar_t>((1ull << bit_index) - 1)};

        key_t local_index = index & modulo_mask; // reduce to index 0..31

        constexpr key_t mask{static_cast<scalar_t>(0xFF)};
        // First, select all low indices, then select the high indices.
        key_t local_index_lo = blend(local_index.lt(key_t{static_cast<scalar_t>(16)}), local_index, mask);
        key_t local_index_hi = (local_index_lo ^ mask) | local_index;
        return std::tuple{index >> bit_index, local_index_lo, local_index_hi};
    }

    static key_t select_rank_for_impl(split_key_t const & ranks, key_t const & key_low, key_t const & key_high)
        noexcept
    {
        auto [ranks_lo, ranks_hi] = ranks;
        return to_packed(_mm256_or_si256(_mm256_shuffle_epi8(to_native(ranks_lo), to_native(key_low)),
                                         _mm256_shuffle_epi8(to_native(ranks_hi), to_native(key_high))));
    }

    static __m256i const & to_native(key_t const & packed) noexcept
    {
        return reinterpret_cast<__m256i const &>(packed);
    }

    static key_t to_packed(__m256i const & native) noexcept
    {
        return reinterpret_cast<key_t const &>(native);
    }
};
} // inline namespace v1
} // namespace seqan::pairwise_aligner
