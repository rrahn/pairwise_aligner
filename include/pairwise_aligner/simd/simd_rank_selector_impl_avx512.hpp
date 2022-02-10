// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::simd_rank_selector_impl_avx512.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <immintrin.h>
#include <array>
#include <seqan3/std/bit>
#include <seqan3/std/ranges>
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
struct _simd_rank_selector_impl_avx512
{
    struct type;
};

template <typename key_t>
    requires (key_t::size == detail::max_simd_size)
using simd_rank_selector_impl_avx512 = typename _simd_rank_selector_impl_avx512<key_t>::type;

template <typename key_t>
    requires (key_t::size == detail::max_simd_size)
struct _simd_rank_selector_impl_avx512<key_t>::type
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

        constexpr std::array<int64_t, 8> shuffle_mask{0, 4, 1, 5, 2, 6, 3, 7};
        for (int32_t i = 0; i < std::ranges::ssize(ranks); ++i) {
            // 1. reorder ranks according to shuffle_mask because unpack16 works on low/high 64 packed operands of
            // 4x128 bit lanes. After the shuffle the indices are back in the correct order for the selection operation.
            __m512i a = _mm512_permutexvar_epi64(to_native(shuffle_mask), to_native(ranks[i]));

            // second unpack low and high values
            tmp[i] = std::pair{to_packed(_mm512_unpacklo_epi8(__m512i{}, a)),
                               to_packed(_mm512_unpackhi_epi8(__m512i{}, a))};
        }
        return tmp;
    }

    static key_t select_rank_for(rank_map_t const & rank_map, key_t const & key) noexcept
    {
        auto [offset, key_low, key_high] = to_offset(key);

        key_t tmp{};
        for (int32_t i = 0; i < std::ranges::ssize(rank_map); ++i) {
            tmp |= blend(offset.eq(key_t{static_cast<scalar_t>(i)}),
                         select_rank_for_impl(rank_map[i], key_low, key_high), key_t{});
        }
        return tmp;
    }

private:
    using offset_type = std::tuple<key_t, key_t, key_t>;

    static offset_type to_offset(key_t const & index) noexcept
    {
        constexpr scalar_t bit_index = std::bit_width(detail::max_simd_size) - 1;
        constexpr key_t modulo_mask{static_cast<scalar_t>((1ull << bit_index) - 1)};

        key_t local_index = index & modulo_mask; // & 63 <=> mod 64.
        // Rotate high keys to the right by 8 bit to put them on the low 8 bits of the 16-bit shuffle mask.
        return std::tuple{index >> bit_index, local_index, to_packed(rotate_right(to_native(local_index)))};
    }

    static __m512i rotate_right(__m512i const & index) noexcept
    {
        return _mm512_ror_epi64(index, seqan3::detail::bits_of<scalar_t>);
    }

    static key_t select_rank_for_impl(split_key_t const & ranks, key_t const & key_low, key_t const & key_high)
        noexcept
    {
        auto [ranks_lo, ranks_hi] = ranks;
        // Load 32x16 bit operands using the low 8 bits of each 16 bit operand (they are rotated right by 8 bits).
        __m512i lo_32x16 = _mm512_permutex2var_epi16(to_native(ranks_lo), to_native(key_low), to_native(ranks_hi));
        // Load 32x16 bit operands using the high 8 bits of each 16 bit operand (unmodified).
        __m512i hi_32x16 = _mm512_permutex2var_epi16(to_native(ranks_lo), to_native(key_high), to_native(ranks_hi));
        // Or the results after rotating low 8 bits of each 16-bit operand of the lo_32x16 register to the right by 8
        // bits, which corresponds to their corres position the final 8-bit packed result vector.
        return to_packed(rotate_right(lo_32x16) | hi_32x16);
    }

    template <typename packed_t>
    static __m512i const & to_native(packed_t const & packed) noexcept
    {
        return reinterpret_cast<__m512i const &>(packed);
    }

    static key_t to_packed(__m512i const & native) noexcept
    {
        return reinterpret_cast<key_t const &>(native);
    }
};
} // inline namespace v1
} // namespace seqan::pairwise_aligner
