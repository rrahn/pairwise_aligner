// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implementation for avx2 selection.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <immintrin.h>

#include <array>
#include <concepts>
#include <type_traits>

#include <pairwise_aligner/simd/simd_base.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

namespace detail {

// ----------------------------------------------------------------------------
// epi8
// ----------------------------------------------------------------------------

template <typename simd_offset_t, size_t operand_count>
    requires (operand_count <= 16)
struct simd_selector<simd_offset_t, selector_tag<operand_count, 8, 256>>
{
    static constexpr bool in_lane_shuffle = true;
    static constexpr size_t max_operand_count = 16;

    template <typename value_t>
    using address_t = std::array<value_t, 1>;

    simd_offset_t const & offsets{};

    template <typename value_t>
    constexpr __m256i operator()(address_t<value_t> const & address) const noexcept {
        return _mm256_shuffle_epi8(reinterpret_cast<__m256i const &>(address[0]),
                                   reinterpret_cast<__m256i const &>(offsets));
    }
};

template <typename simd_offset_t, size_t operand_count>
    requires (operand_count > 16)
struct simd_selector<simd_offset_t, selector_tag<operand_count, 8, 256>>
{
    static constexpr bool in_lane_shuffle = false;
    static constexpr size_t max_operand_count = 32;

    template <typename value_t>
    using address_t = std::array<value_t, 1>;

    using scalar_offset_t = typename simd_offset_t::value_type;

    __m256i keys_low;
    __m256i keys_high;

    simd_selector() = delete;

    /**
     * @brief Construct a new simd selector object
     *
     * @param[in] offsets The offsets to select the values from the input array from.
     *
     * The select operation utilises the AVX2 shuffle instruction for 8-bit operand types to select the values from the
     * input vector. Hence the selector can address at most 32 values in a single operation.
     * Furthermore, the shuffle instruction only selects within a 128-bit lane requiring to split the selection
     * keys into two parts for the two address spaces 0-15 and 16-31. These are represented by `keys_low` and
     * `keys_high` respectively.
     *
     * To initialize these vectors the offsets are first truncated to the lower 5 bits (index range 0-31) and then
     * the MSB is set to 1 for keys greater than 15. Using this mask the keys for the low and high address space are
     * marked accordingly, i.e. in keys_low all entries with an index >= 16 have the MSB set to 1 and in keys_high
     * all entries with an index < 16 have the MSB set to 1.
     * The MSB is later used to set the corresponding entry to 0 (see documentation of `_mm256_shuffle_epi8` for
     * details).
     */
    simd_selector(simd_offset_t const & offsets) noexcept
    {
        __m256i keys = reinterpret_cast<__m256i const &>(offsets);
        // Step 1: truncate the offsets to the lower 5 bits (index range 0-31).
        keys = _mm256_and_si256(keys, _mm256_set1_epi8(uint8_t(31)));
        // Step 2: set the MSB to 1 for keys greater than 15.
        __m256i high_keys_mask = _mm256_cmpgt_epi8(keys, _mm256_set1_epi8(uint8_t(15)));
        high_keys_mask = _mm256_and_si256(high_keys_mask, _mm256_set1_epi8(uint8_t(128)));
        // Step 3: set keys_low to keys where the MSB is set if the corresponding entry is from the high address space.
        keys_low = _mm256_or_si256(keys, high_keys_mask);
        // Step 4: set keys_high to keys where the MSB is set if the corresponding entry is from the low address space.
        keys_high = _mm256_or_si256(keys, _mm256_xor_si256(high_keys_mask,
                                                           _mm256_set1_epi8(static_cast<uint8_t>(128))));

        // std::cout << "high keys mask = "; reinterpret_cast<simd_offset_t const &>(high_keys_mask).print(std::cout) << "\n";
        // std::cout << "keys low = "; reinterpret_cast<simd_offset_t const &>(keys_low).print(std::cout) << "\n";
        // std::cout << "keys high = "; reinterpret_cast<simd_offset_t const &>(keys_high).print(std::cout) << "\n";
    }

    /**
     * @brief Selects the values from the input array using the offsets from the initialization.
     *
     * @tparam value_t The value type for address array.
     * @param value_array An array containing the values to select from.
     * @return constexpr __m256i The selected values in a 256-bit register.
     *
     * To select the values from the input array the `_mm256_shuffle_epi8` shuffle instruction is used.
     * However, this instruction can only select values from within a 128-bit lane. Hence, the given values are first
     * split into two vectors where the first vector stores in both 128-bit lanes the values from the lower 128-bit lane
     * of the original input values and the second vector stores in both 128-bit lanes the values from the higher
     * 128-bit lane of the original input values. This is achieved by two calls of `_mm256_permute4x64_epi64`.
     * Having this, the shuffle instruction can be applied to both vectors separately using the corresponding keys for
     * the low and high address space.
     * The result of the two shuffle instructions is then combined using a bitwise or operation.
     */
    template <typename value_t>
    constexpr __m256i operator()(address_t<value_t> const & value_array) const noexcept {

        __m256i values = reinterpret_cast<__m256i const &>(value_array[0]);
        __m256i values_low = _mm256_permute4x64_epi64(values, 0b0100'0100);
        __m256i values_high = _mm256_permute4x64_epi64(values, 0b1110'1110);


        __m256i tmp = _mm256_or_si256(_mm256_shuffle_epi8(values_low, keys_low),
                                       _mm256_shuffle_epi8(values_high, keys_high));
        // std::cout << "values_low = "; reinterpret_cast<value_t const &>(values_low).print(std::cout) << "\n";
        // std::cout << "values_high = "; reinterpret_cast<value_t const &>(values_high).print(std::cout) << "\n";
        // std::cout << "result = "; reinterpret_cast<value_t const &>(tmp).print(std::cout) << "\n";
        return tmp;
    }
};

// ----------------------------------------------------------------------------
// epi16
// ----------------------------------------------------------------------------

// template <simd_offset_t>
// struct in_lane_selector
// {
//     static constexpr simd_offset_t mask{0x00, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00, 0x01,
//                                      0x00, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00, 0x01,
//                                      0x00, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00, 0x01,
//                                      0x00, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00, 0x01}

//     simd_offset_t offsets{};

//     in_lane_selector() = delete;
//     in_lane_selector(simd_offset_t const & offsets) noexcept : offsets{offsets + mask}
//     {}

//     template <typename address_t>
//     constexpr auto operator()(address_t const & data) const noexcept {
//         return _mm512_shuffle_epi8(data, to_native(offsets));
//     }
// };
} // namespace detail
} // v1
} // namespace seqan::pairwise_aligner
