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

    static constexpr simd_offset_t lane_selector{static_cast<scalar_offset_t>(0x80)};

    __m256i first_mask;
    __m256i second_mask;

    simd_selector() = delete;
    simd_selector(simd_offset_t const & offsets) noexcept
    {
        // Check that no value is greater than or equal to 128 (0x80).
        // assert([&](){
        //     for (size_t i = 0; i < simd_offset_t::size; ++i) {
        //         if (offsets[i] >= 0x80)
        //             return false;
        //     }
        //     return true;
        // }());
        // what happens:
        // index:[0, 1, 2, 3 ... 15][16, 17, 18, ... 31] => one register
        // address -> if (addr[i] < 16) take from low
        //            else take from high
        // high bit is now based on 4-th bit position

        // set bit at position 4 to position 7 and mask out all other bits.
        __m256i tmp = _mm256_and_si256(_mm256_slli_si256(reinterpret_cast<__m256i const &>(offsets), 3),
                                       reinterpret_cast<__m256i const &>(lane_selector));
        // add the offset (assume bits not).
        first_mask = _mm256_or_si256(tmp, reinterpret_cast<__m256i const &>(offsets));
        second_mask = _mm256_xor_si256(first_mask, reinterpret_cast<__m256i const &>(lane_selector));

        std::cout << "offsets = "; offsets.print(std::cout) << "\n";
        std::cout << "first_mask = "; reinterpret_cast<simd_offset_t const &>(first_mask).print(std::cout) << "\n";
        std::cout << "second_mask = "; reinterpret_cast<simd_offset_t const &>(second_mask).print(std::cout) << "\n";
    }

    template <typename value_t>
    constexpr __m256i operator()(address_t<value_t> const & address) const noexcept {
        std::cout << "address = "; address[0].print(std::cout) << "\n";
        __m256i permuted_address = _mm256_permute4x64_epi64(reinterpret_cast<__m256i const &>(address[0]), 0b0100'1110);
        std::cout << "permuted address = "; reinterpret_cast<value_t const &>(permuted_address).print(std::cout) << "\n";

        __m256i tmp1 = _mm256_shuffle_epi8(reinterpret_cast<__m256i const &>(address[0]), first_mask);
        __m256i tmp2 = _mm256_shuffle_epi8(permuted_address, second_mask);
        __m256i tmp = _mm256_add_epi8(_mm256_shuffle_epi8(reinterpret_cast<__m256i const &>(address[0]), first_mask),
                                       _mm256_shuffle_epi8(permuted_address, second_mask));
        std::cout << "tmp1 = "; reinterpret_cast<value_t const &>(tmp1).print(std::cout) << "\n";
        std::cout << "tmp2 = "; reinterpret_cast<value_t const &>(tmp2).print(std::cout) << "\n";
        std::cout << "tmp = "; reinterpret_cast<value_t const &>(tmp).print(std::cout) << "\n";
        return tmp;
    }
};

// template <typename simd_offset_t, size_t operand_count>
//     requires (operand_count > 32 && operand_count <= 64)
// struct simd_selector<simd_offset_t, selector_tag<operand_count, 8, 256>>
// {
//     static constexpr bool in_lane_shuffle = false;
//     static constexpr size_t max_operand_count = 64;

//     template <typename value_t>
//     using address_t = std::array<value_t, 2>;
//     using scalar_offset_t = typename simd_offset_t::value_type;

//     static constexpr simd_offset_t mask1{
//         0x00,0xF0,0x02,0xF0,0x04,0xF0,0x06,0xF0,0x08,0xF0,0x0A,0xF0,0x0C,0xF0,0x0E,0xF0,
//         0x10,0xF0,0x12,0xF0,0x14,0xF0,0x16,0xF0,0x18,0xF0,0x1A,0xF0,0x1C,0xF0,0x1E,0xF0,
//         0x20,0xF0,0x22,0xF0,0x24,0xF0,0x26,0xF0,0x28,0xF0,0x2A,0xF0,0x2C,0xF0,0x2E,0xF0,
//         0x30,0xF0,0x32,0xF0,0x34,0xF0,0x36,0xF0,0x38,0xF0,0x3A,0xF0,0x3C,0xF0,0x3E,0xF0
//     };

//     static constexpr simd_offset_t mask2{
//         0xF0,0x00,0xF0,0x02,0xF0,0x04,0xF0,0x06,0xF0,0x08,0xF0,0x0A,0xF0,0x0C,0xF0,0x0E,
//         0xF0,0x10,0xF0,0x12,0xF0,0x14,0xF0,0x16,0xF0,0x18,0xF0,0x1A,0xF0,0x1C,0xF0,0x1E,
//         0xF0,0x20,0xF0,0x22,0xF0,0x24,0xF0,0x26,0xF0,0x28,0xF0,0x2A,0xF0,0x2C,0xF0,0x2E,
//         0xF0,0x30,0xF0,0x32,0xF0,0x34,0xF0,0x36,0xF0,0x38,0xF0,0x3A,0xF0,0x3C,0xF0,0x3E
//     };

//     __m512i offsets_even{};
//     __m512i offsets_uneven{};
//     __m512i shuffle_mask1{};
//     __m512i shuffle_mask2{};

//     simd_selector() = delete;
//     simd_selector(simd_offset_t const & offsets) noexcept
//     {
//         offsets_even = _mm512_srli_epi16(reinterpret_cast<__m512i const &>(offsets), 1);
//         offsets_uneven = _mm512_srli_epi16(reinterpret_cast<__m512i const &>(offsets), 9);

//         shuffle_mask1 = reinterpret_cast<__m512i &&>((offsets & simd_offset_t{static_cast<scalar_offset_t>(1)}) + mask1);
//         shuffle_mask2 = reinterpret_cast<__m512i &&>((offsets & simd_offset_t{static_cast<scalar_offset_t>(1)}) + mask2);
//     }

//     template <typename value_t>
//     constexpr auto operator()(address_t<value_t> const & address) const noexcept
//     {
//         // Load the scores in 2x8bit packs from the a and b register.
//         __m512i tmp_even = _mm512_permutex2var_epi16(reinterpret_cast<__m512i const &>(address[0]),
//                                                      offsets_even,
//                                                      reinterpret_cast<__m512i const &>(address[1]));

//         __m512i tmp_uneven = _mm512_permutex2var_epi16(reinterpret_cast<__m512i const &>(address[0]),
//                                                        offsets_uneven,
//                                                        reinterpret_cast<__m512i const &>(address[1]));


//         // This step is necessary to put the corresponding word
//         return _mm512_mask_add_epi8(_mm512_shuffle_epi8(tmp_even, shuffle_mask1),
//                                     0xAAAAAAAAAAAAAAAA,
//                                     _mm512_shuffle_epi8(tmp_uneven, shuffle_mask2),
//                                     __m512i{});
//     }
// };

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

// ----------------------------------------------------------------------------
// epi32
// ----------------------------------------------------------------------------

} // namespace detail
} // v1
} // namespace seqan::pairwise_aligner
