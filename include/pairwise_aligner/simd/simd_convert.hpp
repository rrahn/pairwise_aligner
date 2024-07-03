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

#include <immintrin.h>

#include <pairwise_aligner/simd/simd_base.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace detail {

class simd_convert_base
{
private:

    template <typename simd_t>
    static constexpr bool has_8_bit_packing = seqan3::simd_traits<simd_t>::length ==
                                                    seqan3::simd_traits<simd_t>::max_length;

    template <typename simd_t>
    static constexpr bool has_16_bit_packing = seqan3::simd_traits<simd_t>::length ==
                                                    seqan3::simd_traits<simd_t>::max_length / 2;

    template <typename simd_t>
    static constexpr bool has_32_bit_packing = seqan3::simd_traits<simd_t>::length ==
                                                    seqan3::simd_traits<simd_t>::max_length / 4;

public:

    template <typename target_simd_t, typename source_simd_t, size_t source_count>
    constexpr void merge_into(target_simd_t & target,
                              std::array<source_simd_t, source_count> const & source_array) const noexcept
    {
        if constexpr (has_8_bit_packing<target_simd_t> &&
                      (has_32_bit_packing<source_simd_t> || has_32_bit_packing<source_simd_t>)) {
            merge_into_intrinsics(target, source_array);
        } else { // try auto-vectorisation.
            merge_into_auto(target,
                            source_array,
                            std::make_index_sequence<seqan3::simd_traits<target_simd_t>::length>());
        }
    }

    template <typename target_simd_t, size_t target_count, typename source_simd_t>
    constexpr void expand_into(std::array<target_simd_t, target_count> & target_array,
                               source_simd_t const & source) const noexcept
    {
        if constexpr (seqan3::simd_traits<target_simd_t>::max_length >= 32 &&
                      (has_16_bit_packing<target_simd_t> || has_32_bit_packing<target_simd_t>) &&
                      has_8_bit_packing<source_simd_t>) {
            expand_into_intrinsics(target_array, source);
        } else { // try auto-vectorisation.
            expand_into_auto(target_array,
                             source,
                             std::make_index_sequence<seqan3::simd_traits<source_simd_t>::length>());
        }
    }

private:

    // ----------------------------------------------------------------------------
    // AVX512 intrinsics
    // ----------------------------------------------------------------------------

    // From 16 bit to 8 bit
    template <typename target_simd_t, typename source_simd_t, size_t source_count>
        requires (seqan3::simd_traits<target_simd_t>::max_length == 64 && has_16_bit_packing<source_simd_t>)
    constexpr void merge_into_intrinsics(target_simd_t & target,
                                         std::array<source_simd_t, source_count> const & source_array) const noexcept
    {
        target =
            reinterpret_cast<target_simd_t>(
                _mm512_inserti64x4(
                    _mm512_castsi256_si512(_mm512_cvtepi16_epi8(reinterpret_cast<__m512i const &>(source_array[0]))),
                    _mm512_cvtepi16_epi8(reinterpret_cast<__m512i const &>(source_array[1])),
                    1));
    }

    // From 32 bit to 8 bit
    template <typename target_simd_t, typename source_simd_t, size_t source_count>
        requires (seqan3::simd_traits<target_simd_t>::max_length == 64 && has_32_bit_packing<source_simd_t>)
    constexpr void merge_into_intrinsics(target_simd_t & target,
                                         std::array<source_simd_t, source_count> const & source_array) const noexcept
    {
        constexpr auto lo_16_shuffle_mask{[]() -> std::array<int16_t, 32> {
            std::array<int16_t, 32> tmp{};
            for (size_t i = 0; i < tmp.size(); ++i)
                tmp[i] = i * 2;

            return tmp;
        }()};

        __m512i lo_32x16 = _mm512_permutex2var_epi16(reinterpret_cast<__m512i const &>(source_array[0]),
                                                     reinterpret_cast<__m512i const &>(lo_16_shuffle_mask),
                                                     reinterpret_cast<__m512i const &>(source_array[1]));
        __m512i hi_32x16 = _mm512_permutex2var_epi16(reinterpret_cast<__m512i const &>(source_array[2]),
                                                     reinterpret_cast<__m512i const &>(lo_16_shuffle_mask),
                                                     reinterpret_cast<__m512i const &>(source_array[3]));

        target = reinterpret_cast<target_simd_t>(
                _mm512_inserti64x4(_mm512_castsi256_si512(_mm512_cvtepi16_epi8(lo_32x16)),
                                   _mm512_cvtepi16_epi8(hi_32x16),
                                   1));
    }

    // From 8 bit to 16 bit
    template <typename target_simd_t, size_t target_count, typename source_simd_t>
        requires (seqan3::simd_traits<target_simd_t>::max_length == 64 && has_16_bit_packing<target_simd_t>)
    constexpr void expand_into_intrinsics(std::array<target_simd_t, target_count> & target_array,
                                          source_simd_t const & source) const noexcept
    {
        target_array[0] = reinterpret_cast<target_simd_t>(
                _mm512_cvtepi8_epi16(_mm512_castsi512_si256(reinterpret_cast<__m512i const &>(source))));
        target_array[1] =
            reinterpret_cast<target_simd_t>(
                _mm512_cvtepi8_epi16(_mm512_extracti64x4_epi64(reinterpret_cast<__m512i const &>(source), 1)));
    }

    // From 8 bit to 32 bit
    template <typename target_simd_t, size_t target_count, typename source_simd_t>
        requires (seqan3::simd_traits<target_simd_t>::max_length == 64 && has_32_bit_packing<target_simd_t>)
    constexpr void expand_into_intrinsics(std::array<target_simd_t, target_count> & target_array,
                                          source_simd_t const & source) const noexcept
    {
        __m512i lo_32x16_t = _mm512_cvtepi8_epi16(
                _mm512_castsi512_si256(reinterpret_cast<__m512i const &>(source)));
        __m512i hi_32x16_t = _mm512_cvtepi8_epi16(
                _mm512_extracti64x4_epi64(reinterpret_cast<__m512i const &>(source), 1));

        target_array[0] = reinterpret_cast<target_simd_t>(
                _mm512_cvtepi16_epi32(_mm512_castsi512_si256(lo_32x16_t)));
        target_array[1] = reinterpret_cast<target_simd_t>(
                _mm512_cvtepi16_epi32(_mm512_extracti64x4_epi64(lo_32x16_t, 1)));
        target_array[2] = reinterpret_cast<target_simd_t>(
                _mm512_cvtepi16_epi32(_mm512_castsi512_si256(hi_32x16_t)));
        target_array[3] = reinterpret_cast<target_simd_t>(
                _mm512_cvtepi16_epi32(_mm512_extracti64x4_epi64(hi_32x16_t, 1)));
    }

    // ----------------------------------------------------------------------------
    // AVX2 intrinsics
    // ----------------------------------------------------------------------------

    // From 16 bit to 8 bit
    template <typename target_simd_t, typename source_simd_t, size_t source_count>
        requires (seqan3::simd_traits<target_simd_t>::max_length == 32 && has_16_bit_packing<source_simd_t>)
    constexpr void merge_into_intrinsics(target_simd_t & target,
                                         std::array<source_simd_t, source_count> const & source_array) const noexcept
    {
        constexpr auto lo_8_mask{[]() -> std::array<int16_t, 16> {
            std::array<int16_t, 16> tmp{};
            tmp.fill(0xFF); // := 255
            return tmp;
        }()};

        // Zero out the upper 8 bits.
        __m256i lo_16x8 = _mm256_and_si256(reinterpret_cast<__m256i const &>(source_array[0]),
                                           reinterpret_cast<__m256i const &>(lo_8_mask));
        __m256i hi_16x8 = _mm256_and_si256(reinterpret_cast<__m256i const &>(source_array[1]),
                                           reinterpret_cast<__m256i const &>(lo_8_mask));
        // Pack and convert 16 bits (the lower 8 bits) into 8 bits of the target register using unsigned saturation.
        // The pack operation interleaves 4 elements of a with 4 elements of b, so the final result is permuted
        // back into the correct order: imm8 := 3, 1, 2, 0 := 11 01 10 00
        target = reinterpret_cast<target_simd_t>(_mm256_permute4x64_epi64(_mm256_packus_epi16(lo_16x8, hi_16x8),
                                                                          0b1101'1000));
    }

    // From 32 bit to 8 bit
    template <typename target_simd_t, typename source_simd_t, size_t source_count>
        requires (seqan3::simd_traits<target_simd_t>::max_length == 32 && has_32_bit_packing<source_simd_t>)
    constexpr void merge_into_intrinsics(target_simd_t & target,
                                         std::array<source_simd_t, source_count> const & source_array) const noexcept
    {
        constexpr auto lo_16_mask{[]() -> std::array<int32_t, 8> {
            std::array<int32_t, 8> tmp{};
            tmp.fill(0xFFFF); // := 65535
            return tmp;
        }()};

        constexpr auto lo_8_mask{[]() -> std::array<int16_t, 16> {
            std::array<int16_t, 16> tmp{};
            tmp.fill(0xFF); // := 255
            return tmp;
        }()};

        // Zero out the upper 16 bits.
        __m256i lo_32x4 = _mm256_and_si256(reinterpret_cast<__m256i const &>(source_array[0]),
                                           reinterpret_cast<__m256i const &>(lo_16_mask));
        __m256i hi_32x4 = _mm256_and_si256(reinterpret_cast<__m256i const &>(source_array[1]),
                                           reinterpret_cast<__m256i const &>(lo_16_mask));

        // Extract lower 16 bit from lo_32x4 and hi_32x4 forming the values of the lower half of the target vector.
        __m256i lo_16x8 = _mm256_permute4x64_epi64(_mm256_packus_epi32 (lo_32x4, hi_32x4), 0b1101'1000);

        // Zero out the upper 16 bits.
        lo_32x4 = _mm256_and_si256(reinterpret_cast<__m256i const &>(source_array[2]),
                                   reinterpret_cast<__m256i const &>(lo_16_mask));
        hi_32x4 = _mm256_and_si256(reinterpret_cast<__m256i const &>(source_array[3]),
                                   reinterpret_cast<__m256i const &>(lo_16_mask));

        // Extract lower 16 bit from lo_32x4 and hi_32x4 forming the values of the higher half of the target vector.
        __m256i hi_16x8 = _mm256_permute4x64_epi64(_mm256_packus_epi32 (lo_32x4, hi_32x4), 0b1101'1000);

        // Zero out the upper 8 bits.
        lo_16x8 = _mm256_and_si256(reinterpret_cast<__m256i const &>(lo_16x8),
                                   reinterpret_cast<__m256i const &>(lo_8_mask));
        hi_16x8 = _mm256_and_si256(reinterpret_cast<__m256i const &>(hi_16x8),
                                   reinterpret_cast<__m256i const &>(lo_8_mask));

        // Pack and convert 16 bits (the lower 8 bits) into 8 bits of the target register using unsigned saturation.
        // The pack operation interleaves 4 elements of a with 4 elements of b, so the final result is permuted
        // back into the correct order: imm8 := 3, 1, 2, 0 := 11 01 10 00
        target = reinterpret_cast<target_simd_t>(
                _mm256_permute4x64_epi64(_mm256_packus_epi16(lo_16x8, hi_16x8), 0b1101'1000));
    }

    // From 8 bit to 16 bit
    template <typename target_simd_t, size_t target_count, typename source_simd_t>
        requires (seqan3::simd_traits<target_simd_t>::max_length == 32 && has_16_bit_packing<target_simd_t>)
    constexpr void expand_into_intrinsics(std::array<target_simd_t, target_count> & target_array,
                                          source_simd_t const & source) const noexcept
    {
        target_array[0] = reinterpret_cast<target_simd_t>(
                _mm256_cvtepi8_epi16(_mm256_castsi256_si128(reinterpret_cast<__m256i const &>(source))));
        target_array[1] =
            reinterpret_cast<target_simd_t>(
                _mm256_cvtepi8_epi16(_mm256_extracti128_si256(reinterpret_cast<__m256i const &>(source), 1)));
    }

    // From 8 bit to 32 bit
    template <typename target_simd_t, size_t target_count, typename source_simd_t>
        requires (seqan3::simd_traits<target_simd_t>::max_length == 32 && has_32_bit_packing<target_simd_t>)
    constexpr void expand_into_intrinsics(std::array<target_simd_t, target_count> & target_array,
                                          source_simd_t const & source) const noexcept
    {
        __m256i lo_16x16v = _mm256_cvtepi8_epi16(
                _mm256_castsi256_si128(reinterpret_cast<__m256i const &>(source)));
        __m256i hi_16x16v = _mm256_cvtepi8_epi16(
                _mm256_extracti128_si256(reinterpret_cast<__m256i const &>(source), 1));

        target_array[0] = reinterpret_cast<target_simd_t>(
                _mm256_cvtepi16_epi32(_mm256_castsi256_si128(lo_16x16v)));
        target_array[1] = reinterpret_cast<target_simd_t>(
                _mm256_cvtepi16_epi32(_mm256_extracti128_si256(lo_16x16v, 1)));
        target_array[2] = reinterpret_cast<target_simd_t>(
                _mm256_cvtepi16_epi32(_mm256_castsi256_si128(hi_16x16v)));
        target_array[3] = reinterpret_cast<target_simd_t>(
                _mm256_cvtepi16_epi32(_mm256_extracti128_si256(hi_16x16v, 1)));
    }

    // ----------------------------------------------------------------------------
    // SSE4 intrinsics
    // ----------------------------------------------------------------------------

    // From 16 bit to 8 bit
    template <typename target_simd_t, typename source_simd_t, size_t source_count>
        requires (seqan3::simd_traits<target_simd_t>::max_length == 16 && has_16_bit_packing<source_simd_t>)
    constexpr void merge_into_intrinsics(target_simd_t & target,
                                         std::array<source_simd_t, source_count> const & source_array) const noexcept
    {
        constexpr auto lo_8_mask{[]() -> std::array<int16_t, 8> {
            std::array<int16_t, 8> tmp{};
            tmp.fill(0xFF); // := 255
            return tmp;
        }()};

        // Zero out the upper 8 bits.
        __m128i lo_16x4 = _mm_and_si128(reinterpret_cast<__m128i const &>(source_array[0]),
                                           reinterpret_cast<__m128i const &>(lo_8_mask));
        __m128i hi_16x4 = _mm_and_si128(reinterpret_cast<__m128i const &>(source_array[1]),
                                           reinterpret_cast<__m128i const &>(lo_8_mask));
        // Pack and convert 16 bits (the lower 8 bits) into 8 bits of the target register using unsigned saturation.
        target = reinterpret_cast<target_simd_t>(_mm_packus_epi16(lo_16x4, hi_16x4));
    }

    // From 32 bit to 8 bit
    template <typename target_simd_t, typename source_simd_t, size_t source_count>
        requires (seqan3::simd_traits<target_simd_t>::max_length == 16 && has_32_bit_packing<source_simd_t>)
    constexpr void merge_into_intrinsics(target_simd_t & target,
                                         std::array<source_simd_t, source_count> const & source_array) const noexcept
    {
        constexpr auto lo_16_mask{[]() -> std::array<int32_t, 4> {
            std::array<int32_t, 4> tmp{};
            tmp.fill(0xFFFF); // := 65535
            return tmp;
        }()};

        constexpr auto lo_8_mask{[]() -> std::array<int16_t, 8> {
            std::array<int16_t, 8> tmp{};
            tmp.fill(0xFF); // := 255
            return tmp;
        }()};

        // Zero out the upper 16 bits.
        __m128i lo_32x2 = _mm_and_si128(reinterpret_cast<__m128i const &>(source_array[0]),
                                           reinterpret_cast<__m128i const &>(lo_16_mask));
        __m128i hi_32x2 = _mm_and_si128(reinterpret_cast<__m128i const &>(source_array[1]),
                                           reinterpret_cast<__m128i const &>(lo_16_mask));

        // Extract lower 16 bit from lo_32x2 and hi_32x2 forming the values of the lower half of the target vector.
        __m128i lo_16x4 = _mm_packus_epi32(lo_32x2, hi_32x2);

        // Zero out the upper 16 bits.
        lo_32x2 = _mm_and_si128(reinterpret_cast<__m128i const &>(source_array[2]),
                                   reinterpret_cast<__m128i const &>(lo_16_mask));
        hi_32x2 = _mm_and_si128(reinterpret_cast<__m128i const &>(source_array[3]),
                                   reinterpret_cast<__m128i const &>(lo_16_mask));

        // Extract lower 16 bit from lo_32x4 and hi_32x4 forming the values of the higher half of the target vector.
        __m128i hi_16x4 = _mm_packus_epi32(lo_32x2, hi_32x2);

        // Zero out the upper 8 bits.
        lo_16x4 = _mm_and_si128(reinterpret_cast<__m128i const &>(lo_16x4),
                                   reinterpret_cast<__m128i const &>(lo_8_mask));
        hi_16x4 = _mm_and_si128(reinterpret_cast<__m128i const &>(hi_16x4),
                                   reinterpret_cast<__m128i const &>(lo_8_mask));

        // Pack and convert 16 bits (the lower 8 bits) into 8 bits of the target register using unsigned saturation.
        // The pack operation interleaves 4 elements of a with 4 elements of b, so the final result is permuted
        // back into the correct order: imm8 := 3, 1, 2, 0 := 11 01 10 00
        target = reinterpret_cast<target_simd_t>(_mm_packus_epi16(lo_16x4, hi_16x4));
    }

    // From 8 bit to 16 bit
    template <typename target_simd_t, size_t target_count, typename source_simd_t>
        requires (seqan3::simd_traits<target_simd_t>::max_length == 16 && has_16_bit_packing<target_simd_t>)
    constexpr void expand_into_intrinsics(std::array<target_simd_t, target_count> & target_array,
                                          source_simd_t const & source) const noexcept
    {
        target_array[0] = reinterpret_cast<target_simd_t>(
                            _mm_cvtepi8_epi16(reinterpret_cast<__m128i const &>(source)));
        target_array[1] = reinterpret_cast<target_simd_t>(
                _mm_cvtepi8_epi16(_mm_bsrli_si128(reinterpret_cast<__m128i const &>(source), 8)));
    }

    // From 8 bit to 32 bit
    template <typename target_simd_t, size_t target_count, typename source_simd_t>
        requires (seqan3::simd_traits<target_simd_t>::max_length == 16 && has_32_bit_packing<target_simd_t>)
    constexpr void expand_into_intrinsics(std::array<target_simd_t, target_count> & target_array,
                                          source_simd_t const & source) const noexcept
    {
        __m128i lo_16x8v = _mm_cvtepi8_epi16(reinterpret_cast<__m128i const &>(source));
        __m128i hi_16x8v = _mm_cvtepi8_epi16(_mm_bsrli_si128(reinterpret_cast<__m128i const &>(source), 8));

        target_array[0] = reinterpret_cast<target_simd_t>(_mm_cvtepi16_epi32(lo_16x8v));
        target_array[1] = reinterpret_cast<target_simd_t>(_mm_cvtepi16_epi32(_mm_bsrli_si128(lo_16x8v, 8)));
        target_array[2] = reinterpret_cast<target_simd_t>(_mm_cvtepi16_epi32(hi_16x8v));
        target_array[3] = reinterpret_cast<target_simd_t>(_mm_cvtepi16_epi32(_mm_bsrli_si128(hi_16x8v, 8)));
    }

    // ----------------------------------------------------------------------------
    // Auto-vectorisation
    // ----------------------------------------------------------------------------

    template <typename target_simd_t, typename source_simd_t, size_t source_count, size_t ...vec_idx>
    constexpr void merge_into_auto(target_simd_t & target,
                                    std::array<source_simd_t, source_count> const & source_array,
                                    std::index_sequence<vec_idx...> const &) const noexcept
    {
        using scalar_t = typename seqan3::simd_traits<target_simd_t>::scalar_type;
        constexpr size_t simd_size_v = seqan3::simd_traits<source_simd_t>::length;
        target = target_simd_t{
            (static_cast<scalar_t>(std::get<vec_idx / simd_size_v>(source_array)[vec_idx % simd_size_v]))...
        };
    }

    template <typename target_simd_t, size_t target_count, typename source_simd_t, size_t ...idx>
    constexpr void expand_into_auto(std::array<target_simd_t, target_count> & target_array,
                                    source_simd_t const & source,
                                    std::index_sequence<idx...> const &) const noexcept
    {
        using scalar_t = typename seqan3::simd_traits<target_simd_t>::scalar_type;
        constexpr size_t simd_size_v = seqan3::simd_traits<target_simd_t>::length;

        ((target_array[idx / simd_size_v][idx % simd_size_v] = static_cast<scalar_t>(source[idx])), ...);
    }
};

} // namespace detail
} // v1
} // namespace seqan::pairwise_aligner
