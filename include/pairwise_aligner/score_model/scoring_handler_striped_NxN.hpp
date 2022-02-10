// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::scoring_handler_striped_NxN.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <immintrin.h>

#include <pairwise_aligner/simd/simd_index_map.hpp>
#include <pairwise_aligner/simd/simd_score_type.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename substitution_matrix_t, size_t strip_width_v>
struct _scoring_handler_striped_NxN
{
    class type;
};

template <typename substitution_matrix_t, size_t strip_width_v>
using scoring_handler_striped_NxN = typename _scoring_handler_striped_NxN<substitution_matrix_t, strip_width_v>::type;

template <typename substitution_matrix_t, size_t strip_width_v>
class _scoring_handler_striped_NxN<substitution_matrix_t, strip_width_v>::type
{
private:

    static constexpr size_t matrix_size_v = substitution_matrix_t::diagonal_matrix_size_v;
    static constexpr size_t alphabet_size_v = substitution_matrix_t::dimension_v;

    static_assert(alphabet_size_v <= 22, "Only alphabets of size less than or equal to 22 are currently supported.");

    using simd_score_t = typename substitution_matrix_t::score_type;
    using index_t = typename substitution_matrix_t::index_type;
    using uindex_t = detail::make_unsigned_t<index_t>;
    using index_epu16_t = simd_score<uint16_t>;

    using scalar_index_t = typename index_t::value_type;
    using striped_scores_t = std::array<simd_score_t, strip_width_v>;
    using striped_offsets_t = std::array<uindex_t, strip_width_v>;
    using striped_ranks_t = std::array<index_t, strip_width_v>;

    striped_offsets_t _striped_offsets;
    striped_ranks_t _striped_ranks;
    substitution_matrix_t const & _matrix;

public:

    type() = delete;
    template <typename sequence_strip_t>
    explicit type(substitution_matrix_t const & matrix, sequence_strip_t && strip) noexcept : _matrix(matrix)
    {
        assert(std::ranges::size(strip) <= strip_width_v); // can't be larger.

        // Initialise the rank and offsets for the current strip.
        for (size_t index = 0; index < std::ranges::size(strip); ++index) {
                _striped_ranks[index] = strip[index].first;
                _striped_offsets[index] = strip[index].second;
        }
    }

    template <typename simd_symbol_t>
    constexpr auto scores_for(simd_symbol_t const & symbols) const noexcept
    {
        auto const & [ranks, offsets] = symbols;
        return for_each_lane(ranks, offsets, std::make_index_sequence<strip_width_v>());
    }

public:

    // constexpr uindex_t to_offset(index_t const & ranks) const noexcept {
    //     // Original formula to compute the start offset of this symbol in the diagonal score matrix:
    //     //  (sigma * (sigma + 1))/2 - (((sigma - i) * (sigma - i + 1))/2)
    //     // The term in the dividend is always even such that the division by two will always yield an even number,
    //     // which is not truncated by the integer arithmetic.
    //     // However, working with bytes prohibits us to compute the multiplication first, as it yields an arithmetic
    //     // overflow.
    //     // Instead we first detect which of the ranks is uneven and use this bit as a scale.
    //     // Then by adding this scale we always guarantee, that we can divide the even term by two.
    //     // Subsequently, we add the inverse scale by xor-ing it with one and multiply the second term to the result.
    //     // Thus we avoid an arithmetic overflow and still get the correct result.

    //     // TODO: Get rid of multiplication in 8-bit mode, because it is inefficient.
    //     // no 8-bit multiplication available, so the compiler unpacks it to 16 bit integer and then multiplies 2 times
    //     // using 16 bit.
        // index_t diff = index_t{static_cast<scalar_index_t>(alphabet_size_v)} - ranks;
        // index_t scale = diff & 1;
        // index_t offset = reinterpret_cast<index_t &&>((reinterpret_cast<index_epu16_t &&>(diff + scale) >> 1)) * ((diff + (scale ^ 1)));
        // return reinterpret_cast<uindex_t &&>(index_t{static_cast<scalar_index_t>(matrix_size_v)} - offset);
    // }

    template <size_t ...lane_idx>
    constexpr auto for_each_lane(index_t const & ranks,
                                 uindex_t const & offsets,
                                 std::index_sequence<lane_idx...> const &) const noexcept
    {
        auto gather = [&] (index_t const & column_rank,
                           uindex_t const & column_offset,
                           index_t const & row_rank,
                           uindex_t const & row_offset) {
            // either: max(row_rank, column_rank) - min(row_rank, column_rank) // three instructions vs. abs?
            uindex_t matrix_offset = min(row_offset, column_offset) + absolut_difference(row_rank, column_rank);

            return _matrix._data[matrix_offset];
            // #if defined(__AVX512VBMI__)
            //     __m512i a = select_scores(_matrix[0], matrix_offset, _matrix[1]);
            //     __m512i b = select_scores(_matrix[2], matrix_offset, _matrix[3]);
            // #else // require AVX512BW, to use less optimal version using cross lane permute with 16-bit packed operands.

            //     // can we precompute this?
            //     __m512i idx_32x16_even = (_mm512_srli_epi16(to_native(matrix_offset), 1));
            //     __m512i idx_32x16_uneven = (_mm512_srli_epi16(to_native(matrix_offset), 9));

            //     constexpr uindex_t even_mask{
            //         0x00,0xF0,0x02,0xF0,0x04,0xF0,0x06,0xF0,0x08,0xF0,0x0A,0xF0,0x0C,0xF0,0x0E,0xF0,
            //         0x10,0xF0,0x12,0xF0,0x14,0xF0,0x16,0xF0,0x18,0xF0,0x1A,0xF0,0x1C,0xF0,0x1E,0xF0,
            //         0x20,0xF0,0x22,0xF0,0x24,0xF0,0x26,0xF0,0x28,0xF0,0x2A,0xF0,0x2C,0xF0,0x2E,0xF0,
            //         0x30,0xF0,0x32,0xF0,0x34,0xF0,0x36,0xF0,0x38,0xF0,0x3A,0xF0,0x3C,0xF0,0x3E,0xF0
            //     };

            //     constexpr uindex_t uneven_mask{
            //         0xF0,0x00,0xF0,0x02,0xF0,0x04,0xF0,0x06,0xF0,0x08,0xF0,0x0A,0xF0,0x0C,0xF0,0x0E,
            //         0xF0,0x10,0xF0,0x12,0xF0,0x14,0xF0,0x16,0xF0,0x18,0xF0,0x1A,0xF0,0x1C,0xF0,0x1E,
            //         0xF0,0x20,0xF0,0x22,0xF0,0x24,0xF0,0x26,0xF0,0x28,0xF0,0x2A,0xF0,0x2C,0xF0,0x2E,
            //         0xF0,0x30,0xF0,0x32,0xF0,0x34,0xF0,0x36,0xF0,0x38,0xF0,0x3A,0xF0,0x3C,0xF0,0x3E
            //     };

            //     // Less instructions?
            //     uindex_t selector_even = (matrix_offset & uindex_t{(uint8_t)1}) + even_mask;
            //     uindex_t selector_uneven = (matrix_offset & uindex_t{(uint8_t)1}) + uneven_mask;

            //     __m512i a = select_scores(_matrix[0],
            //                               idx_32x16_even,
            //                               idx_32x16_uneven,
            //                               selector_even,
            //                               selector_uneven,
            //                               _matrix[1]);
            //     __m512i b = select_scores(_matrix[2],
            //                               idx_32x16_even,
            //                               idx_32x16_uneven,
            //                               selector_even,
            //                               selector_uneven,
            //                               _matrix[3]);

            // #endif // defined(__AVX512VBMI__)

            // return blend(matrix_offset.lt(uindex_t{static_cast<typename uindex_t::value_type>(0x80)}),
            //             to_packed(a),
            //             to_packed(b));
        };

        return striped_scores_t{gather(_striped_ranks[lane_idx], _striped_offsets[lane_idx], ranks, offsets)...};
    }

    template <typename simd_rank_t>
    constexpr auto absolut_difference(simd_rank_t const & row_ranks, simd_rank_t const & column_ranks) const noexcept
    {
        if constexpr (std::is_signed_v<typename simd_rank_t::value_type>)
            return reinterpret_cast<uindex_t &&>(abs(row_ranks - column_ranks));
        else
            return reinterpret_cast<uindex_t &&>(max(row_ranks, column_ranks) - min(row_ranks, column_ranks));
    }

    // #if defined(__AVX512VBMI__)
    // constexpr auto select_scores(index_t const & a, uindex_t const & idx, index_t const & b) const noexcept
    // {
    //     return _mm512_permutex2var_epi8(to_native(a), to_native(idx), to_native(b));
    // }
    // #else // Alternative because we don't have epi8 cross lane shuffle without AVX512_VBMI.
    // constexpr auto select_scores(index_t const & a,
    //                              __m512i const & idx_32x16_even,
    //                              __m512i const & idx_32x16_uneven,
    //                              uindex_t const & selector_even,
    //                              uindex_t const & selector_uneven,
    //                              index_t const & b) const noexcept
    // {
    //     // Load the scores in 2x8bit packs from the a and b register.
    //     // permute 16_bit is slow! maybe we can improve this?
    //     __m512i tmp_even = _mm512_permutex2var_epi16(to_native(a), (idx_32x16_even), to_native(b));
    //     __m512i tmp_uneven = _mm512_permutex2var_epi16(to_native(a), (idx_32x16_uneven), to_native(b));

    //     // This step is necessary to put the corresponding word
    //     return _mm512_mask_add_epi8(_mm512_shuffle_epi8(tmp_even, to_native(selector_even)),
    //                                 0xAAAAAAAAAAAAAAAA,
    //                                 _mm512_shuffle_epi8(tmp_uneven, to_native(selector_uneven)),
    //                                 __m512i{});
    // }
    // #endif // defined(__AVX512VBMI__)

    template <typename packed_t>
    static __m512i const & to_native(packed_t const & packed) noexcept
    {
        return reinterpret_cast<__m512i const &>(packed);
    }

    static index_t to_packed(__m512i const & native) noexcept
    {
        return reinterpret_cast<index_t const &>(native);
    }
};
} // inline namespace v1
} // namespace seqan::pairwise_aligner
