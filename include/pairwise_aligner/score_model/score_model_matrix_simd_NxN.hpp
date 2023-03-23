// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::score_model_matrix_simd_NxN.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <seqan3/std/concepts>

#include <pairwise_aligner/simd/simd_base.hpp>
#include <pairwise_aligner/simd/simd_index_map.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

struct offset_transform
{
    size_t alphabet_size{};
    size_t matrix_size{};

    template <typename simd_rank_t>
    constexpr auto operator()(simd_rank_t const & simd_ranks) const noexcept
    {
        // Original formula to compute the start offset of this symbol in the diagonal score matrix:
        //  (sigma * (sigma + 1))/2 - (((sigma - i) * (sigma - i + 1))/2)
        // The term in the dividend is always even such that the division by two will always yield an even number,
        // which is not truncated by the integer arithmetic.
        // However, working with bytes prohibits us to compute the multiplication first, as it yields an arithmetic
        // overflow.
        // Instead we first detect which of the ranks is uneven and use this bit as a scale.
        // Then by adding this scale we always guarantee, that we can divide the even term by two.
        // Subsequently, we add the inverse scale by xor-ing it with one and multiply the second term to the result.
        // Thus we avoid an arithmetic overflow and still get the correct result.
        simd_rank_t diff = simd_rank_t{static_cast<typename simd_rank_t::value_type>(alphabet_size)} - simd_ranks;
        simd_rank_t scale = diff & 1;
        return std::pair{
            simd_ranks,
            simd_rank_t{static_cast<typename simd_rank_t::value_type>(matrix_size)} -
                (divide_by_two(diff + scale) * (diff + (scale ^ 1)))
        };
    }

private:
    template <typename simd_rank_t>
    constexpr auto divide_by_two(simd_rank_t && rank) const noexcept
    {
        // There is no `shift right` intrinsic for 8-bit packed simd vectors and the compiler seems to emulate this by
        // unpacking the values to 16-bit packed simd vectors shifting them and packing them pack to 8-bit.
        // It is much faster to just reinterpret the vector as 16-bit packed vector and use the corresponding intrsinic
        // for this.
        // return rank >> 1;
        if constexpr (simd_rank_t::size_v == detail::max_simd_size) { // 8-bit packed simd vector.
            using scalar_rank_t = typename simd_rank_t::value_type;
            using scalar_rank_16_t = std::conditional_t<std::is_signed_v<scalar_rank_t>, int16_t, uint16_t>;
            using upcasted_simd_rank_t = simd_score<scalar_rank_16_t>;
            upcasted_simd_rank_t tmp = reinterpret_cast<upcasted_simd_rank_t const &>(rank) >> 1;
            return *reinterpret_cast<simd_rank_t *>(&tmp[0]);
        } else {
            return rank >> 1;
        }
    }
};
template <typename score_t, typename index_t, size_t dimension>
struct _score_model_matrix_simd_NxN
{
    class type;
};

template <typename score_t, typename index_t, size_t dimension>
using  score_model_matrix_simd_NxN = typename _score_model_matrix_simd_NxN<score_t, index_t, dimension>::type;

template <typename score_t, typename index_t, size_t dimension>
class _score_model_matrix_simd_NxN<score_t, index_t, dimension>::type
{
private:

    static constexpr size_t diagonal_matrix_size = ((dimension * (dimension + 1)) / 2);
    static constexpr size_t data_size_v = (diagonal_matrix_size - 1 + index_t::size_v) / index_t::size_v;

    using scalar_score_t = typename score_t::value_type;
    using scalar_index_t = typename index_t::value_type;
    using map_t = simd_index_map<scalar_score_t, scalar_index_t, diagonal_matrix_size>;

    map_t _data{};

public:
    using score_type = score_t;
    using index_type = index_t;

    type() = default;

    template <typename substitution_matrix_t>
    constexpr explicit type(substitution_matrix_t const & matrix) // assume non-linearised matrix
    {
        size_t position{};
        std::array<scalar_score_t, diagonal_matrix_size> tmp{};
        for (size_t rank_h = 0; rank_h < dimension; ++rank_h) {
            for (size_t rank_v = rank_h; rank_v < dimension; ++rank_v, ++position) {
                tmp[position] = matrix[rank_h][rank_v];
            }
        }
        _data = map_t{tmp};
        assert(position == diagonal_matrix_size);
    }

    template <typename value1_t, typename value2_t>
    score_type score(score_type const & last_diagonal, value1_t const & value1, value2_t const & value2) const noexcept
    {
        return last_diagonal + gather(value1, value2);
    }

    // TODO: Refactor into separate factory CPO.
    constexpr type make_substitution_scheme() const noexcept
    {
        return *this;
    }

private:
    template <typename value1_t, typename value2_t>
    constexpr score_type gather(value1_t const & value1, value2_t const & value2) const noexcept
    {
        auto && [column_rank, column_offset] = value1;
        auto && [row_rank, row_offset] = value2;
        auto matrix_offset = min(row_offset, column_offset) + absolute_difference(row_rank, column_rank);

        return _data[matrix_offset];
    }

    template <typename simd_rank_t>
    constexpr simd_rank_t absolute_difference(simd_rank_t const & row_ranks, simd_rank_t const & column_ranks) const noexcept
    {
        using signed_simd_rank_t = detail::make_signed_t<simd_rank_t>;
        signed_simd_rank_t tmp = abs(reinterpret_cast<signed_simd_rank_t const &>(row_ranks) -
                                     reinterpret_cast<signed_simd_rank_t const &>(column_ranks));
        return *reinterpret_cast<simd_rank_t *>(&tmp[0]);
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
