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

#include <pairwise_aligner/score_model/strip_width.hpp>
#include <pairwise_aligner/score_model/scoring_handler_striped_NxN.hpp>
#include <pairwise_aligner/simd/simd_base.hpp>

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
        if constexpr (simd_rank_t::size == detail::max_simd_size) { // 8-bit packed simd vector.
            using scalar_rank_t = typename simd_rank_t::value_type;
            using scalar_rank_16_t = std::conditional_t<std::is_signed_v<scalar_rank_t>, int16_t, uint16_t>;
            using upcasted_simd_rank_t = simd_score<scalar_rank_16_t>;
            return reinterpret_cast<simd_rank_t &&>(reinterpret_cast<upcasted_simd_rank_t const &>(rank) >> 1);
        } else {
            return rank >> 1;
        }
    }
};
template <typename score_t, size_t dimension>
struct _score_model_matrix_simd_NxN
{
    class type;
};

template <typename score_t, size_t dimension>
using  score_model_matrix_simd_NxN = typename _score_model_matrix_simd_NxN<score_t, dimension>::type;

template <typename score_t, size_t dimension>
class _score_model_matrix_simd_NxN<score_t, dimension>::type
{
private:

    template <typename, size_t>
    friend class _scoring_handler_striped_NxN;

    using index_t = simd_score<uint8_t, score_t::size>;

    // Do we only work on these types?
    static constexpr size_t diagonal_matrix_size = ((dimension * (dimension + 1)) / 2);
    static constexpr size_t data_size_v = (diagonal_matrix_size - 1 + index_t::size) / index_t::size;

    using scalar_score_t = typename score_t::value_type;

    using map_t = simd_index_map<scalar_score_t, uint8_t, diagonal_matrix_size>;
    // std::array<index_t, data_size_v> _data{};
    map_t _data{};

public:

    static constexpr size_t dimension_v = dimension;
    using score_type = score_t;
    using index_type = index_t;
    // using offset_type = int32_t;

    type() = default;

    template <typename substitution_matrix_t>
    constexpr explicit type(substitution_matrix_t const & matrix) // assume non-linearised matrix
    {
        size_t position{};
        std::array<scalar_score_t, diagonal_matrix_size> tmp{};
        for (size_t rank_h = 0; rank_h < dimension_v; ++rank_h) {
            for (size_t rank_v = rank_h; rank_v < dimension_v; ++rank_v, ++position) {
                // _data[position / index_t::size][position % index_t::size] = matrix[rank_h][rank_v];
                tmp[position] = matrix[rank_h][rank_v];
            }
        }
        _data = map_t{tmp};
        assert(position == diagonal_matrix_size);
    }

    template <typename strip_t, size_t width>
        // requires (std::same_as<std::ranges::range_value_t<strip_t>, index_type>)
    constexpr auto initialise_profile(strip_t && sequence_strip, strip_width_t<width> const &) const noexcept
    {
        return scoring_handler_striped_NxN<type, width>{*this, std::forward<strip_t>(sequence_strip)};
    }

    template <typename value1_t>
    score_type score(score_type const & last_diagonal,
                     [[maybe_unused]] value1_t const & value1,
                     score_type const & value2) const noexcept
    {
        // Upcasting the index scores to the score type.
        return last_diagonal + static_cast<score_type>(value2);
    }

    // TODO: Refactor into separate factory CPO.
    constexpr type make_substitution_scheme() const noexcept
    {
        return *this;
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
