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

    using index_t = simd_score<int8_t, score_t::size>;

    // Do we only work on these types?
    static constexpr size_t diagonal_matrix_size_v = ((dimension * (dimension + 1)) / 2);
    static constexpr size_t data_size_v = (diagonal_matrix_size_v - 1 + index_t::size) / index_t::size;

    std::array<index_t, data_size_v> _data{};

public:

    static constexpr size_t dimension_v = dimension;
    using score_type = score_t;
    using index_type = index_t;
    using offset_type = int32_t;

    type() = default;

    template <typename substitution_matrix_t>
    constexpr explicit type(substitution_matrix_t const & matrix) // assume non-linearised matrix
    {
        size_t position{};
        for (size_t rank_h = 0; rank_h < dimension_v; ++rank_h) {
            for (size_t rank_v = rank_h; rank_v < dimension_v; ++rank_v, ++position) {
                _data[position / index_t::size][position % index_t::size] = matrix[rank_h][rank_v];
            }
        }
        assert(position == diagonal_matrix_size_v);
    }

    constexpr index_t const & operator[](offset_type const & offset) const noexcept
    {
        return _data[offset];
    }

    template <typename strip_t, size_t width_v>
        requires (std::same_as<std::ranges::range_value_t<strip_t>, index_type>)
    constexpr auto initialise_profile(strip_t && sequence_strip, strip_width_t<width_v> const &) const noexcept
    {
        return scoring_handler_striped_NxN<type, width_v>{*this, std::forward<strip_t>(sequence_strip)};
    }

    template <typename value1_t>
    score_type score(score_type const & last_diagonal,
                     [[maybe_unused]] value1_t const & value1,
                     index_type const & value2) const noexcept
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
