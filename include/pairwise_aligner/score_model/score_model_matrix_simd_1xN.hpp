// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::score_model_matrix_simd_1xN.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <seqan3/std/concepts>

#include <pairwise_aligner/score_model/strip_width.hpp>
#include <pairwise_aligner/score_model/scoring_handler_striped_1xN.hpp>
#include <pairwise_aligner/simd/simd_base.hpp>
#include <pairwise_aligner/simd/simd_rank_selector.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename score_t, size_t dimension>
struct _score_model_matrix_simd_1xN
{
    class type;
};

template <typename score_t, size_t dimension>
using  score_model_matrix_simd_1xN = typename _score_model_matrix_simd_1xN<score_t, dimension>::type;

template <typename score_t, size_t dimension>
class _score_model_matrix_simd_1xN<score_t, dimension>::type :
    protected detail::simd_rank_selector_t<simd_score<int8_t, score_t::size>>
{
private:

    using index_t = simd_score<int8_t, score_t::size>;
    using simd_rank_selector_t = detail::simd_rank_selector_t<index_t>;
    using typename simd_rank_selector_t::rank_map_t;

    std::array<rank_map_t, dimension> _matrix{};

public:

    static constexpr size_t dimension_v = dimension;
    using score_type = score_t;
    using index_type = index_t;
    using offset_type = std::pair<int32_t, index_t>;

    type() = default;

    template <typename substitution_matrix_t> // TODO: does this remain scalar?
    constexpr explicit type(substitution_matrix_t const & matrix)
    {
        constexpr size_t chunk_size = (dimension_v - 1 + index_type::size) / index_type::size;
        for (size_t symbol_rank = 0; symbol_rank < dimension_v; ++symbol_rank) { // we move over the substitution_matrix
            std::array<index_type, chunk_size> tmp;
            for (size_t i = 0; i < dimension_v; ++i) {
                auto [index, offset] = std::pair{i / index_type::size, i % index_type::size};
                tmp[index][offset] = matrix[symbol_rank][i];
            }
            _matrix[symbol_rank] = simd_rank_selector_t::initialise_rank_map(tmp);
        }
    }

    constexpr index_type operator[](offset_type const & offset) const noexcept
    {
        return simd_rank_selector_t::select_rank_for(_matrix[offset.first], offset.second);
    }

    template <typename strip_t, size_t width_v>
        requires (std::same_as<std::ranges::range_value_t<strip_t>, index_type>)
    constexpr auto initialise_profile(strip_t && sequence_strip, strip_width_t<width_v> const &) const noexcept
    {
        return scoring_handler_striped_1xN<type, width_v>{*this, std::forward<strip_t>(sequence_strip)};
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
