// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::affine_dp_algorithm.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <functional>
#include <seqan3/std/type_traits>

#include <pairwise_aligner/pairwise_aligner.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <template <typename ...> typename pairwise_aligner_t, typename ...policies_t>
class affine_dp_algorithm : public pairwise_aligner_t<affine_dp_algorithm<pairwise_aligner_t, policies_t...>>,
                            protected policies_t...
{
private:

    using base_t = pairwise_aligner_t<affine_dp_algorithm<pairwise_aligner_t, policies_t...>>;

    friend base_t;

    alignas(32) scalar_score_t gap_extend_score = -1;
    alignas(32) scalar_score_t gap_open_score = -10;

    alignas(32) scalar_score_t match_score = 4;
    alignas(32) scalar_score_t mismatch_score = -5;

    simd_score_t gap_extend_score_simd{gap_extend_score};
    simd_score_t gap_open_score_simd{gap_open_score};

    simd_score_t match_score_simd{match_score};
    simd_score_t mismatch_score_simd{mismatch_score};

public:
    affine_dp_algorithm() = default;

protected:

    template <std::ranges::viewable_range sequence_t, typename dp_vector_t>
        requires std::ranges::forward_range<sequence_t>
    auto initialise_row_vector(sequence_t && sequence, dp_vector_t & dp_vector)
    {
        using cell_t = typename std::remove_cvref_t<dp_vector_t>::value_type;
        using score_t = std::tuple_element_t<0, cell_t>;

        size_t column_index = 0;
        auto init_columns_strategy = [&] (cell_t & affine_cell)
        {
            if (column_index == 0)
                affine_cell = cell_t{score_t{0}, score_t{0}};
            else
                affine_cell = cell_t{score_t{static_cast<scalar_score_t>(gap_open_score + column_index * gap_extend_score)},
                                     score_t{static_cast<scalar_score_t>(gap_open_score + column_index * gap_extend_score +
                                                                         gap_open_score+ gap_extend_score)}};
            ++column_index;
        };

        return dp_vector.initialise(std::forward<sequence_t>(sequence), init_columns_strategy);
    }

    template <std::ranges::viewable_range sequence_t, typename dp_vector_t>
        requires std::ranges::forward_range<sequence_t>
    auto initialise_column_vector(sequence_t && sequence, dp_vector_t & dp_vector) const
    {
        using cell_t = typename std::remove_cvref_t<dp_vector_t>::value_type;
        using score_t = std::tuple_element_t<0, cell_t>;

        size_t row_index = 0;
        auto init_rows_strategy = [&] (cell_t & affine_cell)
        {
            if (row_index == 0)
                affine_cell = cell_t{score_t{0}, score_t{0}};
            else
                affine_cell = cell_t{score_t{static_cast<scalar_score_t>(gap_open_score + row_index * gap_extend_score)},
                                     score_t{static_cast<scalar_score_t>(gap_open_score + row_index * gap_extend_score +
                                                                         gap_open_score + gap_extend_score)}};
            ++row_index;
        };

        return dp_vector.initialise(std::forward<sequence_t>(sequence), init_rows_strategy);
    }

    template <typename row_cell_t, typename column_cell_t>
    auto initialise_column(row_cell_t & current_row_cell, column_cell_t & first_column_cell) const noexcept
    {
        using std::max;
        using score_t = std::tuple_element_t<0, row_cell_t>;

        std::pair cache{get<0>(first_column_cell), get<1>(current_row_cell)};
        get<0>(first_column_cell) = get<0>(current_row_cell);

        // horizontal_cache = horizontal_cache + as_derived().gap_extend_score;
        if constexpr (!std::integral<score_t>)
        {
            get<1>(first_column_cell) = max(static_cast<score_t>(cache.first + gap_open_score_simd),
                                            static_cast<score_t>(get<1>(first_column_cell) + gap_extend_score_simd));
        }
        else
        {
            get<1>(first_column_cell) = max(static_cast<score_t>(cache.first + gap_open_score),
                                            static_cast<score_t>(get<1>(first_column_cell) + gap_extend_score));
        }
        // // update the previous cache value.
        // get<0>(first_column_cell) = get<1>(current_row_cell);
        // // initialise the column value.
        // get<1>(first_column_cell) = max(get<1>(first_column_cell) + gap_extend_score,
        //                                 get<0>(first_column_cell) + gap_open_score);
        return cache;
    }

    template <typename row_cell_t, typename column_cell_t, typename cache_t>
    void finalise_column(row_cell_t & current_row_cell,
                         column_cell_t const & last_column_cell,
                         cache_t & cache) const noexcept
    {
        get<0>(current_row_cell) = get<0>(last_column_cell);
        get<1>(current_row_cell) = std::move(cache.second);
    }

    template <typename cache_t, typename seq1_val_t, typename seq2_val_t, typename dp_cell_t>
    auto compute_cell(cache_t & cache,
                      dp_cell_t & column_cell,
                      seq1_val_t const & seq1_val,
                      seq2_val_t const & seq2_val) const noexcept
    {
        using std::max;
        using score_t = std::tuple_element_t<0, dp_cell_t>;

        auto [next_diagonal, horizontal_score] = column_cell;
        // TODO: Should be a score model -> the code depends on the score model.
        if constexpr (std::integral<score_t>) {
            cache.first += (seq1_val == seq2_val) ? match_score : mismatch_score;
        }
        else {
            cache.first += compare_and_blend(seq1_val, seq2_val, match_score_simd, mismatch_score_simd);
        }
        cache.first = max(max(cache.first, cache.second), horizontal_score);
        get<0>(column_cell) = cache.first;
        cache.first += (gap_open_score + gap_extend_score);
        cache.second = max(static_cast<score_t>(cache.second + gap_extend_score), cache.first);
        get<1>(column_cell) = max(static_cast<score_t>(horizontal_score + gap_extend_score), cache.first);
        cache.first = next_diagonal; // cache score
    }
};

// ----------------------------------------------------------------------------
// Define the concrete parwise aligner instances.
// ----------------------------------------------------------------------------

template <typename ...policies_t>
using pairwise_aligner_affine = affine_dp_algorithm<pairwise_aligner, policies_t...>;

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
