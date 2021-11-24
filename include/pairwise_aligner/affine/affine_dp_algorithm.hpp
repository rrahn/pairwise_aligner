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

    int32_t gap_open_score = -10;
    int32_t gap_extend_score = -1;

public:
    affine_dp_algorithm() = default;

protected:

    template <std::ranges::viewable_range sequence_t, typename dp_vector_t>
        requires std::ranges::forward_range<sequence_t>
    auto initialise_row_vector(sequence_t && sequence, dp_vector_t & dp_vector)
    {
        int column_index = 0;
        auto init_columns_strategy = [&] (affine_dp_cell_t & affine_cell)
        {
            if (column_index == 0)
                affine_cell = affine_dp_cell_t{0, 0, 0};
            else
                affine_cell = affine_dp_cell_t{gap_open_score + column_index * gap_extend_score,
                                               gap_open_score + column_index * gap_extend_score,
                                               gap_open_score + column_index * gap_extend_score + gap_open_score + gap_extend_score};
            ++column_index;
        };

        return dp_vector.initialise(std::forward<sequence_t>(sequence), init_columns_strategy);
    }

    template <std::ranges::viewable_range sequence_t, typename dp_vector_t>
        requires std::ranges::forward_range<sequence_t>
    auto initialise_column_vector(sequence_t && sequence, dp_vector_t & dp_vector) const
    {
        int row_index = 0;
        auto init_rows_strategy = [&] (affine_dp_cell_t & affine_cell)
        {
            if (row_index == 0)
                affine_cell = affine_dp_cell_t{0, 0, 0};
            else
                affine_cell = affine_dp_cell_t{gap_open_score + row_index * gap_extend_score,
                                               gap_open_score + row_index * gap_extend_score + gap_open_score + gap_extend_score,
                                               gap_open_score + row_index * gap_extend_score};
            ++row_index;
        };

        return dp_vector.initialise(std::forward<sequence_t>(sequence), init_rows_strategy);
    }

    template <typename dp_column_vector_t, typename dp_row_vector_t>
    auto initialise_cache(dp_column_vector_t & first_row_cell, dp_row_vector_t & current_column_cell) const noexcept
    {
        // int32_t & vertical_score = get<2>(dp_row_vector[j+1]); // get value
        // int32_t diagonal_score = get<0>(dp_column_vector[0]); // cache diagonal
        return std::tuple<int32_t, int32_t &>{get<0>(first_row_cell), get<2>(current_column_cell)};
    }

    template <typename cache_t, typename seq1_val_t, typename seq2_val_t, typename dp_cell_t>
    auto compute_cell(cache_t & cache,
                      dp_cell_t & column_cell,
                      seq1_val_t const & seq1_val,
                      seq2_val_t const & seq2_val) const noexcept
    {
        auto score = [] (auto val1, auto val2)
        {
            return val1 == val2 ? 4 : -5;
        };

        // std::cout << "[DEBUG] ";
        auto && [diagonal_score, vertical_score] = cache;
        // std::cout << "diagonal_score old = " << diagonal_score << " ";
        // call inline derived type function here
        int32_t & horizontal_score = get<1>(column_cell);
        diagonal_score += score(seq1_val, seq2_val);
        // std::cout << "diagonal_score new = " << diagonal_score << " ";
        // std::cout << "horizontal_score = " << horizontal_score << " ";
        // std::cout << "vertical_score = " << vertical_score << " ";
        int32_t best = std::max(std::max(diagonal_score, vertical_score), horizontal_score);
        diagonal_score = get<0>(column_cell); // cache score
        get<0>(column_cell) = best;
        // std::cout << "best = " << best << " ";

        // set the horizontal and vertical score
        best += gap_open_score;
        horizontal_score = std::max(horizontal_score, best) + gap_extend_score;
        vertical_score = std::max(vertical_score, best) + gap_extend_score;
        // std::cout << "\n";
        return cache;
    }
};

// ----------------------------------------------------------------------------
// Define the concrete parwise aligner instances.
// ----------------------------------------------------------------------------

template <typename ...policies_t>
using pairwise_aligner_affine = affine_dp_algorithm<pairwise_aligner, policies_t...>;

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
