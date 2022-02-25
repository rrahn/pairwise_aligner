// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_algorithm_template_standard.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <seqan3/std/ranges>

#include <seqan3/utility/views/slice.hpp>

#include <pairwise_aligner/dp_algorithm_template/dp_algorithm_template_base.hpp>
#include <pairwise_aligner/matrix/dp_matrix_cpo.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename algorithm_impl_t>
struct _dp_algorithm_template_standard
{
    class type;
};

template <typename algorithm_impl_t>
using dp_algorithm_template_standard = typename _dp_algorithm_template_standard<algorithm_impl_t>::type;

template <typename algorithm_impl_t>
class _dp_algorithm_template_standard<algorithm_impl_t>::type : public dp_algorithm_template_base<algorithm_impl_t>
{
protected:

    using base_t = dp_algorithm_template_base<algorithm_impl_t>;

    template <typename sequence1_t, typename sequence2_t, typename dp_column_t, typename dp_row_t>
    auto run(sequence1_t && sequence1, sequence2_t && sequence2, dp_column_t dp_column, dp_row_t dp_row) const
    {
        // ----------------------------------------------------------------------------
        // Initialisation
        // ----------------------------------------------------------------------------
        auto transformed_seq1 = base_t::initialise_column(sequence1, dp_column);
        auto transformed_seq2 = base_t::initialise_row(sequence2, dp_row);

        auto matrix = base_t::initialise_dp_matrix(dp_column, dp_row, transformed_seq1, transformed_seq2);
        // auto tracker = base_t::initialise_tracker();

        // using block_sequence1_t = decltype(seqan3::views::slice(transformed_seq1, 0, 1));
        // using block_sequence1_collection_t = std::vector<block_sequence1_t>;

        // block_sequence1_collection_t block_sequences1{};
        // block_sequences1.reserve(dp_column.size());

        // { // Slice the first sequence according to the column size of each dp block.
        //     size_t offset = 0;
        //     std::ranges::for_each(std::views::iota(0ull, dp_column.size()), [&] (size_t const index) {
        //         size_t const column_size = dp_column[index].size() - 1;
        //         block_sequences1.emplace_back(seqan3::views::slice(transformed_seq1, offset, offset + column_size));
        //         offset += column_size;
        //     });
        // }

        // ----------------------------------------------------------------------------
        // Recursion
        // ----------------------------------------------------------------------------

        for (std::ptrdiff_t column_idx = 0; column_idx < dp_matrix::column_count(matrix); ++column_idx) {
            // size_t const row_size = dp_row[column_idx].size() - 1;
            // auto block_sequence2 = seqan3::views::slice(transformed_seq2, row_offset, row_offset + row_size);
            auto current_column = dp_matrix::column_at(matrix, column_idx);
            // dp_column,
            //                                        dp_row[column_idx],
            //                                        base_t::initialise_substitution_scheme(),
            //                                        tracker,
            //                                        std::move(block_sequence2),
            //                                        base_t::lane_width());
            for (std::ptrdiff_t row_idx = 0; row_idx < dp_matrix::row_count(current_column); ++row_idx) {
                auto dp_block = dp_matrix::row_at(current_column, row_idx);
                base_t::compute_block(dp_block);
            }
            // row_offset += row_size;
        }

        // ----------------------------------------------------------------------------
        // Create result
        // ----------------------------------------------------------------------------

        return base_t::make_result(std::move(dp_matrix::tracker(matrix)),
                                   std::forward<sequence1_t>(sequence1),
                                   std::forward<sequence2_t>(sequence2),
                                   std::move(dp_column),
                                   std::move(dp_row));
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
