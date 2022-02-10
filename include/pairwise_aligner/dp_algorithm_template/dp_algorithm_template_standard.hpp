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

        auto [dp_matrix_column] = base_t::initialise_policies();
        auto tracker = base_t::initialise_tracker();
        auto scorer = base_t::initialise_substitution_scheme();

        using block_sequence1_t = decltype(seqan3::views::slice(transformed_seq1, 0, 1));
        using block_sequence1_collection_t = std::vector<block_sequence1_t>;

        block_sequence1_collection_t block_sequences1{};
        block_sequences1.reserve(dp_column.size());

        { // Slice the first sequence according to the column size of each dp block.
            size_t offset = 0;
            std::ranges::for_each(std::views::iota(0ull, dp_column.size()), [&] (size_t const index) {
                size_t const column_size = dp_column[index].size() - 1;
                block_sequences1.emplace_back(seqan3::views::slice(transformed_seq1, offset, offset + column_size));
                offset += column_size;
            });
        }

        // ----------------------------------------------------------------------------
        // Recursion
        // ----------------------------------------------------------------------------

        size_t row_offset{};
        for (size_t column_idx = 0; column_idx < dp_row.size(); ++column_idx) {
            auto current_column = dp_matrix_column(dp_column, dp_row[column_idx], std::move(scorer), tracker);
            size_t const row_size = current_column.row().size() - 1;
            auto block_sequence2 = seqan3::views::slice(transformed_seq2, row_offset, row_offset + row_size);
            for (size_t block_idx = 0; block_idx < current_column.size(); ++block_idx) {
                auto dp_block = current_column[block_idx];
                base_t::compute_block(block_sequences1[block_idx], block_sequence2, dp_block);
            }
            row_offset += row_size;
        }

        // ----------------------------------------------------------------------------
        // Create result
        // ----------------------------------------------------------------------------

        return base_t::make_result(std::move(tracker),
                                   std::forward<sequence1_t>(sequence1),
                                   std::forward<sequence2_t>(sequence2),
                                   std::move(dp_column),
                                   std::move(dp_row));
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
