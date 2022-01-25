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

        // ----------------------------------------------------------------------------
        // Recursion
        // ----------------------------------------------------------------------------

        auto tracker = base_t::initialise_tracker();
        base_t::initialise_block(dp_column, dp_row);
        base_t::compute_block(transformed_seq1, transformed_seq2, dp_column, dp_row, tracker);
        base_t::postprocess_block(dp_column, dp_row);

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
