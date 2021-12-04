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

#include <pairwise_aligner/result/aligner_result.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename derived_t>
struct _dp_algorithm_template_standard
{
    class type;
};

template <typename derived_t>
using dp_algorithm_template_standard = typename _dp_algorithm_template_standard<derived_t>::type;

template <typename derived_t>
class _dp_algorithm_template_standard<derived_t>::type
{
protected:

    friend derived_t;

    template <std::ranges::forward_range sequence1_t,
              std::ranges::forward_range sequence2_t,
              typename dp_column_vector_t,
              typename dp_row_vector_t>
    auto run(sequence1_t && sequence1,
             sequence2_t && sequence2,
             dp_column_vector_t dp_column_vector,
             dp_row_vector_t dp_row_vector)
    {
        // ----------------------------------------------------------------------------
        // Initialisation
        // ----------------------------------------------------------------------------

        auto transformed_seq1 = as_derived().initialise_column_vector(sequence1, dp_column_vector);
        auto transformed_seq2 = as_derived().initialise_row_vector(sequence2, dp_row_vector);

        // ----------------------------------------------------------------------------
        // Recursion
        // ----------------------------------------------------------------------------

        for (size_t j = 0; j < std::ranges::size(transformed_seq2); ++j)
        {
            auto cache = as_derived().initialise_column(dp_row_vector[j+1], dp_column_vector[0]);

            size_t i = 0;
            for (; i < std::ranges::size(transformed_seq1); ++i)
                as_derived().compute_cell(cache, dp_column_vector[i+1], transformed_seq1[i], transformed_seq2[j]);

            as_derived().finalise_column(dp_row_vector[j+1], dp_column_vector[i], cache);
        }

        // ----------------------------------------------------------------------------
        // Create result
        // ----------------------------------------------------------------------------

        auto best_score = get<0>(dp_column_vector[std::ranges::size(transformed_seq1)]);
        return make_result(std::forward<sequence1_t>(sequence1),
                           std::forward<sequence2_t>(sequence2),
                           std::move(dp_column_vector),
                           std::move(dp_row_vector),
                           std::move(best_score));
    }

    derived_t const & as_derived() const noexcept
    {
        return static_cast<derived_t const &>(*this);
    }

    derived_t & as_derived() noexcept
    {
        return static_cast<derived_t &>(*this);
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
