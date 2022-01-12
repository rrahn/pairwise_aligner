// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::interface_one_to_one_single.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <seqan3/std/ranges>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename dp_algorithm_t>
struct _interface_one_to_one_single
{
    struct type;
};

template <typename dp_algorithm_t>
using interface_one_to_one_single = typename _interface_one_to_one_single<dp_algorithm_t>::type;

template <typename dp_algorithm_t>
struct _interface_one_to_one_single<dp_algorithm_t>::type : protected dp_algorithm_t
{
    using dp_algorithm_t::dp_algorithm_t;

    using dp_algorithm_t::column_vector;
    using dp_algorithm_t::row_vector;

    template <std::ranges::forward_range sequence1_t,
              std::ranges::forward_range sequence2_t>
        requires (std::ranges::viewable_range<sequence1_t> &&
                  std::ranges::viewable_range<sequence2_t>)
    auto compute(sequence1_t && sequence1, sequence2_t && sequence2)
    {
        return compute(std::forward<sequence1_t>(sequence1),
                       std::forward<sequence2_t>(sequence2),
                       column_vector(),
                       row_vector());
    }

    template <std::ranges::forward_range sequence1_t,
              std::ranges::forward_range sequence2_t,
              typename column_vector_t,
              typename row_vector_t>
        requires (std::ranges::viewable_range<sequence1_t> &&
                  std::ranges::viewable_range<sequence2_t>)
    auto compute(sequence1_t && sequence1,
                 sequence2_t && sequence2,
                 column_vector_t first_dp_column,
                 row_vector_t first_dp_row)
    {
        return dp_algorithm_t::run(std::forward<sequence1_t>(sequence1),
                                   std::forward<sequence2_t>(sequence2),
                                   std::move(first_dp_column),
                                   std::move(first_dp_row));
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
