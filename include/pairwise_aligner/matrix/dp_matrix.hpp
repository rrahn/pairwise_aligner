// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_matrix_policies.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename dp_matrix_column_closure_t>
struct _dp_matrix_policies
{
    struct type;
};

template <typename dp_matrix_column_closure_t>
using dp_matrix_policies = typename _dp_matrix_policies<dp_matrix_column_closure_t>::type;

template <typename dp_matrix_column_closure_t>
struct _dp_matrix_policies<dp_matrix_column_closure_t>::type
{
    dp_matrix_column_closure_t dp_matrix_column_closure;

    constexpr auto make_policies() const noexcept {
        return std::forward_as_tuple(dp_matrix_column_closure);
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
