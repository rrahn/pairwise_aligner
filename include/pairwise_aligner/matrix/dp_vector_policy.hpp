// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_vector_policy.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename column_vector_t, typename row_vector_t>
class dp_vector_policy
{
private:

    column_vector_t _column_vector;
    row_vector_t _row_vector;

public:

    dp_vector_policy() = default;
    dp_vector_policy(column_vector_t column_vector, row_vector_t row_vector) noexcept :
        _column_vector{std::move(column_vector)},
        _row_vector{std::move(row_vector)}
    {}

    constexpr column_vector_t column_vector() const noexcept
    {
        return _column_vector;
    }

    constexpr row_vector_t row_vector() const noexcept
    {
        return _row_vector;
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
