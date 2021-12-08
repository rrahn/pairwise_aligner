// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::affine_cell_row and seqan::pairwise_aligned::affine_cell_column.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <utility>

#include <pairwise_aligner/type_traits.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename score_t, dp_vector_order order>
struct affine_cell : public dp_cell_base<order>,
                     public std::pair<score_t, score_t>
{
    using base_t = std::pair<score_t, score_t>;
    using score_type = score_t;

    using base_t::base_t;

    constexpr score_t const & score() const noexcept
    {
        return this->first;
    }
};

// template <typename score_t = int32_t>
// struct affine_cell_row : public affine_cell_base<score_t, dp_vector_order::row>
// {};

// template <typename score_t = int32_t>
// struct affine_cell_column : public affine_cell_base<score_t, dp_vector_order::column>
// {};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
