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

#include <seqan3/std/concepts>
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

    affine_cell() = default;

    template <typename other_score_t>
        requires std::constructible_from<score_t, other_score_t const &>
    explicit affine_cell(affine_cell<other_score_t, order> const & other_cell) :
        base_t{other_cell.first, other_cell.second}
    {}

    template <typename other_score_t>
        requires std::constructible_from<score_t, other_score_t>
    explicit affine_cell(affine_cell<other_score_t, order> && other_cell) :
        base_t{std::move(other_cell.first), std::move(other_cell.second)}
    {}

    template <typename other_score_t>
        requires std::assignable_from<score_t &, other_score_t const &>
    affine_cell & operator=(affine_cell<other_score_t, order> const & other_cell)
    {
        static_cast<base_t &>(*this) = base_t{other_cell.first, other_cell.second};
        return *this;
    }

    template <typename other_score_t>
        requires std::assignable_from<score_t &, other_score_t>
    affine_cell & operator=(affine_cell<other_score_t, order> && other_cell)
    {
        static_cast<base_t &>(*this) = base_t{std::move(other_cell.first), std::move(other_cell.second)};
        return *this;
    }

    constexpr score_t & score() noexcept
    {
        return this->first;
    }

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

namespace std
{

template <typename score_t, seqan::pairwise_aligner::dp_vector_order order>
struct tuple_size<seqan::pairwise_aligner::affine_cell<score_t, order>> :
    tuple_size<typename seqan::pairwise_aligner::affine_cell<score_t, order>::base_t>
{};

template <size_t idx, typename score_t, seqan::pairwise_aligner::dp_vector_order order>
struct tuple_element<idx, seqan::pairwise_aligner::affine_cell<score_t, order>> :
    tuple_element<idx, typename seqan::pairwise_aligner::affine_cell<score_t, order>::base_t>
{};

} // namespace std
