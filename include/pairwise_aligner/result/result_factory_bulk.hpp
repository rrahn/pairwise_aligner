// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::aligner_result.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <pairwise_aligner/result/aligner_result.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace _bulk_factory
{
template <typename base_value_t, typename score_t>
struct _value
{
    struct type;
};

template <typename base_value_t, typename score_t>
using value = typename _value<base_value_t, score_t>::type;

template <typename base_value_t, typename score_t>
struct _value<base_value_t, score_t>::type : public base_value_t
{
    score_t _padding_score;

    constexpr auto score_at(size_t const idx) const noexcept
    {
        using scalar_t = typename score_t::value_type;

        assert(this->sequence1()[idx].has_value());
        assert(this->sequence2()[idx].has_value());


        auto && [row_idx, col_idx, offset] = projected_coordinate(idx);
        assert(row_idx == this->dp_column().size() - 1 || col_idx == this->dp_row().size() - 1);

        scalar_t best_score = 0;
        if (row_idx == this->dp_column().size() - 1)
            best_score = get<0>(this->dp_row()[col_idx])[idx];
        else
            best_score = get<0>(this->dp_column()[row_idx])[idx];

        return best_score - static_cast<scalar_t>(_padding_score[idx] * offset);
    }

    constexpr std::tuple<size_t, size_t, size_t> projected_coordinate(size_t const idx) const noexcept
    {
        size_t const original_row_dim = std::ranges::distance(*(this->sequence1()[idx]));
        size_t const original_column_dim = std::ranges::distance(*(this->sequence2()[idx]));

        size_t offset = std::min<size_t>(this->dp_column().size() - 1 - original_row_dim,
                                         this->dp_row().size() -1 - original_column_dim);
        return std::tuple{original_row_dim + offset, original_column_dim + offset, offset};
    }
};

} // namespace _bulk_factory

template <typename score_t>
struct _result_factory_bulk
{
    struct type;
};

template <typename score_t>
using result_factory_bulk = typename _result_factory_bulk<score_t>::type;

template <typename score_t>
struct _result_factory_bulk<score_t>::type
{
    score_t _padding_score;

    template <typename ...args_t>
    auto operator()(args_t && ...args) const noexcept
    {
        auto base = std::invoke(result_factory_single{}, std::forward<args_t>(args)...);
        return _bulk_factory::value<decltype(base), score_t>{std::move(base), _padding_score};
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
