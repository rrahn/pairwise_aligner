// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::aligner_result_batch.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <memory>
#include <seqan3/std/ranges>
#include <tuple>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename aligner_result_t>
class aligner_result_bulk
{
    using score_type = typename aligner_result_t::score_type;

    std::shared_ptr<aligner_result_t> _result{};
    size_t _index{};

public:

    aligner_result_bulk() = default;
    explicit aligner_result_bulk(std::shared_ptr<aligner_result_t> result, size_t const index) noexcept :
        _result{std::move(result)},
        _index{index}
    {}

    auto const & dp_column() const & noexcept
    {
        return _result->dp_column();
    }

    auto const & dp_row() const & noexcept
    {
        return _result->dp_row();
    }

    auto const & sequence1() const noexcept
    {
        return *(_result->sequence1()[_index]);
    }

    auto const & sequence2() const noexcept
    {
        return *(_result->sequence2()[_index]);
    }

    auto score() const noexcept
    {
        return original_score();
    }

private:

    constexpr std::tuple<size_t, size_t, size_t> projected_coordinate() const noexcept
    {
        size_t const original_row_dim = std::ranges::distance(sequence1());
        size_t const original_column_dim = std::ranges::distance(sequence2());

        size_t offset = std::min<size_t>(dp_column().size() - 1 - original_row_dim,
                                         dp_row().size() -1 - original_column_dim);
        return std::tuple{original_row_dim + offset, original_column_dim + offset, offset};
    }

    constexpr auto original_score() const noexcept
    {
        using scalar_t = typename score_type::value_type;

        auto && [row_idx, col_idx, offset] = projected_coordinate();

        scalar_t best_score = 0;
        if (row_idx == dp_column().size() - 1)
        { // only the gap model knows how to access the correct value
            best_score = get<0>(dp_row()[col_idx])[_index];
        }
        else
        {
            assert(col_idx == dp_row().size() - 1);
            best_score = get<0>(dp_column()[row_idx])[_index];
        }
        return  best_score - static_cast<scalar_t>(_result->_padding_score[_index] * offset);
    }
};

// namespace cpo
// {
// struct fn
// {
//     template <typename aligner_result_t>
//     auto operator()(std::shared_ptr<aligner_result_t> aligner,
//                     index const index) const noexcept
//     {
//         return value<aligner_result>{std::move(aligner), index};
//     }
// };

// } // namespace cpo
// } // namespace _result.

// inline constexpr _aligner_result_bulk::cpo::fn make_result_bulk{};
} // inline namespace v1
}  // namespace seqan::pairwise_aligner
