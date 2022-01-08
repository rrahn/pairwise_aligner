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

#include <seqan3/std/ranges>

#include <pairwise_aligner/dp_trailing_gaps.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace _aligner_result
{
template <typename sequence1_t, typename sequence2_t, typename dp_column_t, typename dp_row_t>
struct _value
{
    struct type;
};

template <typename sequence1_t, typename sequence2_t, typename dp_column_t, typename dp_row_t>
using value = typename _value<sequence1_t, sequence2_t, dp_column_t, dp_row_t>::type;

template <typename sequence1_t, typename sequence2_t, typename dp_column_t, typename dp_row_t>
struct _value<sequence1_t, sequence2_t, dp_column_t, dp_row_t>::type
{
    sequence1_t _sequence1;
    sequence2_t _sequence2;
    dp_column_t _dp_column;
    dp_row_t _dp_row;
    dp_trailing_gaps _column_trailing_gaps;
    dp_trailing_gaps _row_trailing_gaps;

    dp_column_t const & dp_column() const & noexcept
    {
        return _dp_column;
    }

    dp_column_t && dp_column() && noexcept
    {
        return std::move(_dp_column);
    }

    dp_row_t const & dp_row() const & noexcept
    {
        return _dp_row;
    }

    dp_row_t && dp_row() && noexcept
    {
        return std::move(_dp_row);
    }

    sequence1_t const & sequence1() const noexcept
    {
        return _sequence1;
    }

    sequence2_t const & sequence2() const noexcept
    {
        return _sequence2;
    }

    auto score() const noexcept
    {
        return score_impl();
    }

private:

    constexpr auto score_impl() const noexcept
    {
        auto best_score = dp_column()[dp_column().size() - 1].score();

        if (_row_trailing_gaps == dp_trailing_gaps::free)
        {
            for (size_t cell_idx = 0; cell_idx < dp_row().size(); ++cell_idx)
                best_score = std::max(dp_row()[cell_idx].score(), best_score);
        }

        if (_column_trailing_gaps == dp_trailing_gaps::free)
        {
            for (size_t cell_idx = 0; cell_idx < dp_column().size(); ++cell_idx)
                best_score = std::max(dp_column()[cell_idx].score(), best_score);
        }

        return best_score;
    }
};

} // namespace _aligner_result

struct result_factory_single
{
    dp_trailing_gaps _column_trailing_gaps{};
    dp_trailing_gaps _row_trailing_gaps{};

    template <std::ranges::viewable_range sequence1_t,
              std::ranges::viewable_range sequence2_t,
              typename dp_column_t,
              typename dp_row_t>
    auto operator()(sequence1_t && sequence1,
                    sequence2_t && sequence2,
                    dp_column_t dp_column,
                    dp_row_t dp_row,
                    dp_trailing_gaps _column_trailing_gaps = dp_trailing_gaps::regular,
                    dp_trailing_gaps _row_trailing_gaps = dp_trailing_gaps::regular) const noexcept
    {
        using aligner_result_t = _aligner_result::value<sequence1_t, sequence2_t, dp_column_t, dp_row_t>;
        return aligner_result_t{std::forward<sequence1_t>(sequence1),
                                std::forward<sequence2_t>(sequence2),
                                std::move(dp_column),
                                std::move(dp_row),
                                _column_trailing_gaps,
                                _row_trailing_gaps};
    }

private:


};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
