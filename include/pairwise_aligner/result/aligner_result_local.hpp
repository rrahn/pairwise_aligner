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

#include <pairwise_aligner/configuration/end_gap_policy.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace _aligner_result_local
{
template <typename sequence1_t, typename sequence2_t, typename dp_column_t, typename dp_row_t, typename score_t>
struct _value
{
    struct type;
};

template <typename sequence1_t, typename sequence2_t, typename dp_column_t, typename dp_row_t, typename score_t>
using value = typename _value<sequence1_t, sequence2_t, dp_column_t, dp_row_t, score_t>::type;

template <typename sequence1_t, typename sequence2_t, typename dp_column_t, typename dp_row_t, typename score_t>
struct _value<sequence1_t, sequence2_t, dp_column_t, dp_row_t, score_t>::type
{
    sequence1_t _sequence1;
    sequence2_t _sequence2;
    dp_column_t _dp_column;
    dp_row_t _dp_row;
    score_t _score;

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

    score_t const & score() const noexcept
    {
        return _score;
    }
};

} // namespace _aligner_result_local

template <typename score_t>
struct result_factory_single_local
{
    mutable score_t _best_score{std::numeric_limits<score_t>::lowest()};

    template <std::ranges::viewable_range sequence1_t,
              std::ranges::viewable_range sequence2_t,
              typename dp_column_t,
              typename dp_row_t>
    auto operator()(sequence1_t && sequence1,
                    sequence2_t && sequence2,
                    dp_column_t dp_column,
                    dp_row_t dp_row,
                    [[maybe_unused]] cfg::end_gap _column_trailing_gaps = cfg::end_gap::penalised,
                    [[maybe_unused]] cfg::end_gap _row_trailing_gaps = cfg::end_gap::penalised) const noexcept
    {
        using aligner_result_t = _aligner_result_local::value<sequence1_t, sequence2_t, dp_column_t, dp_row_t, score_t>;
        return aligner_result_t{std::forward<sequence1_t>(sequence1),
                                std::forward<sequence2_t>(sequence2),
                                std::move(dp_column),
                                std::move(dp_row),
                                std::move(_best_score)};
    }

    constexpr score_t const & track_best_score(score_t const & score) const noexcept
    {
        using std::max;
        _best_score = max(_best_score, score);
        return score;
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
