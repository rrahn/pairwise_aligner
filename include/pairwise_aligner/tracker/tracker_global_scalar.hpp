// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::tracker::global_scalar.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <pairwise_aligner/configuration/end_gap_policy.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

namespace tracker::global_scalar {

class tracker
{
public:
    cfg::trailing_end_gap _end_gap;

    template <typename score_t>
    score_t const & track(score_t const & score) const noexcept {
        return score; // no-op.
    }

    template <typename sequence1_t, typename sequence2_t, typename dp_column_t, typename dp_row_t>
    constexpr auto max_score([[maybe_unused]] sequence1_t && sequence1,
                             [[maybe_unused]] sequence2_t && sequence2,
                             dp_column_t const & dp_column,
                             dp_row_t const & dp_row) const noexcept
    {
        assert(dp_column.size() == 1);
        assert(dp_row.size() == 1);

        size_t const inner_size = dp_column[0].size();
        auto best_score = dp_column[0][inner_size - 1].score();

        if (_end_gap.last_row == cfg::end_gap::free)
        {
            for (size_t cell_idx = 0; cell_idx < dp_row[0].size(); ++cell_idx)
                best_score = std::max(dp_row[0][cell_idx].score(), best_score);
        }

        if (_end_gap.last_column == cfg::end_gap::free)
        {
            for (size_t cell_idx = 0; cell_idx < dp_column[0].size(); ++cell_idx)
                best_score = std::max(dp_column[0][cell_idx].score(), best_score);
        }

        return best_score;
    }

    // TODO: optimal_coordinate()
};

// this is the capture object which creates a new instance of a tracker every time we call it.
// the client can provide everything here.
// I need the factory here but some value to create it.
struct factory
{
    // params for free end-gaps.
    cfg::trailing_end_gap _end_gap{};

    constexpr auto make_tracker() const noexcept {
        return tracker{_end_gap};
    }
};

} // namespace tracker::global_scalar

} // inline namespace v1
} // namespace seqan::pairwise_aligner
