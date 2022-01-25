// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::tracker_local_scalar.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <pairwise_aligner/configuration/end_gap_policy.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

namespace tracker::local_scalar {

template <typename score_t>
struct _tracker
{
    class type;
};

template <typename score_t>
using tracker = typename _tracker<score_t>::type;

template <typename score_t>
class _tracker<score_t>::type
{
public:

    score_t _max_score{std::numeric_limits<score_t>::lowest()};

    score_t const & track(score_t const & score) noexcept {
        using std::max;
        _max_score = max(_max_score, score);
        return score;
    }

    template <typename ...args_t>
    auto max_score([[maybe_unused]] args_t && ...args) const noexcept {
        return _max_score;
    }

    // TODO: optimal_coordinate()
};

// this is the capture object which creates a new instance of a tracker every time we call it.
// the client can provide everything here.
// I need the factory here but some value to create it.
template <typename score_t>
struct _factory
{
    struct type;
};

template <typename score_t>
using factory = typename _factory<score_t>::type;

template <typename score_t>
struct _factory<score_t>::type
{
    constexpr auto make_tracker() const noexcept {
        return tracker<score_t>{};
    }
};

} // namespace tracker::local_scalar

} // inline namespace v1
} // namespace seqan::pairwise_aligner
