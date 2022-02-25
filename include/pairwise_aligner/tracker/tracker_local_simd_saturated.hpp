// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::tracker_local_simd_saturated.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <pairwise_aligner/tracker/tracker_local_simd_fixed.hpp>
namespace seqan::pairwise_aligner
{
inline namespace v1
{

namespace tracker::local_simd_saturated {

template <typename saturated_score_t, typename regular_score_t>
struct _tracker
{
    class type;
};

template <typename saturated_score_t, typename regular_score_t>
using tracker = typename _tracker<saturated_score_t, regular_score_t>::type;

template <typename saturated_score_t, typename regular_score_t>
class _tracker<saturated_score_t, regular_score_t>::type
{
private:

    regular_score_t _max_score{std::numeric_limits<typename regular_score_t::value_type>::lowest()};

    struct _in_block_tracker : public local_simd_fixed::tracker<saturated_score_t>
    {
    private:
        using base_t = local_simd_fixed::tracker<saturated_score_t>;

        type & _parent_tracker;
        regular_score_t const & _score_offset;

    public:
        _in_block_tracker() = delete;
        _in_block_tracker(type & parent_tracker, regular_score_t const & score_offset) noexcept :
            base_t{},
            _parent_tracker{parent_tracker},
            _score_offset{score_offset}
        {} // initialising the local tracker
        ~_in_block_tracker() noexcept
        {
            _parent_tracker.track(regular_score_t{base_t::max_score()} + _score_offset);
        }
    };

public:

    constexpr auto in_block_tracker(regular_score_t const & score_offset) noexcept {
        return _in_block_tracker{*this, score_offset};
    }

    constexpr void track(regular_score_t const & in_block_score) noexcept {
        _max_score = max(_max_score, in_block_score);
    }

    template <typename ...args_t>
    constexpr auto max_score([[maybe_unused]] args_t && ...args) const noexcept {
        return _max_score;
    }

    // TODO: optimal_coordinate()
};

// this is the capture object which creates a new instance of a tracker every time we call it.
// the client can provide everything here.
// I need the factory here but some value to create it.
template <typename saturated_score_t, typename regular_score_t>
struct _factory
{
    struct type;
};

template <typename saturated_score_t, typename regular_score_t>
using factory = typename _factory<saturated_score_t, regular_score_t>::type;

template <typename saturated_score_t, typename regular_score_t>
struct _factory<saturated_score_t, regular_score_t>::type
{
    constexpr auto make_tracker() const noexcept {
        return tracker<saturated_score_t, regular_score_t>{};
    }
};

} // namespace tracker::local_simd_saturated

} // inline namespace v1
} // namespace seqan::pairwise_aligner
