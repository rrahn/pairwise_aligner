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
public:

    regular_score_t _max_score{std::numeric_limits<typename regular_score_t::value_type>::lowest()};
    saturated_score_t _block_max_score{std::numeric_limits<typename saturated_score_t::value_type>::lowest()};

    template <typename score_t>
    constexpr score_t const & track(score_t const & score) noexcept {
        using std::max;
        if constexpr (std::same_as<score_t, saturated_score_t>)
            _block_max_score = max(_block_max_score, score);
        else
            _max_score = max(_max_score, score);

        return score;
    }

    constexpr void reset_block_max_score() noexcept
    {
        // TODO: move assignment or better fill?
        _block_max_score = saturated_score_t{std::numeric_limits<typename saturated_score_t::value_type>::lowest()};
    }

    template <typename ...args_t>
    constexpr auto block_max_score([[maybe_unused]] args_t && ...args) const noexcept {
        return _block_max_score;
    }

    template <typename ...args_t>
    constexpr auto max_score([[maybe_unused]] args_t && ...args) const noexcept {
        if constexpr (std::same_as<saturated_score_t, regular_score_t>)
            return _block_max_score;
        else
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
