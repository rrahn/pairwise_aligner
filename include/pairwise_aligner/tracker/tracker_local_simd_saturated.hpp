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

template <typename original_score_t, typename score_t>
struct _tracker
{
    class type;
};

template <typename original_score_t, typename score_t>
using tracker = typename _tracker<original_score_t, score_t>::type;

template <typename original_score_t, typename score_t>
class _tracker<original_score_t, score_t>::type
{
private:
    struct _saturated_tracker
    {
        original_score_t _offset{};
        score_t _max_score{std::numeric_limits<typename score_t::value_type>::lowest()};

        score_t const & track(score_t const & score) noexcept {
            using std::max;
            _max_score = max(_max_score, score);
            return score;
        }

        template <typename ...args_t>
        auto max_score([[maybe_unused]] args_t && ...args) const noexcept {
            return original_score_t{_max_score} - _offset;
        }
    };
public:

    original_score_t _max_score{std::numeric_limits<typename original_score_t::value_type>::lowest()};

    auto saturated_tracker(original_score_t const & offset) const noexcept {
        return _saturated_tracker{._offset = offset};
    }

    original_score_t const & track(original_score_t const & score) noexcept {
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

template <typename original_score_t, typename score_t>
struct _factory
{
    struct type;
};

template <typename original_score_t, typename score_t>
using factory = typename _factory<original_score_t, score_t>::type;

template <typename original_score_t, typename score_t>
struct _factory<original_score_t, score_t>::type
{
    constexpr auto make_tracker() const noexcept {
        return tracker<original_score_t, score_t>{};
    }
};

} // namespace tracker::local_simd_saturated

} // inline namespace v1
} // namespace seqan::pairwise_aligner
