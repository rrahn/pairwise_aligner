// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::cfg::detail::saturated_block_handler.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <cassert>
#include <cmath>

namespace seqan::pairwise_aligner {
inline namespace v1
{
namespace cfg
{
namespace detail
{

struct saturated_block_handler
{
public:

    template <typename score_t>
    static auto lowest_viable_local_score(score_t const gap_open, score_t const gap_extension) noexcept
    {
        return std::numeric_limits<int8_t>::lowest() - gap_open - (2 * gap_extension);
    }

    template <typename score_t, typename gap_score_t>
    static constexpr auto compute_max_block_size(score_t const match,
                                                 score_t const mismatch,
                                                 gap_score_t const gap_open,
                                                 gap_score_t const gap_extension) noexcept
    {
        auto [block_size_gap, zero_offset_gap] = max_block_size_with_gaps(match,
                                                                          std::abs(gap_open),
                                                                          std::abs(gap_extension));
        auto [block_size_mismatch, zero_offset_mismatch] = max_block_size_with_mismatch(match, std::abs(mismatch));

        // Choose the max of both block sizes since due to the recursion the largest negative distance is affected
        // by the lowest negative score with the given block size.
        // Also set the corresponding zero offset accordingly.
        int8_t zero_offset{};
        size_t max_block_size = (block_size_gap < block_size_mismatch)
                              ? (zero_offset = zero_offset_mismatch, block_size_mismatch)
                              : (zero_offset = zero_offset_gap, block_size_gap);
        return std::pair{zero_offset, max_block_size};
    }

private:
    static constexpr auto max_block_size_with_gaps(float const match,
                                                   float const gap_open,
                                                   float const gap_extension) noexcept
    {
        // Assume positive gap scores for the following formula.
        assert(gap_open >= 0);
        assert(gap_extension >= 0);

        // This formula computes the maximal block size and the respective zero offset for which the scores during
        // the computation of the alignment blocks in saturated mode are guaranteed to not exceed the value range
        // of the machines int8_t type, under the condition of having only gaps in a single block.
        // It is derived from solving the following equations:
        //  I: 127 >= x + match * block_size + gap_open + gap_extension * block_size
        // II: -128 <= x - 2 * (gap_open + gap_extension * block_size)
        // Here x represents the zero offset for which the block size is maximised.

        float const max_score = std::numeric_limits<int8_t>::max();
        float const min_score = std::numeric_limits<int8_t>::lowest();
        float const upper_score_limit = max_score - gap_open;
        float const block_divisor = match + gap_extension;
        float const block_scale = 2 * gap_extension / block_divisor;

        int8_t const zero_offset = std::ceil((min_score + 2 * gap_open + (block_scale * upper_score_limit)) /
                                             (1 + block_scale));
        size_t const block_size = std::floor((upper_score_limit - zero_offset) / block_divisor);

        return std::pair{block_size, zero_offset};
    }

    static constexpr auto max_block_size_with_mismatch(float const match, float const mismatch) noexcept
    {
        // Assume positive mismatch score for the following equation
        assert(mismatch >= 0);

        // This formula computes the maximal block size and the respective zero offset for which the scores during
        // the computation of the alignment blocks in saturated mode are guaranteed to not exceed the value range
        // of the machines int8_t type, under the condition of having only mismatches in a single block.
        // It is derived from solving the following equations:
        //  I: 127 >= x + match * block_size + mismatch * block_size
        // II: -128 <= x - mismatch * block_size
        // Here x represents the zero offset for which the block size is maximised.

        float const max_score = std::numeric_limits<int8_t>::max();
        float const min_score = std::numeric_limits<int8_t>::lowest();
        float const block_scale = mismatch / (match + mismatch);

        int8_t const zero_offset = std::ceil((min_score + max_score * block_scale) / (1 + block_scale));
        size_t const block_size = std::floor((max_score - zero_offset) / (match + mismatch));

        return std::pair{block_size, zero_offset};
    }
};
} // namespace detail
} // namespace cfg
} // inline namespace v1
} // namespace seqan::pairwise_aligner
