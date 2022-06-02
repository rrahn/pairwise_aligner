// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides affine recursion.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <type_traits>
#include <utility>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename score_t>
    // requires multiplication etc.
struct affine_recursion {

    typename score_type = score_t;

    score_t open_cost;
    score_t extension_cost;

private:

    template <typename entry_t>
    constexpr friend void tag_invoke(tag_t<align::compute>,
                                     affine_recursion & me,
                                     entry_t & entry,
                                     score_t substitution_score) {
        using std::max;

        substitution_score += align::diagonal(entry);
        best = max(max(substitution_score, align::left(up)), align::left(entry));
        align::diagonal(entry) = align::current(entry); // cache next diagonal score!
        align::current(entry) = best;
        // get<0>(column_cell) = tracker.track(best); // get<0>(column_cell) = best;
        best += (me.open_cost + me.extension_cost);
        align::up(entry) = max(align::up(entry) + me.extension_cost, best);
        align::left(entry) = max(align::left(entry) + me.extension_cost, best);
    }

    constexpr friend score_t tag_invoke(tag_t<align::score>, affine_recursion const & me, size_t const gap_length) {
        return (gap_length > 0) * me.open_cost + gap_length * me.extension_cost;
    }
};

} // inline namespace v1
} // namespace seqan::pairwise_aligner
