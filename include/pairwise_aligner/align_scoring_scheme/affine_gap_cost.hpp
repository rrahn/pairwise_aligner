// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides affine gap cost model.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <pairwise_aligner/sender/tag_invoke.hpp>
#include <pairwise_aligner/align_matrix/dp_entry_concept.hpp>
#include <pairwise_aligner/align_matrix/dp_recursion_concept.hpp>
#include <pairwise_aligner/align_configuration/gap_cost_model.hpp>

namespace align {

template <typename score_t>
class affine_cost_model {
private:
    score_t _extension_cost{-1};
    score_t _open_cost{};

public:

    using score_type = score_t;

    affine_cost_model() = default;
    constexpr affine_cost_model(score_t open_cost, score_t extension_cost) noexcept :
        _extension_cost{extension_cost},
        _open_cost{open_cost + extension_cost}
    {}

private:

    template <typename entry_t>
        // requires align::dp_entry<entry_t>
    constexpr friend void tag_invoke(tag_t<align::compute>, affine_cost_model const & me,
                                     entry_t && entry,
                                     score_t substitution_score) noexcept {
        using std::max;

        substitution_score += align::diagonal_score(entry);
        score_type best = max(max(substitution_score, align::left_score(entry)), align::left_score(entry));
        align::diagonal_score(entry) = align::current_score(entry); // cache next diagonal score!
        align::current_score(entry) = best;
        // get<0>(column_cell) = tracker.track(best); // get<0>(column_cell) = best;
        best += me._open_cost;
        align::up_score(entry) = max(align::up_score(entry) + me._extension_cost, best);
        align::left_score(entry) = max(align::left_score(entry) + me._extension_cost, best);
    }

    constexpr friend score_t tag_invoke(tag_t<align::score>, affine_cost_model const & me, size_t const gap_length)
        noexcept {
        return (gap_length > 0) * (me._open_cost + (gap_length - 1) * me._extension_cost);
    }

    // Make to option
    template <typename cost_model_t>
        requires std::same_as<std::remove_cvref_t<cost_model_t>, affine_cost_model>
    constexpr friend cost_model_t tag_invoke(tag_t<align::get_gap_cost_model>, cost_model_t && me) noexcept {
        return me;
    }
};
} // namespace align
