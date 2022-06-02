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

namespace align {

template <typename cost_model_t>
class zero_cost {
public:
    using score_type = typename std::remove_cvref_t<cost_model_t>::score_type;

    constexpr explicit zero_cost(cost_model_t & cost_model, score_type zero = {}) noexcept :
        _cost_model{cost_model},
        _zero(zero)
    {}

private:
    cost_model_t & _cost_model;
    score_type _zero;

    template <typename entry_t>
    constexpr friend void tag_invoke(tag_t<align::compute>, zero_cost const & me,
                                     entry_t & entry,
                                     score_type substitution_score)
        noexcept(noexcept(align::compute(me._cost_model, entry, std::move(substitution_score)))) {
        using std::max;

        align::compute(me._cost_model, entry, std::move(substitution_score));
        align::current_score(entry) = max(align::current_score(entry), me._zero);
    }

    template <typename ...args_t>
    constexpr friend score_type tag_invoke(tag_t<align::score>, zero_cost const & me, args_t && ...) {
        return me._zero;
    }

    // Make to option
    template <typename model_t>
        requires std::same_as<std::remove_cvref_t<model_t>, zero_cost>
    constexpr friend model_t tag_invoke(tag_t<align::get_gap_cost_model>, model_t && me) noexcept {
        return me;
    }
};
} // namespace align
