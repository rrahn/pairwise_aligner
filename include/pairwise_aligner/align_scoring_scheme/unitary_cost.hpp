// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides unitary cost model.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <concepts>

#include <pairwise_aligner/sender/tag_invoke.hpp>
#include <pairwise_aligner/align_matrix/dp_recursion_concept.hpp>
#include <pairwise_aligner/align_configuration/substitution_cost_model.hpp>

namespace align {

template <typename score_t>
class unitary_cost_model {
private:
    score_t _match_cost{0};
    score_t _mismatch_cost{-1};

public:

    using score_type = score_t;

    unitary_cost_model() = default;
    constexpr unitary_cost_model(score_t match_cost, score_t mismatch_cost) noexcept :
        _match_cost{match_cost},
        _mismatch_cost{mismatch_cost}
    {}

private:

    template <typename lhs_t, typename rhs_t>
        requires std::equality_comparable_with<lhs_t const &, rhs_t const &>
    constexpr friend score_t tag_invoke(tag_t<align::score>,
                                        unitary_cost_model const & me,
                                        lhs_t const & lhs,
                                        rhs_t const & rhs)
        noexcept {
        return (lhs == rhs) ? me._match_cost : me._mismatch_cost;
    }

    // Make to option
    template <typename cost_model_t>
        requires std::same_as<std::remove_cvref_t<cost_model_t>, unitary_cost_model>
    constexpr friend cost_model_t tag_invoke(tag_t<align::get_substitution_cost_model>, cost_model_t && me) noexcept {
        return me;
    }
};
} // namespace align
