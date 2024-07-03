// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::affine_initialisation_strategy.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <pairwise_aligner/affine/affine_cell.hpp>
#include <pairwise_aligner/configuration/end_gap_policy.hpp>
#include <pairwise_aligner/type_traits.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <dp_vector_order order, typename affine_gap_model_t>
struct affine_initialisation_strategy
{
    affine_gap_model_t _gap_model;
    cfg::end_gap _rule;

    template <typename score_t>
    struct _op
    {
        using cell_t = affine_cell<score_t, order>;

        affine_gap_model_t _gap_model;
        cfg::end_gap _rule;

        constexpr cell_t operator()(size_t const index) const noexcept
        {
            // initialise the vertical at index 0:
            // causes the vertical value to have gap_open + gap_extension in first cell
            auto [gap_open, gap_extension] = _gap_model;
            score_t first{0};
            score_t second{static_cast<score_t>(gap_open + gap_extension)};
            if (_rule == cfg::end_gap::penalised && index > 0) {
                first = static_cast<score_t>(gap_open + gap_extension * index);
                second += first;
            }

            return cell_t{first, second};
        }
    };

    template <typename score_t>
    constexpr _op<score_t> create() const noexcept
    {
        return _op<score_t>{_gap_model, _rule};
    }
};
} // inline namespace v1
}  // namespace seqan::pairwise_aligner
