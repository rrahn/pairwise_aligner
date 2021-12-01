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

#include <pairwise_aligner/dp_initialisation_rule.hpp>
#include <pairwise_aligner/type_traits.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename affine_gap_model_t>
class affine_initialisation_strategy
{
private:
    affine_gap_model_t _gap_model{};
    dp_initialisation_rule _row_initialisation_rule{};
    dp_initialisation_rule _column_initialisation_rule{};
    size_t _row_index{};
    size_t _column_index{};

public:

    affine_initialisation_strategy() = default;
    explicit affine_initialisation_strategy(
        affine_gap_model_t affine_gap_model,
        dp_initialisation_rule const row_initialisation_rule = dp_initialisation_rule::regular,
        dp_initialisation_rule const column_initialisation_rule = dp_initialisation_rule::regular) noexcept :

        _gap_model{std::move(affine_gap_model)},
        _row_initialisation_rule{row_initialisation_rule},
        _column_initialisation_rule{column_initialisation_rule}
    {}

    template <typename affine_cell_t>
    void operator()(affine_cell_t & affine_cell)
    {
        if constexpr (is_row_cell_v<affine_cell_t>) {
            initialise_cell(affine_cell, _row_index, _row_initialisation_rule);
        } else {
            initialise_cell(affine_cell, _column_index, _column_initialisation_rule);
        }
    }

private:

    template <typename affine_cell_t>
    constexpr void initialise_cell(affine_cell_t & affine_cell,
                                   size_t & index,
                                   dp_initialisation_rule const initialisation_rule) const
    {
        using score_t = affine_cell_t::score_type;

        switch (initialisation_rule) {
            case dp_initialisation_rule::regular: {
                if (index == 0) {
                    affine_cell = affine_cell_t{score_t{0}, score_t{0}};
                } else {
                    auto [gap_open, gap_extension] = _gap_model;
                    auto first = gap_open + index * gap_extension;
                    auto second = gap_open + index * gap_extension + gap_open + gap_extension;
                    affine_cell = affine_cell_t{static_cast<score_t>(first), static_cast<score_t>(second)};
                }
                break;
            } case dp_initialisation_rule::zero: {
                affine_cell = affine_cell_t{score_t{0}, score_t{0}};
                break;
            }
            default: throw std::invalid_argument{"Unknown rule to initialise affine cell."};
        }
        ++index;
    }
};
} // inline namespace v1
}  // namespace seqan::pairwise_aligner
