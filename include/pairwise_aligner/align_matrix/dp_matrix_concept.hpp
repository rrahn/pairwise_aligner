// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides align::.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <pairwise_aligner/sender/tag_invoke.hpp>

namespace align
{
namespace _dp_matrix {
namespace _entry_at {

inline constexpr struct _fn {
    template <typename matrix_t, typename matrix_position_t>
        requires tag_invocable<_fn, matrix_t, matrix_position_t const &, matrix_position_t const &>
    constexpr auto operator()(matrix_t && matrix,
                              matrix_position_t const & row_position,
                              matrix_position_t const & column_position) const
        noexcept(is_nothrow_tag_invocable_v<_fn, matrix_t, matrix_position_t const &, matrix_position_t const &>)
        -> tag_invoke_result_t<_fn, matrix_t, matrix_position_t const &, matrix_position_t const &> {
        return align::tag_invoke(_fn{}, std::forward<matrix_t>(matrix), row_position, column_position);
    }
} entry_at;
} // namespace _entry_at

namespace _initialise_vector {

template <typename>
struct _fn {
    template <typename matrix_t, typename gap_cost_model_t>
        requires tag_invocable<_fn, matrix_t, size_t, gap_cost_model_t>
    constexpr auto operator()(matrix_t && matrix, size_t const dimension, gap_cost_model_t && fn) const
        noexcept(is_nothrow_tag_invocable_v<_fn, matrix_t, size_t, gap_cost_model_t>)
        -> tag_invoke_result_t<_fn, matrix_t, size_t, gap_cost_model_t> {
        return align::tag_invoke(_fn{}, std::forward<matrix_t>(matrix), dimension, std::forward<gap_cost_model_t>(fn));
    }
};
inline constexpr _fn<struct _row> initialise_row;
inline constexpr _fn<struct _column> initialise_column;
} // namespace _initialise_vector
} // namespace _dp_matrix

using _dp_matrix::_entry_at::entry_at;
using _dp_matrix::_initialise_vector::initialise_row;
using _dp_matrix::_initialise_vector::initialise_column;

} // namespace align
