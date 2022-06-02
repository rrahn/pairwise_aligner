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
namespace _dp_recursion {
namespace _compute {

inline constexpr struct _fn {
    template <typename cost_model_t, typename entry_t, typename score_t>
        requires tag_invocable<_fn, cost_model_t, entry_t, score_t>
    constexpr auto operator()(cost_model_t && cost_model, entry_t && entry, score_t && score) const
        noexcept(is_nothrow_tag_invocable_v<_fn, cost_model_t, entry_t, score_t>)
        -> tag_invoke_result_t<_fn, cost_model_t, entry_t, score_t> {
        return align::tag_invoke(_fn{},
                                 std::forward<cost_model_t>(cost_model),
                                 std::forward<entry_t>(entry),
                                 std::forward<score_t>(score));
    }
} compute;
} // namespace _compute

namespace _score {

inline constexpr struct _fn {
    template <typename cost_model_t, typename ...args_t>
        requires tag_invocable<_fn, cost_model_t, args_t...>
    constexpr auto operator()(cost_model_t && cost_model, args_t && ...args) const
        noexcept(is_nothrow_tag_invocable_v<_fn, cost_model_t, args_t...>)
        -> tag_invoke_result_t<_fn, cost_model_t, args_t...> {
        return align::tag_invoke(_fn{}, std::forward<cost_model_t>(cost_model), std::forward<args_t>(args)...);
    }
} score;
} // namespace _score
} // namespace _dp_recursion

using _dp_recursion::_compute::compute;
using _dp_recursion::_score::score;

} // namespace align
