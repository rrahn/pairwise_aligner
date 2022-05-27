// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides .
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <pairwise_aligner/sender/tag_invoke.hpp>

namespace seqan::align
{
inline namespace v1
{

inline constexpr struct initialise_fn {
    template <typename ...args_t>
    constexpr auto operator()(args_t && ... args)
        noexcept(noexcept(tag_invoke(std::declval<initialise_fn const &>(), std::forward<args_t>(args)...)))
        -> tag_invoke_result_t<initialise_fn, args_t...>
    {
        return tag_invoke(*this, std::forward<args_t>(args)...);
    }
} initialise;

inline constexpr struct entry_at_fn {
    template <typename ...args_t>
    constexpr auto operator()(args_t && ... args)
        noexcept(noexcept(tag_invoke(std::declval<entry_at_fn const &>(), std::forward<args_t>(args)...)))
        -> tag_invoke_result_t<score_fn, args_t...>
    {
        return tag_invoke(*this, std::forward<args_t>(args)...);
    }
} entry_at;

inline constexpr struct result_fn {
    template <typename ...args_t>
    constexpr auto operator()(args_t && ... args)
        noexcept(noexcept(tag_invoke(std::declval<result_fn const &>(), std::forward<args_t>(args)...)))
        -> tag_invoke_result_t<result_fn, args_t...>
    {
        return tag_invoke(*this, std::forward<args_t>(args)...);
    }
} result;

} // inline namespace v1
} // namespace seqan::align
