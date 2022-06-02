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

#include <unifex/tag_invoke.hpp>

namespace align {
inline namespace v1 {
using namespace unifex::_tag_invoke_cpo;

template <auto & cpo>
using tag_t = unifex::tag_t<cpo>;

template <typename cpo_t, typename ...args_t>
using tag_invoke_result_t = unifex::tag_invoke_result_t<cpo_t, args_t...>;

template <typename cpo_t, typename... args_t>
using is_tag_invocable = unifex::is_tag_invocable<cpo_t, args_t...>;

template <typename cpo_t, typename ...args_t>
inline constexpr bool is_tag_invocable_v = unifex::is_tag_invocable_v<cpo_t, args_t...>;

template <typename cpo_t, typename ...args_t>
inline constexpr bool is_nothrow_tag_invocable_v = unifex::is_nothrow_tag_invocable_v<cpo_t, args_t...>;

template <typename cpo_t, typename ...args_t>
concept tag_invocable = unifex::tag_invocable<cpo_t, args_t...>;

using unifex::tag_invoke;

} // inline namespace v1
} // namespace align
