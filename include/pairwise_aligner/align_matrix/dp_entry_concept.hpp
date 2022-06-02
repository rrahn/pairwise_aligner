// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides concepts and CPOs related to the DP matrix computation.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <concepts>
#include <type_traits>

#include <pairwise_aligner/sender/tag_invoke.hpp>

namespace align
{

// Accessor functions to DP entry
namespace _dp_entry {

template <typename>
struct _fn {
    template <typename entry_t>
        requires tag_invocable<_fn, entry_t>
    constexpr auto operator()(entry_t && entry) const
        noexcept(is_nothrow_tag_invocable_v<_fn, entry_t>)
        -> tag_invoke_result_t<_fn, entry_t> {
        return align::tag_invoke(_fn{}, std::forward<entry_t>(entry));
    }
};

inline constexpr _fn<struct _current> current_score;
inline constexpr _fn<struct _diagonal> diagonal_score;
inline constexpr _fn<struct _up> up_score;
inline constexpr _fn<struct _left> left_score;
} // namespace _dp_entry

using _dp_entry::current_score;
using _dp_entry::diagonal_score;
using _dp_entry::up_score;
using _dp_entry::left_score;

template <typename t>
concept dp_entry = requires (t & e) {

    typename std::remove_cvref_t<t>::score_type;

    { align::current_score(e) } -> std::same_as<typename std::remove_cvref_t<t>::score_type &>;
    { align::diagonal_score(e) } -> std::same_as<typename std::remove_cvref_t<t>::score_type &>;
    { align::up_score(e) } -> std::same_as<typename std::remove_cvref_t<t>::score_type &>;
    { align::left_score(e) } -> std::same_as<typename std::remove_cvref_t<t>::score_type &>;
};

} // namespace align
