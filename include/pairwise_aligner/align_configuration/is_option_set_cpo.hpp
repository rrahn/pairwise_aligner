// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides has_option CPO.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <pairwise_aligner/sender/tag_invoke.hpp>

namespace align
{

namespace _is_option_set {

inline const struct _fn {

    // template<typename option_t>
    //     requires tag_invocable<_fn, option_t const &, option_t const &>
    // constexpr auto operator()(option_t const & option, option_t const & rhs_option) const
    //     noexcept(is_nothrow_tag_invocable_v<_fn, option_t const &, option_t const &>)
    //     -> tag_invoke_result_t<_fn, option_t const &, option_t const &> {
    //     return tag_invoke(_fn{}, option, rhs_option);
    // }

    template<typename option_t>
        requires tag_invocable<_fn, option_t const &, typename option_t::value_type const &>
    constexpr auto operator()(option_t const & option, typename option_t::value_type const & value) const
        noexcept(is_nothrow_tag_invocable_v<_fn, option_t const &, typename option_t::value_type const &>)
        -> tag_invoke_result_t<_fn, option_t const &, typename option_t::value_type const &> {
        return tag_invoke(_fn{}, option, value);
    }

    // // default implementation returns always false.
    // template<typename option_t>
    //     requires (!tag_invocable<_fn, option_t const &, option_t const &>)
    // constexpr bool operator()(option_t const &, option_t const &) const noexcept {
    //     return false;
    // }

} is_option_set{};

} // namespace _is_option_set

using _is_option_set::is_option_set;

} // namespace align
