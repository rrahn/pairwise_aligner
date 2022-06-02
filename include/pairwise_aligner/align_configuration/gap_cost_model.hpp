// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides configurations for alignment.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <concepts>
#include <type_traits>

#include <pairwise_aligner/sender/tag_invoke.hpp>

namespace align {

namespace _get_gap_cost_model {

inline constexpr struct _fn {

    template<typename configuration_t>
        requires tag_invocable<_fn, configuration_t>
    constexpr auto operator()(configuration_t && configuration) const
        noexcept(is_nothrow_tag_invocable_v<_fn, configuration_t>)
        -> tag_invoke_result_t<_fn, configuration_t> {
        return tag_invoke(_fn{}, std::forward<configuration_t>(configuration));
    }

} get_gap_cost_model{};

} // namespace _get_gap_cost_model

using _get_gap_cost_model::get_gap_cost_model;

} // namespace align
