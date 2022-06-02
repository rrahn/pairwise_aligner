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

#include <pairwise_aligner/sender/tag_invoke.hpp>

namespace align {

enum class alignment_begin_position {
    first_row_and_column,
    first_row,
    first_column,
    first_row_or_column,
    any
};

class alignment_begin;

namespace _get_alignment_begin {

inline constexpr struct _fn {
    template<typename align_config_t>
        requires tag_invocable<_fn, align_config_t>
    constexpr auto operator()(align_config_t && configuration) const
        noexcept(is_nothrow_tag_invocable_v<_fn, align_config_t>)
        -> tag_invoke_result_t<_fn, align_config_t> {
        return tag_invoke(_fn{}, std::forward<align_config_t>(configuration));
    }

    // default implementation
    template <typename align_config_t>
        requires (!tag_invocable<_fn, align_config_t>)
    constexpr alignment_begin_position operator()(align_config_t &&) const noexcept {
        return alignment_begin_position::first_row_and_column;
    }

} get_alignment_begin{};

} // namespace _get_alignment_begin

using _get_alignment_begin::get_alignment_begin;

} // namespace align
