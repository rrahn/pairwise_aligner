// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides configure CPO.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <pairwise_aligner/sender/tag_invoke.hpp>

namespace align
{

namespace _configure {

inline const struct _fn {

    template<typename object_t>
        requires tag_invocable<_fn, object_t>
    constexpr auto operator()(object_t && object) const
        noexcept(is_nothrow_tag_invocable_v<_fn, object_t>)
        -> tag_invoke_result_t<_fn, object_t> {
        return align::tag_invoke(_fn{}, std::forward<object_t>(object));
    }

} configure{};

} // namespace _configure

using _configure::configure;

} // namespace align
