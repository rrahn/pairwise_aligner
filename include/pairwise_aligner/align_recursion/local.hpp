// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides local entry wrapper.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <utility>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

// wrapper for local recursion
template <typename recursion_t>
    // requires multiplication etc.
struct local_recursion : public recursion_t {

    using typename recursion_t::score_type;

    score_type zero{};

    template <typename entry_t, typename ...args_t>
    constexpr friend void tag_invoke(tag_t<align::compute>,
                                     local_recursion & me,
                                     entry_t & entry,
                                     args_t && ...args) {
        using std::max;

        align::compute(static_cast<recursion_t &>(me), entry, std::forward<args_t>(args)...);
        align::current(entry) = max(zero, align::current(entry));
    }

    constexpr friend score_type tag_invoke(tag_t<align::score>, local_recursion const & /*me*/, size_t const) {
        return zero;
    }
};

} // inline namespace v1
} // namespace seqan::pairwise_aligner
