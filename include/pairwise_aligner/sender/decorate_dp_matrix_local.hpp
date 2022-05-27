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

template <typename entry_t>
struct local_matrix_entry {
    entry_t _entry;

private:

    template <typename ...args_t>
    constexpr friend void tag_invoke(tag_t<compute>, local_matrix_entry & me, args_t && ...args) {
        compute(me._entry, std::forward<args_t>(args)...);

        if (score(me._entry) < 0) {
            score(me._entry) = 0;
        }
    }
};

template <typename matrix_t>
struct dp_matrix_local {

    matrix_t _matrix;

    using value_type = local_matrix_entry<typename std::remove_cvref_t<matrix_t>::value_type>;
    using reference = local_matrix_entry<typename std::remove_cvref_t<matrix_t>::reference>;

    template <typename ...args_t>
    constexpr friend reference tag_invoke(tag_t<entry_at>, local_matrix_wrapper & me, args_t && ...args)
    {
        return reference{entry_at(me._matrix), std::forward<args_t>(args)...};
    }

    template <typename CPO, typename ...args_t>
    constexpr friend auto tag_invoke(CPO cpo, local_matrix_wrapper & me, args_t && ...args)
        -> tag_invoke_result_t<CPO, matrix_t &, args_t...>
    {
        return cpo(me._matrix, std::forward<args_t>(args)...);
    }
};

// What about initialisation?
// that happens in some other form much more lower?
inline constexpr struct decorate_dp_matrix_local_fn {

    template <typename matrix_t> // initialised matrix? so we hook the initialisation function!
    constexpr auto operator()(matrix_t && matrix) {
        return dp_matrix_local<matrix_t>{std::forward<matrix_t>(matrix)};
    }
} local_matrix;

// this can then be selected based on certain local types in the configuration

// local_matrix(dp_matrix())

} // inline namespace v1
} // namespace seqan::align
