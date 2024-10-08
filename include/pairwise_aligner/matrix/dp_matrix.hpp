// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_matrix_policies.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <cassert>
#include <tuple>

#include <pairwise_aligner/matrix/dp_matrix_state_handle.hpp>
#include <pairwise_aligner/type_traits.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename dp_matrix_column_closure_t>
struct _dp_matrix_policies
{
    struct type;
};

template <typename dp_matrix_column_closure_t>
using dp_matrix_policies = typename _dp_matrix_policies<dp_matrix_column_closure_t>::type;

template <typename dp_matrix_column_closure_t>
struct _dp_matrix_policies<dp_matrix_column_closure_t>::type
{
    dp_matrix_column_closure_t dp_matrix_column_closure;

    template <typename ...args_t>
    constexpr auto make_policies(args_t && ...args) const noexcept {
        return std::invoke(dp_matrix_column_closure, dp_matrix::detail::make_dp_state(std::forward<args_t>(args)...));
    }
};

namespace dp_matrix {

namespace _dp_matrix {
template <typename dp_column_fn_t, typename dp_state_t>
class _type : public dp_state_t
{
private:

    using base_t = dp_state_t;

    dp_column_fn_t _dp_column_fn;

public:

    _type() = default;
    _type(dp_column_fn_t dp_column_fn, dp_state_t && dp_state) noexcept :
        base_t{std::move(dp_state)},
        _dp_column_fn{std::forward<dp_column_fn_t>(dp_column_fn)}
    {}

    constexpr auto column_at(std::ptrdiff_t const index) noexcept
        -> std::invoke_result_t<dp_column_fn_t,
                                decltype(dp_matrix::detail::make_dp_state(
                                            base_t::dp_column(),
                                            base_t::dp_row()[0],
                                            base_t::column_sequence(),
                                            seqan3::views::slice(base_t::row_sequence(), 0, 1),
                                            base_t::substitution_model(),
                                            base_t::tracker()
                                        ))>
    {
        assert(index < base_t::column_count());

        std::ptrdiff_t const row_chunk_size = base_t::dp_row()[0].size() - 1;
        std::ptrdiff_t const row_offset = row_chunk_size * index;
        return std::invoke(_dp_column_fn, dp_matrix::detail::make_dp_state(
                            base_t::dp_column(),
                            base_t::dp_row()[index],
                            base_t::column_sequence(),
                            seqan3::views::slice(base_t::row_sequence(), row_offset, row_offset + row_chunk_size),
                            base_t::substitution_model(),
                            base_t::tracker()
                        ));
    }
};

// adaptor closure object
struct _fn
{
    template <typename dp_column_fn_t>
    constexpr auto operator()(dp_column_fn_t dp_column_fn) const noexcept
    {
        return [dp_column_fn = std::move(dp_column_fn)] <typename dp_state_t> (dp_state_t && dp_state) {
            using dp_matrix_t = _type<dp_column_fn_t, dp_state_t>;
            return dp_matrix_t{std::move(dp_column_fn), std::move(dp_state)};
        };
    }
};

} // namespace _dp_matrix
inline namespace _cpo {

inline constexpr _dp_matrix::_fn matrix{};

} // inline namespace _cpo
} // namespace dp_matrix
} // inline namespace v1
}  // namespace seqan::pairwise_aligner
