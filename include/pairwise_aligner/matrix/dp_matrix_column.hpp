// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_matrix_column.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <cassert>
#include <concepts>
#include <tuple>
#include <utility>

// #include <pairwise_aligner/matrix/dp_matrix_block.hpp>
#include <pairwise_aligner/matrix/dp_matrix_column_base.hpp>
#include <pairwise_aligner/matrix/dp_matrix_state_handle.hpp>
#include <pairwise_aligner/type_traits.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace dp_matrix {

// template <typename block_closure_t, typename ...dp_data_t>
// struct _column
// {
//     class type;
// };

// template <typename block_closure_t, typename ...dp_data_t>
// using column_t = typename _column<block_closure_t, dp_data_t...>::type;

// template <typename block_closure_t, typename ...dp_data_t>
// class _column<block_closure_t, dp_data_t...>::type :
//     public detail::column_base_t<block_closure_t, detail::dp_matrix_data_handle<dp_data_t...>>
// {
// protected:

//     using base_t = detail::column_base_t<block_closure_t, detail::dp_matrix_data_handle<dp_data_t...>>;
// public:
//     using base_t::base_t;
// // public:
// //     type() = delete;
// //     constexpr explicit type(block_closure_t block_closure, dp_data_t ...dp_data) noexcept :
// //         base_t{std::move(block_closure), std::forward<dp_data_t>(dp_data)...}
// //     {}
// };

// namespace cpo {
// template <typename block_closure_t = dp_matrix::cpo::_block_closure<>>
// struct _column_closure
// {
//     block_closure_t block_closure{};

//     template <typename ...dp_data_t>
//     constexpr auto operator()(dp_data_t && ...dp_data) const noexcept {
//         using dp_column_t = dp_matrix::column_t<block_closure_t, dp_data_t...>;

//         return dp_column_t{block_closure, std::forward<dp_data_t>(dp_data)...};
//     }
// };

// } // namespace cpo

namespace _column {

template <typename block_fn_t, typename dp_state_t>
class _type : public dp_matrix::detail::column_base<block_fn_t, dp_state_t>
{
    using base_t = dp_matrix::detail::column_base<block_fn_t, dp_state_t>;

public:
    using base_t::base_t;

    constexpr auto row_at(std::ptrdiff_t const index) noexcept
        -> decltype(std::declval<_type &>().make_matrix_block(base_t::dp_column()[index],
                                                              base_t::dp_row(),
                                                              base_t::column_slice_at(index),
                                                              base_t::row_sequence(),
                                                              base_t::substitution_model(),
                                                              base_t::tracker()))
    {
        assert(index < base_t::row_count());

        return base_t::make_matrix_block(base_t::dp_column()[index],
                                         base_t::dp_row(),
                                         base_t::column_slice_at(index),
                                         base_t::row_sequence(),
                                         base_t::substitution_model(),
                                         base_t::tracker());
    }
};

struct _fn
{
    template <typename dp_block_fn_t>
    constexpr auto operator()(dp_block_fn_t && dp_block_fn) const noexcept
    {
        return [fwd_capture = std::tuple<dp_block_fn_t>(std::forward<dp_block_fn_t>(dp_block_fn))]
                <typename dp_state_t> (dp_state_t && dp_state) {
            using fwd_dp_block_fn_t = std::tuple_element_t<0, decltype(fwd_capture)>;
            // static_assert(std::same_as<fwd_dp_block_fn_t, void>, "fwd_dp_block_fn_t");                         // seqan::pairwise_aligner::v1::dp_matrix::_block::_fn::operator()::._anon_186&&
            // static_assert(std::same_as<decltype(std::get<0>(fwd_capture)), void>, "std::get<0>(fwd_capture)"); // seqan::pairwise_aligner::v1::dp_matrix::_block::_fn::operator()::._anon_186&
            using column_t = _type<fwd_dp_block_fn_t, dp_state_t>;
            return column_t{std::forward<fwd_dp_block_fn_t>(std::get<0>(fwd_capture)), std::move(dp_state)};
        };
    }
};
} // namespace _column


inline namespace _cpo {
inline constexpr dp_matrix::_column::_fn column{};

} // inline namespace _cpo
} // namespace dp_matrix
} // inline namespace v1
}  // namespace seqan::pairwise_aligner
