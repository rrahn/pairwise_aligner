// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_matrix_block.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <cassert>

#include <pairwise_aligner/matrix/dp_matrix_block_base.hpp>
#include <pairwise_aligner/matrix/dp_matrix_state_handle.hpp>
// #include <pairwise_aligner/matrix/dp_matrix_lane.hpp>
#include <pairwise_aligner/matrix/dp_matrix_lane_width.hpp>
#include <pairwise_aligner/type_traits.hpp>
// #include <pairwise_aligner/simd/simd_base.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace dp_matrix {

// template <typename lane_closure_t, typename ...dp_data_t>
// struct _block
// {
//     class type;
// };

// template <typename lane_closure_t, typename ...dp_data_t>
// using block_t = typename _block<lane_closure_t, dp_data_t...>::type;

// template <typename lane_closure_t, typename ...dp_data_t>
// class _block<lane_closure_t, dp_data_t...>::type : public detail::dp_matrix_data_handle<dp_data_t...>
// {
//     using base_t = detail::dp_matrix_data_handle<dp_data_t...>;

//     lane_closure_t _lane_closure;

// public:

//     type() = delete;
//     constexpr explicit type(lane_closure_t lane_closure, dp_data_t ...dp_data) noexcept :
//         base_t{std::forward<dp_data_t>(dp_data)...},
//         _lane_closure{std::move(lane_closure)}
//     {
//         // Note the first column/row is not computed again, as they were already initialised.
//         base_t::column()[0].score() = base_t::row()[0].score();
//     }

//     ~type() noexcept
//     {
//         // Store score of last column in first cell of row.
//         base_t::row()[0].score() = base_t::column()[base_t::column().size() - 1].score();
//     }

//     constexpr auto operator[](size_t const index) noexcept
//     {
//         return make_block_lane(*this, index * base_t::lane_width_v, base_t::lane_width(), std::false_type{});
//     }

//     constexpr auto last_lane() noexcept
//     {
//         return make_block_lane(*this,
//                                ((base_t::row().size() - 1) / base_t::lane_width_v) * base_t::lane_width_v,
//                                base_t::lane_width(),
//                                std::true_type{});
//     }

//     constexpr size_t size() const noexcept
//     {
//         return base_t::lanes_per_block();
//     }

// protected:
//     template <typename ...args_t>
//     constexpr auto make_block_lane(args_t && ...args) const noexcept
//         -> std::invoke_result_t<lane_closure_t, args_t...>
//     {
//         return std::invoke(_lane_closure, std::forward<args_t>(args)...);
//     }

// };

// namespace cpo {

// template <typename lane_closure_t = dp_matrix::cpo::_lane_closure>
// struct _block_closure
// {
//     lane_closure_t lane_closure{};

//     template <typename ...dp_data_t>
//     constexpr auto operator()(dp_data_t && ...dp_data) const noexcept {
//         using dp_block_t = dp_matrix::block_t<lane_closure_t, dp_data_t...>;

//         return dp_block_t{lane_closure, std::forward<dp_data_t>(dp_data)...};
//     }
// };

// } // namespace cpo

namespace _block {

template <typename lane_fn_t, typename lane_width_t, typename ...dp_state_t>
class _type : public dp_matrix::detail::block_base<lane_fn_t, lane_width_t, dp_state_t...>
{
    using base_t = dp_matrix::detail::block_base<lane_fn_t, lane_width_t, dp_state_t...>;

public:
    using base_t::base_t;

    constexpr auto column_at(std::ptrdiff_t const index) noexcept
        -> decltype(base_t::make_lane(dp_matrix::detail::static_lane<lane_width_t>,
                                      base_t::lane_offset(index),
                                      base_t::dp_column(),
                                      base_t::dp_row(),
                                      base_t::column_sequence(),
                                      base_t::row_slice_at(index),
                                      base_t::substitution_model(),
                                      base_t::tracker()))
    {
        assert(index < base_t::column_count());

        return base_t::make_lane(dp_matrix::detail::static_lane<lane_width_t>,
                                 base_t::lane_offset(index),
                                 base_t::dp_column(),
                                 base_t::dp_row(),
                                 base_t::column_sequence(),
                                 base_t::row_slice_at(index),
                                 base_t::substitution_model(),
                                 base_t::tracker());
    }

    constexpr auto final_lane() noexcept
        -> decltype(base_t::make_lane((dp_matrix::detail::dynamic_lane<lane_width_t>),
                                      (base_t::lane_offset(base_t::column_count() - 1)),
                                      (base_t::dp_column()),
                                      (base_t::dp_row()),
                                      (base_t::column_sequence()),
                                      (base_t::row_slice_at(0)),
                                      (base_t::substitution_model()),
                                      (base_t::tracker())))
    {
        std::ptrdiff_t const last_index = base_t::column_count() - 1;
        return base_t::make_lane(dp_matrix::detail::dynamic_lane<lane_width_t>,
                                 base_t::lane_offset(last_index),
                                 base_t::dp_column(),
                                 base_t::dp_row(),
                                 base_t::column_sequence(),
                                 base_t::row_slice_at(last_index),
                                 base_t::substitution_model(),
                                 base_t::tracker());
    }
};

struct _fn
{
    template <typename dp_lane_fn_t, size_t lane_width>
    constexpr auto operator()(dp_lane_fn_t && dp_lane_fn, dp_matrix::lane_width_t<lane_width> const &) const noexcept
    {
        std::tuple<dp_lane_fn_t> tmp{std::forward<dp_lane_fn_t>(dp_lane_fn)};
        return [fwd_capture = std::move(tmp)] (auto && ...dp_state) {
            using fwd_dp_lane_fn_t = std::tuple_element_t<0, decltype(fwd_capture)>;
            using block_t = _type<fwd_dp_lane_fn_t,
                                  dp_matrix::lane_width_t<lane_width>,
                                  remove_rvalue_reference_t<decltype(dp_state)>...>;
            return block_t{std::forward<fwd_dp_lane_fn_t &&>(get<0>(fwd_capture)),
                           std::forward<decltype(dp_state)>(dp_state)...};
        };
    }

    template <typename dp_lane_fn_t>
    constexpr auto operator()(dp_lane_fn_t && dp_lane_fn) const noexcept
        -> decltype(std::declval<_fn const &>()(std::declval<dp_lane_fn_t&&>(), lane_width<>))
    {
        return (*this)(std::forward<dp_lane_fn_t>(dp_lane_fn), lane_width<>);
    }
};
} // namespace _block

inline namespace _cpo {
inline constexpr dp_matrix::_block::_fn block{};

} // inline namespace _cpo
} // namespace dp_matrix

// inline constexpr dp_matrix::cpo::_block_closure<> dp_matrix_block{};
} // inline namespace v1
} // namespace seqan::pairwise_aligner
