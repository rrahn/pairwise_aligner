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

#include <pairwise_aligner/matrix/dp_matrix_data_handle.hpp>
#include <pairwise_aligner/matrix/dp_matrix_lane.hpp>
#include <pairwise_aligner/simd/simd_base.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace dp_matrix {

template <typename lane_closure_t, size_t _lane_width, typename ...dp_data_t>
struct _block
{
    class type;
};

template <typename lane_closure_t, size_t _lane_width, typename ...dp_data_t>
using block_t = typename _block<lane_closure_t, _lane_width, dp_data_t...>::type;

template <typename lane_closure_t, size_t _lane_width, typename ...dp_data_t>
class _block<lane_closure_t, _lane_width, dp_data_t...>::type : public detail::dp_matrix_data_handle<dp_data_t...>
{
    using base_t = detail::dp_matrix_data_handle<dp_data_t...>;

    lane_closure_t _lane_closure;

public:

    type() = delete;
    constexpr explicit type(lane_closure_t lane_closure, dp_data_t ...dp_data) noexcept :
        base_t{std::forward<dp_data_t>(dp_data)...},
        _lane_closure{std::move(lane_closure)}
    {
        // Note the first column/row is not computed again, as they were already initialised.
        base_t::column()[0].score() = base_t::row()[0].score();
    }

    ~type() noexcept
    {
        // Store score of last column in first cell of row.
        base_t::row()[0].score() = base_t::column()[base_t::column().size() - 1].score();
    }

    constexpr auto operator[](size_t const index) noexcept
    {
        return make_block_lane(*this, index * lane_width(), lane_width_t<type::lane_width()>{}, std::false_type{});
    }

    constexpr auto last_lane() noexcept
    {
        return make_block_lane(*this,
                               ((base_t::row().size() - 1) / lane_width()) * lane_width(),
                               lane_width_t<type::lane_width()>{},
                               std::true_type{});
    }

    static constexpr size_t lane_width() noexcept
    {
        return _lane_width;
    }

    constexpr size_t size() const noexcept
    {
        return (base_t::row().size() - 1 + lane_width()) / lane_width();
    }
protected:
    template <typename ...args_t>
    constexpr auto make_block_lane(args_t && ...args) const noexcept
        -> std::invoke_result_t<lane_closure_t, args_t...>
    {
        return std::invoke(_lane_closure, std::forward<args_t>(args)...);
    }

};

namespace cpo {

template <typename lane_closure_t = dp_matrix::cpo::_lane_closure,
          size_t lane_width = ((seqan::pairwise_aligner::detail::max_simd_size == 64) ? 8 : 4)>
struct _block_closure
{
    lane_closure_t lane_closure{};

    template <typename ...dp_data_t>
    constexpr auto operator()(dp_data_t && ...dp_data) const noexcept {
        using dp_block_t = dp_matrix::block_t<lane_closure_t, lane_width, dp_data_t...>;

        return dp_block_t{lane_closure, std::forward<dp_data_t>(dp_data)...};
    }
};

} // namespace cpo
} // namespace dp_matrxix

inline constexpr dp_matrix::cpo::_block_closure<> dp_matrix_block{};
} // inline namespace v1
} // namespace seqan::pairwise_aligner
