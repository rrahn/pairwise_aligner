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

#include <pairwise_aligner/matrix/dp_matrix_state_handle.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace dp_matrix::detail {

template <typename lane_fn_t, typename lane_width_t, typename ...dp_state_t>
class block_base : public state_handle<dp_state_t...>
{
    using base_t = state_handle<dp_state_t...>;

    lane_fn_t _lane_fn;

public:

    static constexpr std::ptrdiff_t lane_width = lane_width_t::value;

    block_base() = delete;
    constexpr explicit block_base(lane_fn_t lane_fn, dp_state_t ...dp_state) noexcept :
        base_t{std::forward<dp_state_t>(dp_state)...},
        _lane_fn{std::move(lane_fn)}
    {
        // Note the first column/row is not computed again, as they were already initialised.
        base_t::dp_column()[0].score() = base_t::dp_row()[0].score();
    }

    ~block_base() noexcept
    {
        // Store score of last column in first cell of row.
        base_t::dp_row()[0].score() = base_t::dp_column()[base_t::dp_column().size() - 1].score();
    }

    constexpr std::ptrdiff_t column_count() const noexcept
    {
        return (base_t::column_count() - 2 + lane_width) / lane_width;
    }

protected:

    constexpr auto row_slice_at(std::ptrdiff_t const index) const noexcept
    {
        std::ptrdiff_t const row_offset = lane_width * std::max<std::ptrdiff_t>(0, index);
        return seqan3::views::slice(base_t::row_sequence(), row_offset, row_offset + lane_width);
    }

    template <typename ...args_t>
    constexpr auto make_lane(args_t && ...args) const noexcept
        -> std::invoke_result_t<lane_fn_t, args_t...>
    {
         return std::invoke(_lane_fn, std::forward<args_t>(args)...);
    }

    constexpr std::ptrdiff_t lane_offset(std::ptrdiff_t const index) const noexcept
    {
        return index * lane_width;
    }
};

} // namespace dp_matrxix::detail
} // inline namespace v1
} // namespace seqan::pairwise_aligner
