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

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace dp_matrix {

template <typename ...dp_data_t>
struct _block
{
    class type;
};

template <typename ...dp_data_t>
using block_t = typename _block<dp_data_t...>::type;

template <typename ...dp_data_t>
class _block<dp_data_t...>::type : public detail::dp_matrix_data_handle<dp_data_t...>
{
    static constexpr size_t _lane_width = 8; //TODO: make copnfigurable?

    using base_t = detail::dp_matrix_data_handle<dp_data_t...>;

public:

    type() = delete;
    constexpr explicit type(dp_data_t ...dp_data) noexcept : base_t{std::forward<dp_data_t>(dp_data)...}
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
        return dp_matrix_lane(*this, index * lane_width(), lane_width_t<type::lane_width()>{});
    }

    constexpr auto last_lane() noexcept
    {
        return dp_matrix_last_lane(*this,
                                   ((base_t::row().size() - 1) / lane_width()) * lane_width(),
                                   lane_width_t<type::lane_width()>{});
    }

    static constexpr size_t lane_width() noexcept
    {
        return _lane_width;
    }

    constexpr size_t size() const noexcept
    {
        return (base_t::row().size() - 1 + lane_width()) / lane_width();
    }
};

namespace cpo {

struct _block_closure
{
    template <typename ...dp_data_t>
    constexpr auto operator()(dp_data_t && ...dp_data) const noexcept {
        using dp_block_t = dp_matrix::block_t<dp_data_t...>;

        return dp_block_t{std::forward<dp_data_t>(dp_data)...};
    }
};

} // namespace cpo
} // namespace dp_matrxix

inline constexpr dp_matrix::cpo::_block_closure dp_matrix_block{};
} // inline namespace v1
} // namespace seqan::pairwise_aligner
