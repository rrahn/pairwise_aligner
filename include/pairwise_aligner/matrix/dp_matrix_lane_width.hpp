// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_matrix::lane_width.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <type_traits>

#include <pairwise_aligner/simd/simd_base.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace dp_matrix {

template <std::size_t width>
using lane_width_t = std::integral_constant<std::size_t, width>;

template <std::size_t width = (seqan::pairwise_aligner::detail::max_simd_size == 64 ? 8 : 4)>
inline constexpr lane_width_t<width> lane_width;

namespace detail {

template <typename width_t, bool is_dynamic>
struct lane_tag : public std::bool_constant<is_dynamic>
{
    static constexpr size_t width = width_t::value;
};

template <typename width_t>
inline constexpr lane_tag<width_t, false> static_lane{};

template <typename width_t>
inline constexpr lane_tag<width_t, true> dynamic_lane{};

} // mamespace detail
} // namespace dp_matrix

template <size_t width = ((seqan::pairwise_aligner::detail::max_simd_size == 64) ? 8 : 4)>
struct lane_width_policy
{
    constexpr dp_matrix::lane_width_t<width> make_lane_width() const noexcept
    {
        return dp_matrix::lane_width<width>;
    }
};

} // inline namespace v1
} // namespace seqan::pairwise_aligner
