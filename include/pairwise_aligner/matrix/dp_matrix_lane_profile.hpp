// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_matrix_lane_profile.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <seqan3/utility/views/slice.hpp>

#include <pairwise_aligner/matrix/dp_matrix_lane.hpp>
#include <pairwise_aligner/score_model/strip_width.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace dp_matrix {

template <typename dp_block_t, size_t lane_width, bool is_last_lane>
struct _lane_profile
{
    class type;
};

template <typename dp_block_t, size_t lane_width, bool is_last_lane>
using lane_profile_t = typename _lane_profile<dp_block_t, lane_width, is_last_lane>::type;

template <typename dp_block_t, size_t lane_width, bool is_last_lane>
class _lane_profile<dp_block_t, lane_width, is_last_lane>::type :
    public dp_matrix::lane_t<dp_block_t, lane_width, is_last_lane>
{
private:

    using base_t = dp_matrix::lane_t<dp_block_t, lane_width, is_last_lane>;
    using substitution_model_t = typename std::remove_cvref_t<dp_block_t>::substitution_model_type;
    using profile_t = typename substitution_model_t::template profile_type<lane_width>;

    profile_t _profile;

public:

    type() = delete;
    constexpr explicit type(dp_block_t dp_block, size_t const row_offset) noexcept :
        base_t{std::forward<dp_block_t>(dp_block), row_offset}
    {
        _profile = base_t::dp_block().substitution_model().initialise_profile(base_t::row_sequence(),
                                                                              strip_width<lane_width>);
    }
    ~type() = default;

    constexpr profile_t const & row_sequence() const noexcept
    {
        return _profile;
    }
};

namespace cpo {

struct _lane_profile_closure
{
    template <typename dp_block_t, size_t lane_width, bool is_last_lane>
    constexpr auto operator()(dp_block_t && dp_block,
                              size_t const row_offset,
                              lane_width_t<lane_width>,
                              std::bool_constant<is_last_lane>) const noexcept
    {
        using dp_lane_t = dp_matrix::lane_profile_t<dp_block_t, lane_width, is_last_lane>;

        return dp_lane_t{std::forward<dp_block_t>(dp_block), row_offset};
    }
};

} // namespace cpo
} // namespace dp_matrxix

inline constexpr dp_matrix::cpo::_lane_profile_closure dp_matrix_lane_profile{};
} // inline namespace v1
}  // namespace seqan::pairwise_aligner
