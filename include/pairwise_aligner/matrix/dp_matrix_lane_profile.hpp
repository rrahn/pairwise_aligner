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
#include <pairwise_aligner/matrix/dp_matrix_lane_width.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace dp_matrix {

namespace _lane_profile {

template <typename wrappee_t>
class _type : public wrappee_t
{
    using base_t = wrappee_t;
    using substitution_model_t = typename base_t::substitution_model_type;
    using profile_t = typename substitution_model_t::profile_type;

    profile_t _profile;

public:

    _type() = delete;
    _type(wrappee_t wrappee) noexcept :
        base_t{std::move(wrappee)},
        _profile{base_t::substitution_model().initialise_profile(base_t::row_sequence())}
    {}

    using base_t::base_t;

    constexpr profile_t const & row_sequence() const noexcept
    {
        return _profile;
    }
};

struct _fn
{
    template <typename wrappee_fn_t>
    constexpr auto operator()(wrappee_fn_t wrappee_fn) const noexcept
    {
        return [wrappee_fn = std::move(wrappee_fn)] <typename dp_state_t> (dp_state_t && dp_state, auto && ...args) {

            using substitution_model_t = typename dp_state_t::substitution_model_type;
            using profile_t = typename substitution_model_t::profile_type;

            profile_t profile{std::forward<dp_state_t>(dp_state)
                                    .substitution_model()
                                    .initialise_profile(std::forward<dp_state_t>(dp_state).row_sequence())};

            return std::invoke(
                    std::move(wrappee_fn),
                    dp_matrix::detail::make_dp_state(
                        std::forward<dp_state_t>(dp_state).dp_column(),
                        std::forward<dp_state_t>(dp_state).dp_row(),
                        std::forward<dp_state_t>(dp_state).column_sequence(),
                        std::move(profile),
                        std::forward<dp_state_t>(dp_state).substitution_model(),
                        std::forward<dp_state_t>(dp_state).tracker()
                    ),
                    std::forward<decltype(args)>(args)...);
        };
    }
};
} // namespace _lane_profile

inline namespace _cpo {
inline constexpr dp_matrix::_lane_profile::_fn lane_profile{};

} // inline namespace _cpo
} // namespace dp_matrix

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
