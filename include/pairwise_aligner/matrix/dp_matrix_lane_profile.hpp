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
#include <pairwise_aligner/score_model/strip_width.hpp>

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
    using last_lane_tag_t = typename base_t::last_lane_tag_type;
    using profile_t = typename substitution_model_t::template profile_type<last_lane_tag_t::width>;

    profile_t _profile;

public:

    _type() = delete;
    _type(wrappee_t wrappee) noexcept :
        base_t{std::move(wrappee)},
        _profile{base_t::substitution_model().initialise_profile(base_t::row_sequence(),
                                                                 strip_width<last_lane_tag_t::width>)}
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
    constexpr auto operator()(wrappee_fn_t && wrappee_fn) const noexcept
    {
        std::tuple<wrappee_fn_t> tmp{std::forward<wrappee_fn_t>(wrappee_fn)};
        return [fwd_capture = std::move(tmp)] (auto && ...args) {
            using fwd_wrappee_fn_t = std::tuple_element_t<0, decltype(fwd_capture)>;
            auto wrappee = std::invoke(std::forward<fwd_wrappee_fn_t &&>(get<0>(fwd_capture)),
                                       std::forward<decltype(args)>(args)...);

            return _type{std::move(wrappee)};
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
