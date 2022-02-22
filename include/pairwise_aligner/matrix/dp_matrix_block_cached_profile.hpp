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

#include <pairwise_aligner/matrix/dp_matrix_block.hpp>
#include <pairwise_aligner/matrix/dp_matrix_lane_profile.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace dp_matrix {

template <typename lane_closure_t, typename ...dp_data_t>
struct _block_cached_profile
{
    class type;
};

template <typename lane_closure_t, typename ...dp_data_t>
using block_cached_profile_t = typename _block_cached_profile<lane_closure_t, dp_data_t...>::type;

template <typename lane_closure_t, typename ...dp_data_t>
class _block_cached_profile<lane_closure_t, dp_data_t...>::type :
    public dp_matrix::block_t<lane_closure_t, dp_data_t...>
{
    using base_t = dp_matrix::block_t<lane_closure_t, dp_data_t...>;

    struct proxy_reference
    {
        using row_type = typename base_t::row_type;
        using column_type = typename base_t::column_type;
        using substitution_model_type = typename base_t::substitution_model_type;

        struct substitution_model_wrapper
        {
            substitution_model_type & _substitution_model;
            size_t index{};

            template <typename ...args_t>
            constexpr auto const & initialise_profile(args_t && ...args) const noexcept
            {
                return _substitution_model.initialise_profile(index, std::forward<args_t>(args)...);
            }

            template <typename ...args_t>
            constexpr auto & initialise_profile(args_t && ...args) noexcept
            {
                return _substitution_model.initialise_profile(index, std::forward<args_t>(args)...);
            }
        };

        type & _parent;
        size_t const _lane_index;

        constexpr typename base_t::column_type & column() noexcept
        {
            return _parent.column();
        }

        constexpr typename base_t::column_type const & column() const noexcept
        {
            return _parent.column();
        }

        constexpr typename base_t::row_type & row() noexcept
        {
            return _parent.row();
        }

        constexpr typename base_t::row_type const & row() const noexcept
        {
            return _parent.row();
        }

        constexpr auto row_sequence() const noexcept
        {
            return _parent.row_sequence();
        }

        constexpr auto substitution_model() noexcept
        {
            return substitution_model_wrapper{_parent.substitution_model(), _lane_index};
        }

        constexpr auto substitution_model() const noexcept
        {
            return substitution_model_wrapper{_parent.substitution_model(), _lane_index};
        }
    };

public:

    using typename base_t::row_type;
    using typename base_t::column_type;

    type() = delete;
    constexpr explicit type(lane_closure_t lane_closure, dp_data_t ...dp_data) noexcept :
        base_t{std::move(lane_closure), std::forward<dp_data_t>(dp_data)...}
    {}

    constexpr auto operator[](size_t const index) noexcept
    {
        return base_t::make_block_lane(proxy_reference{*this, index},
                                       index * base_t::lane_width_v,
                                       base_t::lane_width(),
                                       std::false_type{});
    }

    constexpr auto last_lane() noexcept
    {
        return base_t::make_block_lane(proxy_reference{*this, base_t::lanes_per_block() - 1},
                                       ((base_t::row().size() - 1) / base_t::lane_width_v) * base_t::lane_width_v,
                                       base_t::lane_width(),
                                       std::true_type{});
    }
};

namespace cpo {

template <typename lane_closure_t = dp_matrix::cpo::_lane_profile_closure>
struct _block_cached_profile_closure
{
    lane_closure_t lane_closure{};

    template <typename ...dp_data_t>
    constexpr auto operator()(dp_data_t && ...dp_data) const noexcept {
        using dp_block_t = dp_matrix::block_cached_profile_t<lane_closure_t, dp_data_t...>;

        return dp_block_t{lane_closure, std::forward<dp_data_t>(dp_data)...};
    }
};

} // namespace cpo
} // namespace dp_matrxix

inline constexpr dp_matrix::cpo::_block_cached_profile_closure<> dp_matrix_block_cached_profile{};
} // inline namespace v1
} // namespace seqan::pairwise_aligner
