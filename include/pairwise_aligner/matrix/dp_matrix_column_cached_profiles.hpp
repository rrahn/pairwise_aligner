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

#include <pairwise_aligner/matrix/dp_matrix_block_cached_profile.hpp>
#include <pairwise_aligner/matrix/dp_matrix_column_local.hpp>
#include <pairwise_aligner/matrix/dp_matrix_column.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace dp_matrix {

namespace detail {

template <typename profile_t>
class cached_profile : public profile_t
{
    bool _is_initialsed{};

public:
    using profile_t::profile_t;
    constexpr cached_profile & operator=(profile_t profile) noexcept
    {
        static_cast<profile_t &>(*this) = std::move(profile);
        _is_initialsed = true;
        return *this;
    }

    constexpr operator bool() const noexcept
    {
        return _is_initialsed;
    }
};

template <typename substitution_model_t, size_t lane_width>
class substitution_model_profile_lane_cache
{
    using profile_t = typename substitution_model_t::profile_type<lane_width>;
    using cached_profile_t = cached_profile<profile_t>;

    substitution_model_t const & _substitution_model;
    std::vector<cached_profile_t> _cached_profiles;

public:

    using score_type = typename substitution_model_t::score_type;

    template <size_t _lane_width>
        requires (_lane_width == lane_width)
    using profile_type = std::add_lvalue_reference_t<cached_profile_t const>;

    substitution_model_profile_lane_cache() = delete;
    constexpr explicit substitution_model_profile_lane_cache(substitution_model_t const & substitution_model,
                                                             size_t const lanes_per_block) :
        _substitution_model{substitution_model}
    {
        _cached_profiles.resize(lanes_per_block);
    }

    template <typename ...args_t>
    constexpr auto score(args_t && ...args) const noexcept
    {
        return _substitution_model.score(std::forward<args_t>(args)...);
    }

    template <typename ...args_t>
    constexpr cached_profile_t const & initialise_profile(size_t const index, args_t && ...args) noexcept
    {
        assert(index < _cached_profiles.size());
        if (!_cached_profiles[index])
            _cached_profiles[index] = _substitution_model.initialise_profile(std::forward<args_t>(args)...);

        assert(_cached_profiles[index]); // profile must be initialised
        return _cached_profiles[index];
    }
};

} // namespace detail

template <bool is_local, typename block_closure_t, typename ...dp_data_t>
struct _column_cached_profiles
{
    class type;
};

template <typename block_closure_t, typename ...dp_data_t>
using column_cached_profiles_t = typename _column_cached_profiles<false, block_closure_t, dp_data_t...>::type;

template <typename block_closure_t, typename ...dp_data_t>
using column_cached_profiles_local_t = typename _column_cached_profiles<true, block_closure_t, dp_data_t...>::type;

template <bool is_local, typename block_closure_t, typename ...dp_data_t>
class _column_cached_profiles<is_local, block_closure_t, dp_data_t...>::type :
    public std::conditional_t<is_local,
                              dp_matrix::column_local_t<block_closure_t, dp_data_t...>,
                              dp_matrix::column_t<block_closure_t, dp_data_t...>>
{
protected:
    using base_t = std::conditional_t<is_local,
                              dp_matrix::column_local_t<block_closure_t, dp_data_t...>,
                              dp_matrix::column_t<block_closure_t, dp_data_t...>>;

    using typename base_t::substitution_model_type;
    using cached_profiles_t =
        detail::substitution_model_profile_lane_cache<substitution_model_type, base_t::lane_width_v>;

    cached_profiles_t _cached_profiles;

public:
    type() = delete;
    constexpr explicit type(block_closure_t block_closure, dp_data_t ...dp_data) noexcept :
        base_t{std::move(block_closure), std::forward<dp_data_t>(dp_data)...},
        _cached_profiles{base_t::substitution_model(), base_t::lanes_per_block()}
    {}

    ~type() noexcept
    {}

    // So different implementations can have different columns
    constexpr auto operator[](size_t const index) noexcept
    {
        assert(index < base_t::size());
        return base_t::make_matrix_block(base_t::column()[index],
                                         base_t::row(),
                                         substitution_model(),
                                         base_t::tracker(),
                                         base_t::row_sequence(),
                                         base_t::lane_width());
    }

    // overwrites substitution model!?
    constexpr cached_profiles_t & substitution_model() noexcept
    {
        return _cached_profiles;
    }

    constexpr cached_profiles_t const & substitution_model() const noexcept
    {
        return _cached_profiles;
    }
};

// namespace cpo {
// template <typename block_closure_t = dp_matrix::cpo::_block_cached_profile_closure<>>
// struct _column_cached_profiles_closure
// {
//     block_closure_t block_closure{};

//     // One of the other columns need to be set as the derived class from this.
//     // Then they invoke
//     template <typename ...dp_data_t>
//     constexpr auto operator()(dp_data_t && ...dp_data) const noexcept {
//         using dp_column_t = dp_matrix::column_cached_profiles_t<block_closure_t, dp_data_t...>;

//         return dp_column_t{block_closure, std::forward<dp_data_t>(dp_data)...};
//     }
// };

// } // namespace cpo
} // namespace dp_matrix

// inline constexpr dp_matrix::cpo::_column_cached_profiles_closure<> dp_matrix_column_cached_profiles{};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
