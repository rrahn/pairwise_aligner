// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_matrix_column_local.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <pairwise_aligner/matrix/dp_matrix_data_handle.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace dp_matrix {
namespace detail {

template <typename base_substitution_model_t>
class local_substitution_model
{
    using score_t = typename base_substitution_model_t::score_type;

    base_substitution_model_t const & _model;
    score_t _zero{};

public:
    using score_type = score_t;

    local_substitution_model() = delete;
    explicit local_substitution_model(base_substitution_model_t const & model) : _model{model}
    {}

    template <typename mask_t>
    explicit local_substitution_model(base_substitution_model_t const & model, mask_t const & is_local_block) :
        _model{model},
        _zero{blend(is_local_block, model.zero(), score_t{std::numeric_limits<typename score_t::value_type>::lowest()})}
    {}

    template <typename ...args_t>
    constexpr auto score(args_t && ...args) const noexcept
    {
        using std::max;
        return max(_model.score(std::forward<args_t>(args)...), _zero);
    }

    template <typename ...args_t>
    constexpr decltype(auto) initialise_profile(args_t && ...args) noexcept
    {
        return _model.initialise_profile(std::forward<args_t>(args)...);
    }

    template <typename ...args_t>
    constexpr decltype(auto) initialise_profile(args_t && ...args) const noexcept
    {
        return _model.initialise_profile(std::forward<args_t>(args)...);
    }

    constexpr operator base_substitution_model_t & () noexcept
    {
        return reinterpret_cast<base_substitution_model_t &>(*this);
    }

    constexpr operator base_substitution_model_t const & () const noexcept
    {
        return reinterpret_cast<base_substitution_model_t const &>(*this);
    }
};

template <typename ...dp_data_t>
struct _dp_matrix_data_handle_local
{
    class type;
};

template <typename ...dp_data_t>
using dp_matrix_data_handle_local = typename _dp_matrix_data_handle_local<dp_data_t...>::type;

template <typename ...dp_data_t>
class _dp_matrix_data_handle_local<dp_data_t...>::type : public dp_matrix_data_handle<dp_data_t...>
{
    using base_t = dp_matrix_data_handle<dp_data_t...>;

public:
    using base_t::base_t;

    using substitution_model_type = local_substitution_model<typename base_t::substitution_model_type>;

    constexpr auto substitution_model() noexcept
    {
        return local_substitution_model{base_t::substitution_model()};
    }

    constexpr auto substitution_model() const noexcept
    {
        return local_substitution_model{base_t::substitution_model()};
    }

    template <typename mask_t>
    constexpr auto substitution_model(mask_t const & mask) noexcept
    {
        return local_substitution_model{base_t::substitution_model(), mask};
    }

    template <typename mask_t>
    constexpr auto substitution_model(mask_t const & mask) const noexcept
    {
        return local_substitution_model{base_t::substitution_model(), mask};
    }
};

} // namespace detail
} // namespace dp_matrix
} // inline namespace v1
}  // namespace seqan::pairwise_aligner
