// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_matrix_policies.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <concepts>

#include <pairwise_aligner/matrix/dp_matrix.hpp>
#include <pairwise_aligner/utility/type_list.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace dp_matrix {

namespace _dp_matrix_local {
namespace detail {

template <std::move_constructible substitution_model_t>
class local_substitution_model : public substitution_model_t
{
public:

    using substitution_model_t::substitution_model_t;
    local_substitution_model() = default;
    constexpr local_substitution_model(substitution_model_t model) : substitution_model_t{std::move(model)}
    {}

    template <typename ...args_t>
    constexpr auto score(args_t && ...args) const noexcept
        -> decltype(substitution_model_t::score(std::forward<args_t>(args)...))
    {
        using score_t = decltype(substitution_model_t::score(std::forward<args_t>(args)...));

        using std::max;
        return max(substitution_model_t::score(std::forward<args_t>(args)...), static_cast<score_t>(0));
    }
};
} // namespace detail

// adaptor closure object
struct _fn
{
    template <typename dp_column_fn_t>
    constexpr auto operator()(dp_column_fn_t dp_column_fn) const noexcept
    {
        return [dp_column_fn = std::move(dp_column_fn)] <typename dp_state_t>(dp_state_t && dp_state) {

            using original_substitution_model_t = typename dp_state_t::substitution_model_type;
            using local_substitution_model_t = detail::local_substitution_model<original_substitution_model_t>;

            auto modified_state = dp_matrix::detail::make_dp_state(
                std::forward<dp_state_t>(dp_state).dp_column(),
                std::forward<dp_state_t>(dp_state).dp_row(),
                std::forward<dp_state_t>(dp_state).column_sequence(),
                std::forward<dp_state_t>(dp_state).row_sequence(),
                local_substitution_model_t{std::forward<dp_state_t>(dp_state).substitution_model()},
                std::forward<dp_state_t>(dp_state).tracker()
            );
            using dp_matrix_t = _dp_matrix::_type<dp_column_fn_t, decltype(modified_state)>;
            return dp_matrix_t{std::move(dp_column_fn), std::move(modified_state)};
        };
    }
};

} // namespace _dp_matrix_local
inline namespace _cpo {

inline constexpr _dp_matrix_local::_fn matrix_local{};

} // inline namespace _cpo
} // namespace dp_matrix
} // inline namespace v1
}  // namespace seqan::pairwise_aligner
