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

#include <seqan3/std/concepts>

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
    constexpr auto operator()(dp_column_fn_t && dp_column_fn) const noexcept
    {
        std::tuple<dp_column_fn_t> tmp{std::forward<dp_column_fn_t>(dp_column_fn)};
        return [fwd_capture = std::move(tmp)] (auto && ...dp_state) {
            constexpr size_t idx =
                dp_matrix::detail::dp_state_accessor_id_v<dp_matrix::detail::dp_state_accessor::id_substitution_model>;

            using substitution_model_t = std::remove_reference_t<seqan3::pack_traits::at<idx, decltype(dp_state)...>>;
            using modified_pack_list_t =
                    seqan3::list_traits::transform<remove_rvalue_reference_t,
                        seqan3::pack_traits::replace_at<detail::local_substitution_model<substitution_model_t>,
                                                        idx,
                                                        decltype(dp_state)...
                        >
                    >;

            using fwd_dp_column_fn_t = std::tuple_element_t<0, decltype(fwd_capture)>;
            using dp_matrix_t = apply_t<_dp_matrix::_type,
                                        concat_type_lists_t<type_list<fwd_dp_column_fn_t>, modified_pack_list_t>>;
            return dp_matrix_t{std::forward<fwd_dp_column_fn_t>(get<0>(fwd_capture)),
                               std::forward<decltype(dp_state)>(dp_state)...};
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
