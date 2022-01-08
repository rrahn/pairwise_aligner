// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::cfg::configure_aligner_saturated.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <iostream>
#include <seqan3/std/concepts>
#include <seqan3/std/type_traits>

#include <pairwise_aligner/affine/affine_initialisation_strategy.hpp>
#include <pairwise_aligner/dp_algorithm_template/dp_algorithm_template_saturated.hpp>
#include <pairwise_aligner/matrix/dp_vector_grouped.hpp>
#include <pairwise_aligner/matrix/dp_vector_saturated.hpp>
#include <pairwise_aligner/matrix/dp_vector_single.hpp>
#include <pairwise_aligner/result/result_factory_chunk.hpp>
#include <pairwise_aligner/dp_trailing_gaps.hpp>
#include <pairwise_aligner/simd/simd_score_type.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace cfg
{
namespace _configure_aligner_saturated
{

// ----------------------------------------------------------------------------
// traits
// ----------------------------------------------------------------------------

template <typename configurator_types, size_t score_model_position, size_t gap_model_position>
struct traits
{
    // Extract the respective configuration traits.
    using score_model_config_traits_t = std::tuple_element_t<score_model_position, configurator_types>;
    using gap_model_config_traits_t = std::tuple_element_t<gap_model_position, configurator_types>;

    // Define all algorithm entities.
    // here we don't want to have a max score type.
    static constexpr size_t max_size = simd_score<int8_t>::size;
    using saturated_score_t = simd_score<int8_t, max_size>;

    // part of configurator
    // banded (bulk or single) + regular (bulk or single) + saturated (implicit bulk) -> cannot work with non-bulk interface?
    // we keep the original and need to put some type inside!
    using original_score_t = typename score_model_config_traits_t::template score_type<max_size>;
    // using original_score_t = saturated_score_t;

    using dp_cell_column_t = typename gap_model_config_traits_t::template dp_cell_column_type<saturated_score_t>;
    using dp_cell_column_original_t = typename gap_model_config_traits_t::template dp_cell_column_type<original_score_t>;
    using dp_column_saturated_t = dp_vector_saturated<dp_vector_single<dp_cell_column_t>, dp_cell_column_original_t>;
    using dp_column_chunked_t = dp_vector_grouped<dp_column_saturated_t>;
    using dp_vector_column_t = typename score_model_config_traits_t::template dp_vector_column_type<dp_column_chunked_t,
                                                                                                    saturated_score_t>;

    using dp_cell_row_t = typename gap_model_config_traits_t::template dp_cell_row_type<saturated_score_t>;
    using dp_cell_row_original_t = typename gap_model_config_traits_t::template dp_cell_row_type<original_score_t>;
    using dp_row_saturated_t = dp_vector_saturated<dp_vector_single<dp_cell_row_t>, dp_cell_row_original_t>;
    using dp_row_chunked_t = dp_vector_grouped<dp_row_saturated_t>;
    using dp_vector_row_t = typename score_model_config_traits_t::template dp_vector_row_type<dp_row_chunked_t,
                                                                                              saturated_score_t>;

    // Get the instantiated model types.
    using score_model_t = typename score_model_config_traits_t::score_model_type<saturated_score_t>;
    using result_factory_t = result_factory_chunk<typename score_model_config_traits_t::result_factory_type<original_score_t>>;

    // we need to set some parameters to the template without removing its types
    using dp_kernel_t = typename gap_model_config_traits_t:: template dp_kernel_type<dp_algorithm_template_saturated,
                                                                                     score_model_t,
                                                                                     result_factory_t>;
    // define the pairwise aligner type.
    using aligner_type = typename score_model_config_traits_t::template dp_interface_type<dp_kernel_t,
                                                                                          dp_vector_column_t,
                                                                                          dp_vector_row_t,
                                                                                          max_size>;
};

// ----------------------------------------------------------------------------
// configurator
// ----------------------------------------------------------------------------

template <typename pairwise_aligner_ref_t,
          typename score_t,
          typename bulk_score_t,
          int32_t score_model_index,
          int32_t gap_model_index,
          int32_t method_index>
struct _configurator
{
    struct type;
};

template <typename pairwise_aligner_ref_t,
          typename score_t,
          typename bulk_score_t,
          int32_t score_model_index,
          int32_t gap_model_index,
          int32_t method_index>
using configurator_t = typename _configurator<pairwise_aligner_ref_t,
                                              score_t,
                                              bulk_score_t,
                                              score_model_index,
                                              gap_model_index,
                                              method_index>::type;

template <typename pairwise_aligner_ref_t,
          typename score_t,
          typename bulk_score_t,
          int32_t score_model_index,
          int32_t gap_model_index,
          int32_t method_index>
struct _configurator<pairwise_aligner_ref_t, score_t, bulk_score_t, score_model_index, gap_model_index, method_index>::type
{
    pairwise_aligner_ref_t _aligner_ref;
    using pairwise_aligner_t = typename pairwise_aligner_ref_t::type;

    type() = delete;
    explicit type(pairwise_aligner_ref_t aligner_ref) : _aligner_ref{aligner_ref}
    {}

    template <typename ...actions_t>
    void set_config(actions_t && ...actions) && noexcept
    {
        std::tuple tpl_values{std::forward<actions_t>(actions)...};

        initialisation_rule init_rule{};
        trailing_gap_setting trailing_rule{};

        if constexpr (method_index != -1)
        {
            std::tie(init_rule, trailing_rule) = get<method_index>(tpl_values).create();
        }

        auto [score_model, result_factory] = get<score_model_index>(tpl_values).template create<bulk_score_t, score_t>();
        auto [gap_params, initialisation_rule] = get<gap_model_index>(tpl_values).create(init_rule);

        using result_factory_t = result_factory_chunk<std::remove_cvref_t<decltype(result_factory)>>;
        _aligner_ref.get() = pairwise_aligner_t{std::move(score_model),
                                                result_factory_t{std::move(result_factory)},
                                                std::move(gap_params),
                                                std::move(initialisation_rule),
                                                std::move(trailing_rule)};
    }
};

// ----------------------------------------------------------------------------
// implementation
// ----------------------------------------------------------------------------

template <typename predecessor_t>
inline constexpr auto _impl(predecessor_t && predecessor)
{
    // Define the value types.
    using pure_predecessor_t = std::remove_cvref_t<predecessor_t>;
    using configurator_types = typename pure_predecessor_t::template configurator_types<std::tuple>;

    // Get configuration traits positions.
    constexpr size_t score_model_position = pure_predecessor_t::configurator_index_of[0];
    constexpr size_t gap_model_position = pure_predecessor_t::configurator_index_of[1];
    constexpr size_t method_position = pure_predecessor_t::configurator_index_of[2];

    static_assert(score_model_position != -1, "The score model category is required!");
    static_assert(gap_model_position != -1, "The gap model category is required!");

    using traits_t = traits<configurator_types, score_model_position, gap_model_position>;
    using original_score_t = typename traits_t::original_score_t;
    using saturated_score_t = typename traits_t::saturated_score_t;
    using aligner_t = typename traits_t::aligner_type;

    // TODO: convert positions to list of indices.

    aligner_t aligner{};
    configurator_t<std::reference_wrapper<aligner_t>,
                   original_score_t,
                   saturated_score_t,
                   score_model_position,
                   gap_model_position,
                   method_position> cfg{std::ref(aligner)};

    // Store state for the operation on the stack.
    auto operation = std::forward<predecessor_t>(predecessor).apply(std::move(cfg));

    operation.start();

    return aligner;
}

// ----------------------------------------------------------------------------
// CPO
// ----------------------------------------------------------------------------

namespace _cpo
{

struct _fn
{
    template <typename predecessor_t>
    auto operator()(predecessor_t && predecessor) const noexcept
    {
        return _configure_aligner_saturated::_impl(std::forward<predecessor_t>(predecessor));
    }
};

} // namespace _cpo
} // namespace _configure_aligner_saturated

inline constexpr _configure_aligner_saturated::_cpo::_fn configure_aligner_saturated{};

} // namespace cfg
} // inline namespace v1
}  // namespace seqan::pairwise_aligner
