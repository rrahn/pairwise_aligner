// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::cfg::configure_aligner.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <iostream>
#include <seqan3/std/concepts>
#include <seqan3/std/type_traits>

#include <pairwise_aligner/affine/affine_initialisation_strategy.hpp>
#include <pairwise_aligner/dp_algorithm_template/dp_algorithm_template_standard.hpp>


namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace cfg
{
namespace _configure_aligner
{

// ----------------------------------------------------------------------------
// configurator
// ----------------------------------------------------------------------------

template <typename pairwise_aligner_ref_t>
struct _configurator
{
    struct type;
};

template <typename pairwise_aligner_ref_t>
using configurator_t = typename _configurator<pairwise_aligner_ref_t>::type;

template <typename pairwise_aligner_ref_t>
struct _configurator<pairwise_aligner_ref_t>::type
{
    pairwise_aligner_ref_t _aligner_ref;
    using pairwise_aligner_t = typename pairwise_aligner_ref_t::type;

    type() = delete;
    explicit type(pairwise_aligner_ref_t aligner_ref) : _aligner_ref{aligner_ref}
    {}

    template <typename ...values_t>
    void set_config(values_t && ...values) && noexcept
    {
        std::tuple tpl_values{std::forward<values_t>(values)...};

        _aligner_ref.get() = pairwise_aligner_t{get<0>(tpl_values),
                                                get<1>(tpl_values),
                                                affine_initialisation_strategy{get<1>(tpl_values)}};
    }
};

// ----------------------------------------------------------------------------
// implementation
// ----------------------------------------------------------------------------

template <typename predecessor_t>
inline constexpr auto _impl(predecessor_t && predecessor)
{
    // Define the value types.
    using configurator_types = typename predecessor_t::template configurator_types<std::tuple>;

    // Get configuration traits positions.
    constexpr size_t score_model_position = std::remove_cvref_t<predecessor_t>::configurator_index_of[0];
    constexpr size_t gap_model_position = std::remove_cvref_t<predecessor_t>::configurator_index_of[1];

    // Extract the respective configuration traits.
    using score_model_config_traits_t = std::tuple_element_t<score_model_position, configurator_types>;
    using gap_model_config_traits_t = std::tuple_element_t<gap_model_position, configurator_types>;

    // Define all algorithm entities.
    using score_t = typename score_model_config_traits_t::score_type;
    using dp_cell_row_t = typename gap_model_config_traits_t::template dp_cell_row_type<score_t>;
    using dp_cell_column_t = typename gap_model_config_traits_t::template dp_cell_column_type<score_t>;
    using dp_vector_row_t = typename score_model_config_traits_t::template dp_vector_row_type<dp_cell_row_t>;
    using dp_vector_column_t = typename score_model_config_traits_t::template dp_vector_column_type<dp_cell_column_t>;

    // Get the instantiated model types.
    using score_model_t = typename score_model_config_traits_t::score_model_type;
    using gap_model_t = typename gap_model_config_traits_t::template gap_model_type<score_t>;
    using initialisation_t = typename gap_model_config_traits_t::template dp_initialisation_type<gap_model_t>;

    // Define the kernel type.
    using dp_kernel_t = typename gap_model_config_traits_t:: template dp_kernel_type<dp_algorithm_template_standard,
                                                                                     score_model_t,
                                                                                     gap_model_t,
                                                                                     initialisation_t>;
    // define the pairwise aligner type.
    using pairwise_aligner_t = typename score_model_config_traits_t::template dp_interface_type<dp_kernel_t,
                                                                                                dp_vector_column_t,
                                                                                                dp_vector_row_t>;

    static_assert(std::same_as<score_t, int32_t>, "Not the right score type.");

    pairwise_aligner_t aligner{};

    _configure_aligner::configurator_t<std::reference_wrapper<pairwise_aligner_t>> cfg{std::ref(aligner)};

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
        return _configure_aligner::_impl(std::forward<predecessor_t>(predecessor));
    }
};

} // namespace _cpo
} // namespace _configure_aligner

inline constexpr _configure_aligner::_cpo::_fn configure_aligner{};

} // namespace cfg
} // inline namespace v1
}  // namespace seqan::pairwise_aligner
