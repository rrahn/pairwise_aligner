// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::cfg::gap_model_affine.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <seqan3/std/type_traits>

#include <pairwise_aligner/configuration/initial.hpp>
#include <pairwise_aligner/configuration/rule_gap_model.hpp>
#include <pairwise_aligner/affine/affine_cell.hpp>
#include <pairwise_aligner/affine/affine_dp_algorithm.hpp>
#include <pairwise_aligner/affine/affine_gap_model.hpp>
#include <pairwise_aligner/affine/affine_initialisation_strategy.hpp>
#include <pairwise_aligner/type_traits.hpp>
#include <pairwise_aligner/utility/type_list.hpp>

namespace seqan::pairwise_aligner {
inline namespace v1
{
namespace cfg
{
namespace _gap_model_affine
{

// ----------------------------------------------------------------------------
// traits
// ----------------------------------------------------------------------------

struct traits
{
    template <typename score_t>
    using gap_model_type = affine_gap_model<score_t>;

    template <typename gap_model_t>
    using dp_initialisation_type = affine_initialisation_strategy<gap_model_t>;

    // Offer the score type here.
    template <typename score_t>
    using dp_cell_column_type = affine_cell<score_t, dp_vector_order::column>;

    template <typename score_t>
    using dp_cell_row_type = affine_cell<score_t, dp_vector_order::row>;

    // Offer some overload for the column type.
    template <template <typename > typename dp_template_t,
              typename dp_score_model_t,
              typename dp_gap_model_t,
              typename dp_initialisation_t>
    using dp_kernel_type = affine_dp_algorithm<dp_template_t,
                                               dp_score_model_t,
                                               dp_gap_model_t,
                                               dp_initialisation_t>;
};

// ----------------------------------------------------------------------------
// configurator
// ----------------------------------------------------------------------------

template <typename next_configurator_t, typename gap_model_t>
struct _configurator
{
    struct type;
};

template <typename next_configurator_t, typename gap_model_t>
using configurator_t = typename _configurator<next_configurator_t, gap_model_t>::type;

template <typename next_configurator_t, typename gap_model_t>
struct _configurator<next_configurator_t, gap_model_t>::type
{
    next_configurator_t _next_configurator;
    gap_model_t _gap_model;

    template <typename ...values_t>
    void set_config(values_t && ... values) && noexcept
    {
        std::forward<next_configurator_t>(_next_configurator).set_config(std::forward<values_t>(values)...,
                                                                         std::forward<gap_model_t>(_gap_model));
    }
};

// ----------------------------------------------------------------------------
// rule
// ----------------------------------------------------------------------------

template <typename predecessor_t, typename gap_model_t>
struct _rule
{
    struct type;
};

template <typename predecessor_t, typename gap_model_t>
using rule = typename _rule<predecessor_t, gap_model_t>::type;

template <typename predecessor_t, typename gap_model_t>
struct _rule<predecessor_t, gap_model_t>::type : cfg::gap_model::rule<predecessor_t>
{
    predecessor_t _predecessor;
    gap_model_t _gap_model;

    using traits_type = type_list<traits>;

    template <template <typename ...> typename type_list_t>
    using configurator_types = typename concat_type_lists_t<configurator_types_t<std::remove_cvref_t<predecessor_t>,
                                                                                 type_list>,
                                                            traits_type>::template apply<type_list_t>;

    template <typename next_configurator_t>
    auto apply(next_configurator_t && next_configurator)
    {
        return std::forward<predecessor_t>(_predecessor).apply(
                configurator_t<std::remove_cvref_t<next_configurator_t>, gap_model_t>{
                                std::forward<next_configurator_t>(next_configurator),
                                std::forward<gap_model_t>(_gap_model)});
    }
};

// ----------------------------------------------------------------------------
// CPO
// ----------------------------------------------------------------------------

namespace _cpo
{
struct _fn
{
    // implementation of function style connection
    template <typename predecessor_t, typename score_t>
    auto operator()(predecessor_t && predecessor,
                    score_t const gap_open_score,
                    score_t const gap_extension_score) const
    {
        using gap_model_t = pairwise_aligner::affine_gap_model<score_t>;
        return _gap_model_affine::rule<predecessor_t, gap_model_t>{{},
                                                                   std::forward<predecessor_t>(predecessor),
                                                                   gap_model_t{gap_open_score, gap_extension_score}};
    }

    template <typename score_t>
    auto operator()(score_t const gap_open_score, score_t const gap_extension_score) const
    {
        return this->operator()(cfg::initial, gap_open_score, gap_extension_score);
    }
};
} // namespace _cpo
} // namespace _gap_model_affine

inline constexpr _gap_model_affine::_cpo::_fn gap_model_affine{};

} // namespace cfg
} // inline namespace v1
}  // namespace seqan::pairwise_aligner
