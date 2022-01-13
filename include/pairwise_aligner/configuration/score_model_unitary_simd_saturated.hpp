// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::cfg::score_model_unitary_simd_saturated.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <pairwise_aligner/configuration/score_model_unitary_simd.hpp>
#include <pairwise_aligner/dp_algorithm_template/dp_algorithm_template_saturated.hpp>
#include <pairwise_aligner/matrix/dp_vector_bulk.hpp>
#include <pairwise_aligner/matrix/dp_vector_chunk.hpp>
#include <pairwise_aligner/matrix/dp_vector_policy.hpp>
#include <pairwise_aligner/matrix/dp_vector_saturated.hpp>
#include <pairwise_aligner/matrix/dp_vector_single.hpp>
#include <pairwise_aligner/result/result_factory_chunk.hpp>

namespace seqan::pairwise_aligner {
inline namespace v1
{
namespace cfg
{
namespace _score_model_unitary_simd_saturated
{

// ----------------------------------------------------------------------------
// traits
// ----------------------------------------------------------------------------

template <typename score_t>
struct traits
{
    static constexpr cfg::detail::rule_category category = cfg::detail::rule_category::score_model;

    score_t _match_score;
    score_t _mismatch_score;

    // Offer the score type here.
    using score_type = simd_score<int8_t>;
    using original_score_type = simd_score<score_t, score_type::size>;

    // substitution policy configurator
    using score_model_type = seqan::pairwise_aligner::score_model_unitary<score_type>;

    // result_factory configurator
    using result_factory_type = result_factory_bulk<original_score_type>;

    constexpr auto configure_substitution_policy() const noexcept
    {
        return score_model_type{static_cast<score_type>(_match_score), static_cast<score_type>(_mismatch_score)};
    }

    constexpr auto configure_result_factory_policy() const noexcept
    {
        return result_factory_chunk<result_factory_type>{result_factory_type{static_cast<original_score_type>(_match_score)}};
    }

    template <typename common_configurations_t>
    constexpr auto configure_dp_vector_policy(common_configurations_t const & configuration) const
    {
        size_t max_score = std::numeric_limits<int8_t>::max();
        size_t block_size_mismatch = max_score / (_match_score - _mismatch_score);
        size_t block_size_gap = (max_score + configuration._gap_open_score) / (_match_score - configuration._gap_extension_score);
        size_t block_size = std::max(block_size_mismatch, block_size_gap);

        using column_cell_t = typename common_configurations_t::dp_cell_column_type<score_type>;
        using original_column_cell_t = typename common_configurations_t::dp_cell_column_type<original_score_type>;

        auto column_vector =
            dp_vector_bulk_factory<score_type>(
                dp_vector_chunk_factory(
                    dp_vector_saturated_factory<original_column_cell_t>(dp_vector_single<column_cell_t>{}),
                    block_size
                )
            );

        using row_cell_t = typename common_configurations_t::dp_cell_row_type<score_type>;
        using original_row_cell_t = typename common_configurations_t::dp_cell_row_type<original_score_type>;

        auto row_vector =
            dp_vector_bulk_factory<score_type>(
                dp_vector_chunk_factory(
                    dp_vector_saturated_factory<original_row_cell_t>(dp_vector_single<row_cell_t>{}),
                    block_size
                )
            );

        return dp_vector_policy{std::move(column_vector), std::move(row_vector)};
    }

    template <typename configuration_t, typename ...policies_t>
    constexpr auto configure_algorithm(configuration_t const &, policies_t && ...policies) const noexcept
    {
        using algorithm_t = typename configuration_t::algorithm_type<dp_algorithm_template_saturated,
                                                                     std::remove_cvref_t<policies_t>...>;

        return interface_one_to_one_bulk<algorithm_t, score_type::size>{algorithm_t{std::move(policies)...}};
    }
};

// ----------------------------------------------------------------------------
// configurator
// ----------------------------------------------------------------------------

template <typename next_configurator_t, typename traits_t>
struct _configurator
{
    struct type;
};

template <typename next_configurator_t, typename traits_t>
using configurator_t = typename _configurator<next_configurator_t, traits_t>::type;

template <typename next_configurator_t, typename traits_t>
struct _configurator<next_configurator_t, traits_t>::type
{
    next_configurator_t _next_configurator;
    traits_t _traits;

    template <typename ...values_t>
    void set_config(values_t && ... values) noexcept
    {
        std::forward<next_configurator_t>(_next_configurator).set_config(std::forward<values_t>(values)...,
                                                                         std::forward<traits_t>(_traits));
    }
};

// ----------------------------------------------------------------------------
// rule
// ----------------------------------------------------------------------------

template <typename predecessor_t, typename traits_t>
struct _rule
{
    struct type;
};

template <typename predecessor_t, typename traits_t>
using rule = typename _rule<predecessor_t, traits_t>::type;

template <typename predecessor_t, typename traits_t>
struct _rule<predecessor_t, traits_t>::type : cfg::score_model::rule<predecessor_t>
{
    predecessor_t _predecessor;
    traits_t _traits;

    using traits_type = type_list<traits_t>;

    template <template <typename ...> typename type_list_t>
    using configurator_types = typename concat_type_lists_t<configurator_types_t<std::remove_cvref_t<predecessor_t>,
                                                                                 type_list>,
                                                            traits_type>::template apply<type_list_t>;

    template <typename next_configurator_t>
    auto apply(next_configurator_t && next_configurator)
    {
        return std::forward<predecessor_t>(_predecessor).apply(
                configurator_t<next_configurator_t, traits_t>{
                                std::forward<next_configurator_t>(next_configurator),
                                std::forward<traits_t>(_traits)});
    }
};

// ----------------------------------------------------------------------------
// CPO
// ----------------------------------------------------------------------------

namespace _cpo
{
struct _fn
{
    template <typename predecessor_t, typename score_t>
    auto operator()(predecessor_t && predecessor, score_t const match_score, score_t const mismatch_score) const
    {
        using traits_t = traits<score_t>;
        return _score_model_unitary_simd_saturated::rule<predecessor_t, traits_t>{{},
                                                                                  std::forward<predecessor_t>(predecessor),
                                                                                  traits_t{match_score, mismatch_score}};
    }

    template <typename score_t>
    auto operator()(score_t const match_score, score_t const mismatch_score) const
    {
        return this->operator()(cfg::initial, match_score, mismatch_score);
    }
};
} // namespace _cpo
} // namespace _score_model_unitary_simd_saturated

inline constexpr _score_model_unitary_simd_saturated::_cpo::_fn score_model_unitary_simd_saturated{};

} // namespace cfg
} // inline namespace v1
}  // namespace seqan::pairwise_aligner
