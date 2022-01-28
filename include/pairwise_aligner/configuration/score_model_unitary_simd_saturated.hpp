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

#include <cmath>

#include <pairwise_aligner/configuration/score_model_unitary_simd.hpp>
#include <pairwise_aligner/dp_algorithm_template/dp_algorithm_template_saturated.hpp>
#include <pairwise_aligner/matrix/dp_vector_bulk.hpp>
#include <pairwise_aligner/matrix/dp_vector_chunk.hpp>
#include <pairwise_aligner/matrix/dp_vector_policy.hpp>
#include <pairwise_aligner/matrix/dp_vector_saturated.hpp>
#include <pairwise_aligner/matrix/dp_vector_single.hpp>
#include <pairwise_aligner/tracker/tracker_global_simd_saturated.hpp>

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
    using result_factory_type = tracker::global_simd_saturated::factory<original_score_type>;

    constexpr auto configure_substitution_policy() const noexcept
    {
        return score_model_type{static_cast<score_type>(_match_score), static_cast<score_type>(_mismatch_score)};
    }

    template <typename configuration_t>
    constexpr auto configure_result_factory_policy(configuration_t const & configuration) const noexcept
    {
        return result_factory_type{static_cast<original_score_type>(_match_score),
                                   configuration.trailing_gap_setting()};
    }

    template <typename common_configurations_t>
    constexpr auto configure_dp_vector_policy(common_configurations_t const & configuration) const
    {
        auto [zero_offset, max_block_size] = compute_max_block_size(configuration);

        using column_cell_t = typename common_configurations_t::dp_cell_column_type<score_type>;
        using original_column_cell_t = typename common_configurations_t::dp_cell_column_type<original_score_type>;

        auto column_vector =
            dp_vector_bulk_factory<score_type>(
                dp_vector_chunk_factory(
                    dp_vector_saturated_factory<original_column_cell_t>(dp_vector_single<column_cell_t>{}, zero_offset),
                    max_block_size
                )
            );

        using row_cell_t = typename common_configurations_t::dp_cell_row_type<score_type>;
        using original_row_cell_t = typename common_configurations_t::dp_cell_row_type<original_score_type>;

        auto row_vector =
            dp_vector_bulk_factory<score_type>(
                dp_vector_chunk_factory(
                    dp_vector_saturated_factory<original_row_cell_t>(dp_vector_single<row_cell_t>{}, zero_offset),
                    max_block_size
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

private:
    template <typename configuration_t>
    constexpr auto compute_max_block_size(configuration_t const & configuration) const noexcept
    {
        auto [block_size_gap, zero_offset_gap] = max_block_size_with_gaps(std::abs(configuration._gap_open_score),
                                                                          std::abs(configuration._gap_extension_score));
        auto [block_size_mismatch, zero_offset_mismatch] = max_block_size_with_mismatch();

        // Choose the max of both block sizes since due to the recursion the largest negative distance is affected
        // by the lowest negative score with the given block size.
        // Also set the corresponding zero offset accordingly.
        int8_t zero_offset{};
        size_t max_block_size = (block_size_gap < block_size_mismatch)
                              ? (zero_offset = zero_offset_mismatch, block_size_mismatch)
                              : (zero_offset = zero_offset_gap, block_size_gap);
        return std::pair{zero_offset, max_block_size};
    }

    constexpr auto max_block_size_with_gaps(float const gap_open, float const gap_extension) const noexcept
    {
        // Assume positive gap scores for the following formula.
        assert(gap_open >= 0);
        assert(gap_extension >= 0);

        // This formula computes the maximal block size and the respective zero offset for which the scores during
        // the computation of the alignment blocks in saturated mode are guaranteed to not exceed the value range
        // of the machines int8_t type, under the condition of having only gaps in a single block.
        // It is derived from solving the following equations:
        //  I: 127 >= x + match * block_size + gap_open + gap_extension * block_size
        // II: -128 <= x - 2 * (gap_open + gap_extension * block_size)
        // Here x represents the zero offset for which the block size is maximised.

        float const max_score = std::numeric_limits<int8_t>::max();
        float const min_score = std::numeric_limits<int8_t>::lowest();
        float const upper_score_limit = max_score - gap_open;
        float const block_divisor = _match_score + gap_extension;
        float const block_scale = 2 * gap_extension / block_divisor;

        int8_t const zero_offset = std::ceil((min_score + 2 * gap_open + (block_scale * upper_score_limit)) /
                                             (1 + block_scale));
        size_t const block_size = std::floor((upper_score_limit - zero_offset) / block_divisor);

        return std::pair{block_size, zero_offset};
    }

    constexpr auto max_block_size_with_mismatch() const noexcept
    {
        // Assume positive mismatch score for the following equation
        float const mismatch = std::abs(_mismatch_score);

        // This formula computes the maximal block size and the respective zero offset for which the scores during
        // the computation of the alignment blocks in saturated mode are guaranteed to not exceed the value range
        // of the machines int8_t type, under the condition of having only mismatches in a single block.
        // It is derived from solving the following equations:
        //  I: 127 >= x + match * block_size + mismatch * block_size
        // II: -128 <= x - mismatch * block_size
        // Here x represents the zero offset for which the block size is maximised.

        float const max_score = std::numeric_limits<int8_t>::max();
        float const min_score = std::numeric_limits<int8_t>::lowest();
        float const block_scale = mismatch / (_match_score + mismatch);

        int8_t const zero_offset = std::ceil((min_score + max_score * block_scale) / (1 + block_scale));
        size_t const block_size = std::floor((max_score - zero_offset) / (_match_score + mismatch));

        return std::pair{block_size, zero_offset};
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
        std::forward<next_configurator_t>(_next_configurator).set_config(std::forward<values_t>(values)..., _traits);
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
    auto apply(next_configurator_t && next_configurator) const
    {
        return _predecessor.apply(configurator_t<next_configurator_t, traits_t>{
                    std::forward<next_configurator_t>(next_configurator),
                    _traits
                });
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
    constexpr auto operator()(predecessor_t && predecessor,
                              score_t const match_score,
                              score_t const mismatch_score) const
    {
        using traits_t = traits<score_t>;
        return _score_model_unitary_simd_saturated::
            rule<predecessor_t, traits_t>{{},
                                          std::forward<predecessor_t>(predecessor),
                                          traits_t{match_score, mismatch_score}};
    }

    template <typename score_t>
    constexpr auto operator()(score_t const match_score, score_t const mismatch_score) const
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
