// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::cfg::score_model_matrix_simd_saturated_NxN.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <ranges>
#include <type_traits>
#include <utility>

#include <pairwise_aligner/configuration/initial.hpp>
#include <pairwise_aligner/configuration/rule_score_model.hpp>
#include <pairwise_aligner/configuration/saturated_block_handler.hpp>
#include <pairwise_aligner/alphabet_conversion/alphabet_rank_map_simd.hpp>
#include <pairwise_aligner/dp_algorithm_template/dp_algorithm_template_standard.hpp>
#include <pairwise_aligner/interface/interface_one_to_one_bulk.hpp>
// #include <pairwise_aligner/matrix/dp_matrix_block_cached_profile.hpp>
#include <pairwise_aligner/matrix/dp_matrix_block.hpp>
#include <pairwise_aligner/matrix/dp_matrix_column_saturated_local.hpp>
#include <pairwise_aligner/matrix/dp_matrix_column_saturated.hpp>
#include <pairwise_aligner/matrix/dp_matrix_lane.hpp>
// #include <pairwise_aligner/matrix/dp_matrix_lane_profile.hpp>
// #include <pairwise_aligner/matrix/dp_matrix_lane_width.hpp>
#include <pairwise_aligner/matrix/dp_matrix.hpp>
#include <pairwise_aligner/matrix/dp_vector_bulk.hpp>
#include <pairwise_aligner/matrix/dp_vector_chunk.hpp>
#include <pairwise_aligner/matrix/dp_vector_policy.hpp>
#include <pairwise_aligner/matrix/dp_vector_offset_transformation.hpp>
#include <pairwise_aligner/matrix/dp_vector_rank_transformation.hpp>
#include <pairwise_aligner/matrix/dp_vector_saturated_local.hpp>
#include <pairwise_aligner/matrix/dp_vector_saturated.hpp>
#include <pairwise_aligner/matrix/dp_vector_single.hpp>

#include <pairwise_aligner/score_model/score_model_matrix_simd_NxN.hpp>
#include <pairwise_aligner/tracker/tracker_global_simd_saturated.hpp>
#include <pairwise_aligner/tracker/tracker_local_simd_saturated.hpp>
#include <pairwise_aligner/type_traits.hpp>
#include <pairwise_aligner/utility/type_list.hpp>
#include <pairwise_aligner/simd/simd_score_type.hpp>

namespace seqan::pairwise_aligner {
inline namespace v1
{
namespace cfg
{
namespace _score_model_matrix_simd_saturated_NxN
{

// ----------------------------------------------------------------------------
// traits
// ----------------------------------------------------------------------------

template <typename substitution_matrix_t>
struct traits
{
    static constexpr cfg::detail::rule_category category = cfg::detail::rule_category::score_model;

    // extend the dimension to handle padding symbol.
    static constexpr size_t dimension = std::tuple_size_v<substitution_matrix_t> + 1;
    static constexpr size_t matrix_size = (dimension * (dimension + 1)) / 2;
    static constexpr size_t lane_width = 8;

    using matrix_row_t = typename substitution_matrix_t::value_type;
    using symbol_t = std::tuple_element_t<0, matrix_row_t>;
    using score_t = typename std::tuple_element_t<1, matrix_row_t>::value_type;
    using block_handler_t = detail::saturated_block_handler;

    // simd score typ: use full range?
    using score_type = simd_score<int8_t>;
    using index_type = seqan::pairwise_aligner::detail::make_unsigned_t<score_type>;

    using original_score_type = simd_score<score_t, score_type::size_v>;

    substitution_matrix_t _substitution_matrix;
    score_t _match_padding_score{1};
    score_t _mismatch_padding_score{-1};

    using score_model_type = score_model_matrix_simd_NxN<score_type, index_type, dimension>;

    template <typename cell_t>
    using buffer_t = std::vector<cell_t, seqan3::aligned_allocator<cell_t, alignof(cell_t)>>;

    template <typename dp_vector_t>
    using dp_vector_column_type = dp_vector_bulk<dp_vector_t, score_type>;

    template <typename dp_vector_t>
    using dp_vector_row_type = dp_vector_bulk<dp_vector_t, score_type>;

    template <typename configuration_t>
    constexpr auto configure_substitution_policy([[maybe_unused]] configuration_t const & configuration) const noexcept
    {
        // TODO: Run this mode!
        std::array<std::array<score_t, dimension>, dimension> tmp{};

        std::ranges::for_each(std::views::iota(size_t(0), dimension), [&] (size_t const i)
        {
            tmp[i].fill((configuration_t::is_local ? _mismatch_padding_score : _match_padding_score));
            // Overwrite the other values with data from the original substitution matrix.
            if (i < (dimension - 1))
                std::ranges::copy(_substitution_matrix[i].second, tmp[i].data());
        });

        return score_model_type{tmp};
    }

    template <typename configuration_t>
    constexpr auto configure_result_factory_policy([[maybe_unused]] configuration_t const & configuration)
        const noexcept
    {
        score_t max_match_score = std::numeric_limits<score_t>::lowest();
        score_t min_mismatch_score = std::numeric_limits<score_t>::max();
        std::ranges::for_each(_substitution_matrix | std::views::elements<1>, [&] (auto const & matrix_row) {
            auto [local_min, local_max] = std::ranges::minmax(matrix_row);
            min_mismatch_score = std::min(local_min, min_mismatch_score);
            max_match_score = std::max(local_max, max_match_score);
        });

        auto [global_zero, max_block_size] =
            block_handler_t::compute_max_block_size(max_match_score,
                                                    min_mismatch_score,
                                                    configuration._gap_open_score,
                                                    configuration._gap_extension_score);
        if constexpr (configuration_t::is_local) {
            return tracker::local_simd_saturated::factory<score_type, original_score_type>{};
        } else {
            return tracker::global_simd_saturated::factory<original_score_type>{
                static_cast<original_score_type>(_match_padding_score),
                configuration.trailing_gap_setting(),
                max_block_size
            };
        }
    }

    template <typename configuration_t>
    constexpr auto configure_dp_vector_policy([[maybe_unused]] configuration_t const & configuration) const noexcept
    {
        using column_cell_t = typename configuration_t::dp_cell_column_type<score_type>;
        using original_column_cell_t = typename configuration_t::dp_cell_column_type<original_score_type>;

        using row_cell_t = typename configuration_t::dp_cell_row_type<score_type>;
        using original_row_cell_t = typename configuration_t::dp_cell_row_type<original_score_type>;

        score_t max_match_score = std::numeric_limits<score_t>::lowest();
        score_t min_mismatch_score = std::numeric_limits<score_t>::max();
        std::ranges::for_each(_substitution_matrix | std::views::elements<1>, [&] (auto const & matrix_row) {
            auto [local_min, local_max] = std::ranges::minmax(matrix_row);
            min_mismatch_score = std::min(local_min, min_mismatch_score);
            max_match_score = std::max(local_max, max_match_score);
        });

        auto [global_zero, max_block_size] =
            block_handler_t::compute_max_block_size(max_match_score,
                                                    min_mismatch_score,
                                                    configuration._gap_open_score,
                                                    configuration._gap_extension_score);
        // Add padding symbol: Assuming the symbols are sorted lexicographically, take the last symbol plus one
        // (only char?).
        // TODO: Find some unused values between ranks/symbols!
        assert(static_cast<uint8_t>(_substitution_matrix.back().first) < 255);
        using scalar_rank_t = typename index_type::value_type;
        assert(static_cast<scalar_rank_t>(_substitution_matrix.back().first) < 255);
        scalar_rank_t const padding_symbol = _substitution_matrix.back().first + 1;
        std::vector<scalar_rank_t> extended_symbol_list{};
        extended_symbol_list.resize(dimension, padding_symbol);
        std::ranges::copy(_substitution_matrix | std::views::elements<0>, extended_symbol_list.begin());

        // Initialise the scale for the column sequence map.
        alphabet_rank_map_simd<index_type> rank_map{std::move(extended_symbol_list)};

        auto saturated_vector = [&] (auto original_cell_type, auto && base_vector)
        {
            using orginal_cell_t = typename decltype(original_cell_type)::type;
            using base_vector_t = decltype(base_vector);
            if constexpr (configuration_t::is_local) {
                int8_t local_zero = block_handler_t::lowest_viable_local_score(configuration._gap_open_score,
                                                                               configuration._gap_extension_score);
                int8_t threshold = std::numeric_limits<int8_t>::max() - (max_block_size * _match_padding_score);
                return dp_vector_saturated_local_factory<orginal_cell_t>(std::forward<base_vector_t>(base_vector),
                                                                         local_zero,
                                                                         global_zero,
                                                                         threshold);
            } else {
                return dp_vector_saturated_factory<orginal_cell_t>(std::forward<base_vector_t>(base_vector),
                                                                   global_zero);
            }
        };

        return dp_vector_policy{
                    dp_vector_bulk_factory(
                        dp_vector_rank_transformation_factory(
                            dp_vector_offset_transformation(
                                dp_vector_chunk_factory(
                                    saturated_vector(std::type_identity<original_column_cell_t>{},
                                                     dp_vector_single<column_cell_t>{}),
                                    max_block_size),
                                offset_transform{dimension, matrix_size}),
                            rank_map),
                        index_type{padding_symbol}),
                    dp_vector_bulk_factory(
                        dp_vector_rank_transformation_factory(
                            dp_vector_offset_transformation(
                                dp_vector_chunk_factory(
                                    saturated_vector(std::type_identity<original_row_cell_t>{},
                                                     dp_vector_single<row_cell_t>{}),
                                    max_block_size),
                                offset_transform{dimension, matrix_size}),
                            rank_map),
                        index_type{padding_symbol})
        };
    }

    template <typename configuration_t, typename ...policies_t>
    constexpr auto configure_algorithm(configuration_t const &, policies_t && ...policies) const noexcept
    {
        auto make_dp_matrix_policy = [&] () constexpr {
            if constexpr (configuration_t::is_local)
                return dp_matrix::matrix(dp_matrix::column_saturated_local(dp_matrix::block(dp_matrix::lane)));
            else
                return dp_matrix::matrix(dp_matrix::column_saturated(dp_matrix::block(dp_matrix::lane)));
        };
        // auto make_dp_matrix_policy = [&] () constexpr {
        //     if constexpr (configuration_t::is_local)
        //         return dp_matrix_column_saturated_local;
        //     else
        //         return dp_matrix_column_saturated;
        // };

        using dp_matrix_policy_t = dp_matrix_policies<std::invoke_result_t<decltype(make_dp_matrix_policy)>>;
        using algorithm_t = typename configuration_t::algorithm_type<dp_algorithm_template_standard,
                                                                     dp_matrix_policy_t,
                                                                     lane_width_policy<4>,
                                                                     std::remove_cvref_t<policies_t>...>;

        return interface_one_to_one_bulk<algorithm_t, score_type::size_v>{
                algorithm_t{dp_matrix_policy_t{make_dp_matrix_policy()},
                            lane_width_policy<4>{},
                            std::move(policies)...}};
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
    template <typename predecessor_t, typename alphabet_t, typename score_t, size_t dimension>
    constexpr auto operator()(predecessor_t && predecessor,
                              std::array<std::pair<alphabet_t, std::array<score_t, dimension>>, dimension>
                                    substitution_matrix) const
    {
        using substitution_matrix_t = std::array<std::pair<alphabet_t, std::array<score_t, dimension>>, dimension>;

        using traits_t = _score_model_matrix_simd_saturated_NxN::traits<substitution_matrix_t>;
        using rule_t = _score_model_matrix_simd_saturated_NxN::rule<predecessor_t, traits_t>;
        return rule_t{{}, std::forward<predecessor_t>(predecessor), traits_t{std::move(substitution_matrix)}};
    }

    template <typename alphabet_t, typename score_t, size_t dimension>
    constexpr auto operator()(std::array<std::pair<alphabet_t, std::array<score_t, dimension>>, dimension> const &
                                substitution_matrix)
        const
    {
        return this->operator()(cfg::initial, substitution_matrix);
    }
};

} // namespace _cpo
} // namespace _score_model

inline constexpr _score_model_matrix_simd_saturated_NxN::_cpo::_fn score_model_matrix_simd_saturated_NxN{};

} // namespace cfg
} // inline namespace v1
}  // namespace seqan::pairwise_aligner
