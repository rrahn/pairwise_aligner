// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::cfg::score_model_unitary_simd.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <seqan3/std/type_traits>

#include <pairwise_aligner/configuration/initial.hpp>
#include <pairwise_aligner/configuration/rule_score_model.hpp>
#include <pairwise_aligner/pairwise_aligner.hpp>
#include <pairwise_aligner/interface/interface_one_to_one_bulk.hpp>
#include <pairwise_aligner/score_model/score_model_unitary.hpp>
#include <pairwise_aligner/type_traits.hpp>
#include <pairwise_aligner/utility/type_list.hpp>
#include <pairwise_aligner/simd_score_type.hpp>

namespace seqan::pairwise_aligner {
inline namespace v1
{
namespace cfg
{
namespace _score_model_unitary_simd
{

// ----------------------------------------------------------------------------
// traits
// ----------------------------------------------------------------------------

template <typename score_t>
struct traits
{
    score_t _match_score;
    score_t _mismatch_score;

    // Offer the score type here.
    using score_type = simd_score<score_t>;

    using score_model_type = seqan::pairwise_aligner::score_model_unitary<score_type>;

    // Offer some overload for the column type.
    template <typename dp_cell_t>
    using dp_vector_column_type = simd_intermediate_dp_vector<dp_cell_t>;

    // Offer some overload for the column type.
    template <typename dp_cell_t>
    using dp_vector_row_type = simd_intermediate_dp_vector<dp_cell_t>;

    template <typename dp_algorithm_t, typename dp_vector_column_t, typename dp_vector_row_t>
    using dp_interface_type = interface_one_to_one_bulk<dp_algorithm_t, dp_vector_column_t, dp_vector_row_t>;

    constexpr score_model_type create() const
    {
        return score_model_type{static_cast<score_type>(_match_score), static_cast<score_type>(_mismatch_score)};
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
    void set_config(values_t && ... values) && noexcept
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
                configurator_t<std::remove_cvref_t<next_configurator_t>, traits_t>{
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
        return _score_model_unitary_simd::rule<predecessor_t, traits_t>{{},
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
} // namespace _score_model

inline constexpr _score_model_unitary_simd::_cpo::_fn score_model_unitary_simd{};

} // namespace cfg
} // inline namespace v1
}  // namespace seqan::pairwise_aligner
