// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::cfg::score_model_unitary.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <seqan3/std/type_traits>

#include <pairwise_aligner/configuration/initial.hpp>
#include <pairwise_aligner/pairwise_aligner.hpp>
#include <pairwise_aligner/interface/interface_one_to_one_single.hpp>
#include <pairwise_aligner/score_model/score_model_unitary.hpp>
#include <pairwise_aligner/type_traits.hpp>
#include <pairwise_aligner/utility/type_list.hpp>
#include <pairwise_aligner/configuration/rule_score_model.hpp>

namespace seqan::pairwise_aligner {
inline namespace v1
{
namespace cfg
{
namespace _score_model_unitary
{

// ----------------------------------------------------------------------------
// traits
// ----------------------------------------------------------------------------

template <typename score_t>
struct traits
{
    score_t _match_score;
    score_t _mismatch_score;

    using score_model_type = score_model_unitary<score_t>;

    // Offer the score type here.
    using score_type = score_t;

    // Offer some overload for the column type.
    template <typename dp_cell_t>
    using dp_vector_column_type = intermediate_dp_vector<dp_cell_t>;

    // Offer some overload for the column type.
    template <typename dp_cell_t>
    using dp_vector_row_type = intermediate_dp_vector<dp_cell_t>;

    template <typename dp_algorithm_t, typename dp_vector_column_t, typename dp_vector_row_t>
    using dp_interface_type = interface_one_to_one_single<dp_algorithm_t, dp_vector_column_t, dp_vector_row_t>;

    constexpr score_model_type create() const
    {
        return score_model_type{_match_score, _mismatch_score};
    }
};

// ----------------------------------------------------------------------------
// configurator
// ----------------------------------------------------------------------------

template <typename next_configurator_t, typename score_model_t>
struct _configurator
{
    struct type;
};

template <typename next_configurator_t, typename score_model_t>
using configurator_t = typename _configurator<next_configurator_t, score_model_t>::type;

template <typename next_configurator_t, typename score_model_t>
struct _configurator<next_configurator_t, score_model_t>::type
{
    next_configurator_t _next_configurator;
    score_model_t _score_model;

    template <typename ...values_t>
    void set_config(values_t && ... values) && noexcept
    {
        std::forward<next_configurator_t>(_next_configurator).set_config(std::forward<values_t>(values)...,
                                                                         std::forward<score_model_t>(_score_model));
    }
};

// ----------------------------------------------------------------------------
// rule
// ----------------------------------------------------------------------------

template <typename predecessor_t, typename score_model_t>
struct _rule
{
    struct type;
};

template <typename predecessor_t, typename score_model_t>
using rule = typename _rule<predecessor_t, score_model_t>::type;

template <typename predecessor_t, typename score_model_t>
struct _rule<predecessor_t, score_model_t>::type : cfg::score_model::rule<predecessor_t>
{
    predecessor_t _predecessor;
    score_model_t _score_model;

    using score_t = typename score_model_t::score_type;
    using traits_type = type_list<traits<score_t>>;

    template <template <typename ...> typename type_list_t>
    using configurator_types = typename concat_type_lists_t<configurator_types_t<std::remove_cvref_t<predecessor_t>,
                                                                                 type_list>,
                                                            traits_type>::template apply<type_list_t>;

    template <typename next_configurator_t>
    auto apply(next_configurator_t && next_configurator)
    {
        return std::forward<predecessor_t>(_predecessor).apply(
                configurator_t<std::remove_cvref_t<next_configurator_t>, score_model_t>{
                                std::forward<next_configurator_t>(next_configurator),
                                std::forward<score_model_t>(_score_model)});
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
        using score_model_t = pairwise_aligner::score_model_unitary<score_t>;
        return _score_model_unitary::rule<predecessor_t, score_model_t>{{},
                                                                        std::forward<predecessor_t>(predecessor),
                                                                        score_model_t{match_score, mismatch_score}};
    }

    template <typename score_t>
    auto operator()(score_t const match_score, score_t const mismatch_score) const
    {
        return this->operator()(cfg::initial, match_score, mismatch_score);
        return this->operator()(cfg::initial, match_score, mismatch_score);
    }
};
} // namespace _cpo
} // namespace _score_model

inline constexpr _score_model_unitary::_cpo::_fn score_model_unitary{};

} // namespace cfg
} // inline namespace v1
}  // namespace seqan::pairwise_aligner
