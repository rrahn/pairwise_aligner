// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::cfg::method_global.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <seqan3/std/type_traits>

#include <pairwise_aligner/affine/affine_initialisation_strategy.hpp>
#include <pairwise_aligner/configuration/initial.hpp>
#include <pairwise_aligner/configuration/rule_method.hpp>
#include <pairwise_aligner/dp_trailing_gaps.hpp>
#include <pairwise_aligner/type_traits.hpp>
#include <pairwise_aligner/utility/type_list.hpp>

namespace seqan::pairwise_aligner {
inline namespace v1
{
namespace cfg
{
namespace _method_global
{

// ----------------------------------------------------------------------------
// traits
// ----------------------------------------------------------------------------

struct traits
{
    static constexpr cfg::detail::rule_category category = cfg::detail::rule_category::method;

    initialisation_rule _leading_param;
    trailing_gap_setting _trailing_param;

    constexpr std::pair<initialisation_rule, trailing_gap_setting> create() const noexcept
    {
        return std::pair{std::move(_leading_param), std::move(_trailing_param)};
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
struct _rule<predecessor_t, traits_t>::type : cfg::method::rule<predecessor_t>
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
    // implementation of function style connection
    template <typename predecessor_t>
    auto operator()(predecessor_t && predecessor, initialisation_rule init_rule, trailing_gap_setting final_rule) const
    {
        return _method_global::rule<predecessor_t, traits>{{},
                                                           std::forward<predecessor_t>(predecessor),
                                                           traits{init_rule, final_rule}};
    }

    template <typename score_t>
    auto operator()(initialisation_rule init_rule, trailing_gap_setting final_rule) const
    {
        return this->operator()(cfg::initial, std::move(init_rule), std::move(final_rule));
    }
};
} // namespace _cpo
} // namespace _method_global

inline constexpr _method_global::_cpo::_fn method_global{};

} // namespace cfg
} // inline namespace v1
}  // namespace seqan::pairwise_aligner
