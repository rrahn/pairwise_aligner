// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::cfg::rule_base.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <seqan3/std/type_traits>

#include <pairwise_aligner/configuration/rule_category.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace cfg
{

struct rule_base
{};

namespace _base
{

template <typename rule_t, cfg::detail::rule_category category>
struct _rule
{
    struct type;
};

template <typename rule_t, cfg::detail::rule_category category>
using rule = typename _rule<rule_t, category>::type;

template <typename rule_t, cfg::detail::rule_category category>
struct _rule<rule_t, category>::type : rule_base
{
    using pure_rule_t = std::remove_cvref_t<rule_t>;
    inline static constexpr uint8_t _index = static_cast<uint8_t>(category);
    inline static constexpr bool already_applied = pure_rule_t::configurator_index_of[_index] != -1;

    using configurator_indices_t = typename pure_rule_t::configurator_indices_t;

    inline constexpr static int8_t configurator_index = pure_rule_t::configurator_index + 1;

    inline constexpr static configurator_indices_t configurator_index_of{[] (auto configs) -> configurator_indices_t
    {
        configs[_index] = configurator_index;
        return configs;
    }(pure_rule_t::configurator_index_of)};
};
} // namespace _base
} // namespace cfg
} // inline namespace v1
} // namespace seqan::pairwise_aligner
