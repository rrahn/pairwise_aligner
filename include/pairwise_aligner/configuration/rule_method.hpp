// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::cfg::gap_model::rule.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <seqan3/std/type_traits>

#include <pairwise_aligner/configuration/rule_base.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace cfg::method
{

template <typename rule_t>
struct _rule
{
    struct type;
};

template <typename rule_t>
using rule = typename _rule<rule_t>::type;

template <typename rule_t>
struct _rule<rule_t>::type : _base::rule<rule_t, cfg::detail::rule_category::method>
{
    using rule_base_t = _base::rule<rule_t, cfg::detail::rule_category::method>;
    static_assert(!rule_base_t::already_applied, "The method category was already configured by another rule!");
};
} // namespace cfg::method
} // inline namespace v1
} // namespace seqan::pairwise_aligner
