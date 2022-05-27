// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides .
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

namespace seqan::align
{
inline namespace v1
{

// More thought?
// Attribute CPOs?
// inline constexpr query_fn {
//     template <typename property_t, typename configuration_t>
//     constexpr auto operator()(property_t const &, configuration_t && configuration) const noexcept
//         ->
//     {
//         return tag_invoke(property_t, )
//     }
// } query;




template <typename configuration_t>
struct dp_traits{
    using score_type = int; // decltype(align::query(dp_attribute<substituion_cost>, std::declval<configuration_t>()));

};

template <typename configuration_t>
using dp_traits_t = dp_traits<configuration_t>;

} // inline namespace v1
} // namespace seqan::align
