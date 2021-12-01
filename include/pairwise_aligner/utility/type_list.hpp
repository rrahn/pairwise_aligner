// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::type_list and utility type traits.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename... types>
struct type_list
{
    // Invoke the template metafunction with the type-list's elements
    // as arguments.
    template <template <typename...> typename function_t>
    using apply = function_t<types...>;
};

template <typename... list_types>
struct concat_type_lists;

template <>
struct concat_type_lists<>
{
    using type = type_list<>;
};

template <typename list_type>
struct concat_type_lists<list_type>
{
    using type = list_type;
};

template <typename... first_list_types, typename... second_list_types>
struct concat_type_lists<type_list<first_list_types...>, type_list<second_list_types...>>
{
    using type = type_list<first_list_types..., second_list_types...>;
};

template <typename... list_types>
using concat_type_lists_t = typename concat_type_lists<list_types...>::type;

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
