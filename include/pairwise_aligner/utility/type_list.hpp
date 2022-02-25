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

#include <seqan3/utility/type_list/type_list.hpp>
#include <seqan3/utility/type_list/traits.hpp>
#include <seqan3/utility/type_pack/traits.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
template <typename...types>
struct type_list
{
    // Invoke the template metafunction with the type-list's elements
    // as arguments.
    template <template <typename...> typename function_t>
    using apply = function_t<types...>;
};

namespace detail {
template <typename ...types>
struct apply_type_list;

template <template <typename ...> typename list_t, typename ...types>
struct apply_type_list<list_t<types...>>
{
    template <template <typename...> typename function_t>
    using type = typename type_list<types...>::apply<function_t>;
};

} // namespace detail

template <template <typename...> typename function_t, typename list_t>
using apply_t = typename detail::apply_type_list<list_t>::type<function_t>;

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

template <template <typename ...> typename first_list_t,
          template <typename ...> typename second_list_t,
          typename... first_list_types,
          typename... second_list_types>
struct concat_type_lists<first_list_t<first_list_types...>, second_list_t<second_list_types...>>
{
    using type = type_list<first_list_types..., second_list_types...>;
};

template <typename... list_types>
using concat_type_lists_t = typename concat_type_lists<list_types...>::type;

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
