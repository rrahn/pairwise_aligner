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

#include <type_traits>

#include <seqan3/utility/type_pack/traits.hpp>

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

// concat_type_lists_unique<unique_lists_t...>
//
// Result is produced via '::type' which will contain a type_list<Ts...> that
// contains the unique elements from the input type_list types.
// Assumes that the input lists already
template <typename ...unique_lists_t>
struct concat_type_lists_unique;

template <>
struct concat_type_lists_unique<> {
    using type = type_list<>;
};

template <typename unique_list_t>
struct concat_type_lists_unique<unique_list_t> {
    using type = unique_list_t;
};

template <typename... Ts, typename... Us, typename... other_lists_t>
struct concat_type_lists_unique<type_list<Ts...>, type_list<Us...>, other_lists_t...>
    : concat_type_lists_unique<
        typename concat_type_lists<type_list<Ts...>,
                                   std::conditional_t<seqan3::pack_traits::contains<Us, Ts...>,
                                                      type_list<>,
                                                      type_list<Us>>...>::type,
        other_lists_t...>
{};

template <typename... unique_lists_t>
using concat_type_lists_unique_t = typename concat_type_lists_unique<unique_lists_t...>::type;

namespace detail
{
template <template <typename...> class Outer,
          template <typename...> class Inner>
struct type_list_nested_apply_impl {
    template <typename... lists_t>
    using apply = Outer<typename lists_t::template apply<Inner>...>;
};

} // namespace detail

template <typename list_of_lists_t,
          template <typename...> class Outer,
          template <typename...> class Inner>
using type_list_nested_apply_t =
    typename list_of_lists_t::template apply<
        detail::type_list_nested_apply_impl<Outer, Inner>::template apply>;

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
