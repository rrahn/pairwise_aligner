// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides various helper templates to access type properties.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <cstdint>

#include <seqan3/std/concepts>
#include <seqan3/std/type_traits>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

/*!\brief The order of the accessed vector.
 *
 * This enum is used to disambiguate cells from a column or a row vector of the underlying DP matrix.
 */
enum struct dp_vector_order
{
    row,
    column
};

//!\brief To check if a type is a dp cell.
template <dp_vector_order _order>
struct dp_cell_base
{};

// CPO?
template <typename dp_cell_t>
inline constexpr bool is_row_cell_v = false;

template <std::derived_from<dp_cell_base<dp_vector_order::row>> dp_cell_t>
inline constexpr bool is_row_cell_v<dp_cell_t> = true;

template <typename dp_cell_t>
inline constexpr bool is_column_cell_v = false;

template <std::derived_from<dp_cell_base<dp_vector_order::column>> dp_cell_t>
inline constexpr bool is_column_cell_v<dp_cell_t> = true;

template <typename rule_t, template <typename...> typename tuple_t>
using configurator_types_t = typename rule_t::template configurator_types<tuple_t>;

template <typename t>
struct remove_rvalue_reference : std::type_identity<t>
{};

template <typename t>
    requires (std::is_rvalue_reference_v<t>)
struct remove_rvalue_reference<t> : std::type_identity<std::remove_reference_t<t>>
{};

template <typename t>
using remove_rvalue_reference_t = typename remove_rvalue_reference<t>::type;

template <template <typename> typename deferred_t>
struct lazy_type
{
    template <typename ...types>
    using type = deferred_t<types...>;
};

template <typename lazy_type, typename ...types>
using instantiate_t = typename lazy_type::type<types...>;

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
