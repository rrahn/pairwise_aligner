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

#include <type_traits>

namespace seqan::align
{
inline namespace v1
{

// get the attribute for the passed attribute tag
template <auto & attribute_tag_t>
using attribute = std::decay_t<decltype(attribute_tag_t)>;

template <typename tag_t, typename tag_value_t>
struct attribute_tag {
private:
    tag_value_t _value;

    constexpr friend auto tag_invoke(tag_t<query>, attribute_tag & me) noexcept
        -> tag_value_t &
    {
        return _value;
    }

    constexpr friend void tag_invoke(tag_t<set>, attribute_tag & me, tag_value_t value) noexcept
    {
        _value = std::move(value);
    }
};

// we need to define something that we can get the attribute for later!
template <typename gap_cost_model_t>
inline constexpr attribute_tag<struct gap_cost_tag, gap_cost_model_t> gap_cost;



// gap cost attribute:
// template<std::size_t N>
// struct AttributeTag
// {
//     char _tag_name[N+1]{};

//     constexpr AttributeTag(char const(&tag_name)[N])
//     {
//         std::ranges::copy(tag_name, _attribute_tag);
//         // what do we want to customise?
//     };
// };

// template <AttributeTag tag>
// auto operator"" _attribute ()
// {   // what do we want to return?
//     return
// }

// auto attr = "GapCost"_attribute; // now we have a gap cost attribute
// tag_t<"GapCost"_attribute> // this binds the tag_t in a value.
// // So I want to query: query("GapCost"_attribute, configuration);

// template <auto & ...attributes>
// struct Configuraition : std::tuple<>{

//     template <>
//     constexpr friend auto tag_invoke()
// };
    // but what can it do and how is it specialised?

} // inline namespace v1
} // namespace seqan::align
