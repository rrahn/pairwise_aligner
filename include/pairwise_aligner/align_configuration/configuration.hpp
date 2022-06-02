// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides align configuration object.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <tuple>
#include <type_traits>

#include <seqan3/utility/type_pack/traits.hpp>
#include <pairwise_aligner/utility/pack_algorithms.hpp>
#include <pairwise_aligner/sender/tag_invoke.hpp>

namespace align {

template <typename ...options_t>
struct configuration : public std::tuple<options_t...> {
private:
    using base_t = std::tuple<options_t...>;

    template <typename cpo_t>
    struct is_option_set {
        template <typename option_t>
        static constexpr bool invoke = tag_invocable<cpo_t, option_t>;
    };

public:

    using base_t::base_t;

    template <typename cpo_t, typename _configuration_t>
        requires std::same_as<std::remove_cvref_t<_configuration_t>, configuration>
    friend constexpr auto tag_invoke(cpo_t cpo, _configuration_t && me) noexcept
    {
        using config_t = seqan::pairwise_aligner::find_if_t<is_option_set<cpo_t>, options_t...>;
        return cpo(get<config_t>(me));
    }
};

template <typename ...args_t>
configuration(args_t && ...) -> configuration<args_t...>;

} // namespace align

namespace std {

template <typename ...options_t>
struct tuple_size<align::configuration<options_t...>> :
    public std::integral_constant<size_t, sizeof...(options_t)>
{};

template <size_t index, typename ...options_t>
struct tuple_element<index, align::configuration<options_t...>> :
    public tuple_element<index, std::tuple<options_t...>>
{};

} // namespace std
