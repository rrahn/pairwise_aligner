// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides configurations for alignment.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <concepts>
#include <type_traits>

#include <pairwise_aligner/align_configuration/enum_option.hpp>
#include <pairwise_aligner/sender/tag_invoke.hpp>

namespace align {
// must be something more rich so we can add pipeable configuration.
enum class alignment_kind_option {
    unknown = 1,
    score = 2,
    positions = 4,
    alignment = 8,
};

class alignment_kind;

namespace _get_alignment_kind {

inline constexpr struct _fn {
    template<typename align_config_t>
        requires tag_invocable<_fn, align_config_t>
    constexpr auto operator()(align_config_t && configuration) const
        noexcept(is_nothrow_tag_invocable_v<_fn, align_config_t>)
        -> tag_invoke_result_t<_fn, align_config_t> {
        return tag_invoke(_fn{}, std::forward<align_config_t>(configuration));
    }

    // default implementation
    template <typename align_config_t>
        requires (!tag_invocable<_fn, align_config_t>)
    constexpr alignment_kind_option operator()(align_config_t &&) const noexcept {
        return alignment_kind_option::unknown;
    }

} get_alignment_kind{};

} // namespace _get_alignment_kind

using _get_alignment_kind::get_alignment_kind;

template <typename source_t, typename target_t>
struct adopt_reference_from : public std::type_identity<target_t>
{};

template <typename source_t, typename target_t>
    requires std::is_lvalue_reference_v<source_t>
struct adopt_reference_from<source_t, target_t> : public std::add_lvalue_reference<std::remove_reference_t<target_t>>
{};

template <typename source_t, typename target_t>
    requires std::is_rvalue_reference_v<source_t>
struct adopt_reference_from<source_t, target_t> : public std::remove_reference_t<target_t>
{};

template <typename source_t, typename target_t>
using adopt_reference_from_t = typename adopt_reference_from<source_t, target_t>::type;

class alignment_kind : public enum_option<alignment_kind, alignment_kind_option> {
private:
    using base_t = enum_option<alignment_kind, alignment_kind_option>;

public:
    using value_type = alignment_kind_option;

    alignment_kind() : base_t{alignment_kind_option::unknown}
    {}
    constexpr explicit alignment_kind(alignment_kind_option option) noexcept : base_t{option}
    {}
    constexpr alignment_kind(base_t base) noexcept : base_t{base}
    {}

private:

    template <typename alignment_kind_t>
        requires std::same_as<std::remove_cvref_t<alignment_kind_t>, alignment_kind>
    constexpr friend alignment_kind_t tag_invoke(tag_t<get_alignment_kind>, alignment_kind_t && me) noexcept {
        return std::forward<alignment_kind_t>(me);
    }

    constexpr friend bool tag_invoke(tag_t<align::is_option_set>, alignment_kind const & me, value_type const & value)
        noexcept {
        return align::is_option_set(static_cast<base_t const &>(me), value);
    }

    // template <typename cpo_t, typename alignment_kind_t, typename ...args_t>
    //     requires std::same_as<std::remove_cvref_t<alignment_kind_t>, alignment_kind>
    // constexpr friend auto tag_invoke(cpo_t cpo, alignment_kind_t && me, args_t && ...args)
    //     noexcept(noexcept(cpo(std::forward<alignment_kind_t>(me), std::forward<args_t>(args)...)))
    //     -> std::invoke_result_t<cpo_t, adopt_reference_from_t<alignment_kind_t, base_t>, args_t...> {
    //     return cpo(static_cast<adopt_reference_from_t<alignment_kind_t, base_t>>(me), std::forward<args_t>(args)...);
    // }
};

} // namespace align
