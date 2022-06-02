// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides option wrapper.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <concepts>
#include <type_traits>

#include <pairwise_aligner/align_configuration/is_option_set_cpo.hpp>
#include <pairwise_aligner/sender/tag_invoke.hpp>

namespace align
{

template <typename derived_t, typename option_value_t>
    requires std::is_enum_v<option_value_t>
class enum_option {
private:

    friend derived_t;

    using underlying_t = std::underlying_type_t<option_value_t>;
    option_value_t _value;

    enum_option() = default;
    // implicit construction from option
    constexpr explicit enum_option(option_value_t value) noexcept : _value{value}
    {}

public:

    using value_type = option_value_t;

    constexpr derived_t operator|(enum_option const & rhs) const {
        return derived_t{static_cast<underlying_t const &>(_value) | static_cast<underlying_t const &>(rhs._value)};
    }

private:
    constexpr friend bool tag_invoke(tag_t<align::is_option_set>, enum_option const & me, value_type const & value)
        noexcept
    {
        return static_cast<underlying_t>(me._value) & static_cast<underlying_t>(value);
    }
};

} // namespace align
