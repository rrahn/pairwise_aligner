// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::cfg::_initial::rule.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <cstdint>
#include <array>

#include <pairwise_aligner/configuration/rule_base.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace cfg
{

namespace _initial
{

// ----------------------------------------------------------------------------
// rule
// ----------------------------------------------------------------------------

struct rule : rule_base
{
    using configurator_indices_t = std::array<int8_t, static_cast<uint8_t>(detail::rule_category::size)>;
    inline constexpr static configurator_indices_t configurator_index_of = [] () -> configurator_indices_t
    {
        configurator_indices_t tmp{};
        tmp.fill(-1);
        return tmp;
    }();

    inline constexpr static int8_t configurator_index{-1};

    template <template <typename ...> typename tuple_t>
    using configurator_types = tuple_t<>;

    template <typename configurator_t>
    struct _op
    {
        configurator_t _configurator;

        void start() noexcept
        {
            std::forward<configurator_t>(_configurator).set_config();
        }
    };

    template <typename configurator_t>
    auto apply(configurator_t && configurator) const
        -> _op<configurator_t>
    {
        return _op<configurator_t>{std::forward<configurator_t>(configurator)};
    }
};
} // namespace _initial

// ----------------------------------------------------------------------------
// CPO
// ----------------------------------------------------------------------------

using initial_t = _initial::rule;
inline constexpr initial_t initial{};

} // namespace cfg
} // inline namespace v1
}  // namespace seqan::pairwise_aligner
