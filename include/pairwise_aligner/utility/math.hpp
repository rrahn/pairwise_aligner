// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::add.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <pairwise_aligner/utility/priority_tag.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace _add {

constexpr auto add(...) noexcept = delete;

struct _fn
{
    template <typename t1, typename t2>
    constexpr auto operator()(t1 && a, t2 && b, pairwise_aligner::detail::priority_tag<0> const &)
        const noexcept(noexcept(a + b))
        -> decltype(a + b)
    {
        return a + b;
    }

    template <typename t1, typename t2>
    constexpr auto operator()(t1 && a, t2 && b, pairwise_aligner::detail::priority_tag<1> const &)
        const noexcept(noexcept(add(a, b)))
        -> decltype (add(a, b))
    {
        return add(a, b);
    }

    template <typename t1, typename t2>
    constexpr auto operator()(t1 && a, t2 && b, pairwise_aligner::detail::priority_tag<2> const &)
        const noexcept(noexcept(a.add(b)))
        -> decltype (a.add(b))
    {
        return a.add(b);
    }
};

} // namespace _add

inline namespace _cpo
{

inline constexpr auto add = [] (auto const & a, auto const & b)
    noexcept(noexcept(std::invoke(_add::_fn{}, a, b, pairwise_aligner::detail::priority_tag<2>{})))
    -> std::invoke_result_t<_add::_fn, decltype(a), decltype(b), pairwise_aligner::detail::priority_tag<2>>
{
    return std::invoke(_add::_fn{}, a, b, pairwise_aligner::detail::priority_tag<2>{});
};
} // inline namespace _cpo

namespace _subtract {

constexpr auto subtract(...) noexcept = delete;

struct _fn
{
    template <typename t1, typename t2>
    constexpr auto operator()(t1 && a, t2 && b, pairwise_aligner::detail::priority_tag<0> const &)
        const noexcept(noexcept(a - b))
        -> decltype(a - b)
    {
        return a - b;
    }

    template <typename t1, typename t2>
    constexpr auto operator()(t1 && a, t2 && b, pairwise_aligner::detail::priority_tag<1> const &)
        const noexcept(noexcept(subtract(a, b)))
        -> decltype (subtract(a, b))
    {
        return subtract(a, b);
    }

    template <typename t1, typename t2>
    constexpr auto operator()(t1 && a, t2 && b, pairwise_aligner::detail::priority_tag<2> const &)
        const noexcept(noexcept(a.subtract(b)))
        -> decltype (a.subtract(b))
    {
        return a.subtract(b);
    }
};

} // namespace _subtract

inline namespace _cpo
{

inline constexpr auto subtract = [] (auto const & a, auto const & b)
    noexcept(noexcept(std::invoke(_subtract::_fn{}, a, b, pairwise_aligner::detail::priority_tag<2>{})))
    -> std::invoke_result_t<_subtract::_fn, decltype(a), decltype(b), pairwise_aligner::detail::priority_tag<2>>
{
    return std::invoke(_subtract::_fn{}, a, b, pairwise_aligner::detail::priority_tag<2>{});
};
} // inline namespace _cpo
} // inline namespace v1
} // namespace seqan::pairwise_aligner
