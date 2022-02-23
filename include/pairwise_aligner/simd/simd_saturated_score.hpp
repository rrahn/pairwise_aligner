// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::saturated_score.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <pairwise_aligner/simd/add_saturated.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace detail {
template <typename derived_t>
struct _saturated_score
{
    class type;
};

template <typename derived_t>
class _saturated_score<derived_t>::type
{
public:

    template <typename _derived_t = derived_t>
        requires requires (typename _derived_t::simd_type::value_type const & v) {
            { simd::add_saturated(v, v) } -> std::same_as<typename _derived_t::simd_type::value_type>;
        }
    constexpr _derived_t add(_derived_t const & rhs) const noexcept
    {
        _derived_t tmp{};
        as_derived().apply([] (auto & result, auto const & a, auto const & b) {
            result = simd::add_saturated(a, b);
        }, tmp.values, as_derived().values, rhs.values);
        return tmp;
    }

    template <typename _derived_t = derived_t>
        requires requires (typename _derived_t::simd_type::value_type const & v) {
            { simd::add_saturated(v, v) } -> std::same_as<typename _derived_t::simd_type::value_type>;
        }
    constexpr _derived_t add(typename _derived_t::value_type const & rhs) const noexcept
    {
        _derived_t tmp{};
        as_derived().apply([] (auto & result, auto const & a, auto const & b) {
            result = simd::add_saturated(a, b);
        }, tmp.values, as_derived().values, _derived_t{rhs}.values);
        return tmp;
    }

private:

    template <typename t = void>
    constexpr derived_t const & as_derived() const noexcept
    {
        return static_cast<derived_t const &>(*this);
    }

};
} // namespace detail

template <typename derived_t>
using saturated_score = typename detail::_saturated_score<derived_t>::type;

} // inline namespace v1
} // namespace seqan::pairwise_aligner
