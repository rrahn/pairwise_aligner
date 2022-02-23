// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::simd::add_saturated.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <immintrin.h>

#include <seqan3/std/concepts>
#include <seqan3/std/type_traits>

#include <pairwise_aligner/simd/host_simd_tag.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace simd {
namespace detail {

struct add_saturated
{
    template <typename simd_t, size_t size, size_t width>
        requires ((size == width) || (size == (width >> 1)))
    simd_t operator()(simd_t const & a, simd_t const & b, host_simd_tag<size, width> const & tag) const noexcept
    {
        using scalar_t = typename seqan3::simd::simd_traits<simd_t>::scalar_type;
        constexpr bool is_signed = std::is_signed_v<scalar_t>;
        if constexpr (is_signed)
            return add_saturated_signed(a, b, tag);
        else
            return add_saturated_unsigned(a, b, tag);
    }

private:

    // ----------------------------------------------------------------------------
    // Signed saturated add
    // ----------------------------------------------------------------------------

    template <typename simd_t, size_t size>
    simd_t add_saturated_signed(simd_t const & a, simd_t const & b, host_simd_tag<size, 64> const & tag) const noexcept
    {
        if constexpr (size == 64)
            return to_packed<simd_t>(_mm512_adds_epi8(to_native(a, tag), to_native(b, tag)));
        else if constexpr (size == 32)
            return to_packed<simd_t>(_mm512_adds_epi16(to_native(a, tag), to_native(b, tag)));
    }

    template <typename simd_t, size_t size>
    simd_t add_saturated_signed(simd_t const & a, simd_t const & b, host_simd_tag<size, 32> const & tag) const noexcept
    {
        if constexpr (size == 32)
            return to_packed<simd_t>(_mm256_adds_epi8(to_native(a, tag), to_native(b, tag)));
        else if constexpr (size == 16)
            return to_packed<simd_t>(_mm256_adds_epi16(to_native(a, tag), to_native(b, tag)));
    }

    template <typename simd_t, size_t size>
    simd_t add_saturated_signed(simd_t const & a, simd_t const & b, host_simd_tag<size, 16> const & tag) const noexcept
    {
        if constexpr (size == 16)
            return to_packed<simd_t>(_mm_adds_epi8(to_native(a, tag), to_native(b, tag)));
        else if constexpr (size == 8)
            return to_packed<simd_t>(_mm_adds_epi16(to_native(a, tag), to_native(b, tag)));
    }

    // ----------------------------------------------------------------------------
    // Unsigned saturated add
    // ----------------------------------------------------------------------------

    template <typename simd_t, size_t size>
    simd_t add_saturated_unsigned(simd_t const & a, simd_t const & b, host_simd_tag<size, 64> const & tag) const noexcept
    {
        if constexpr (size == 64)
            return to_packed<simd_t>(_mm512_adds_epu8(to_native(a, tag), to_native(b, tag)));
        else if constexpr (size == 32)
            return to_packed<simd_t>(_mm512_adds_epu16(to_native(a, tag), to_native(b, tag)));
    }

    template <typename simd_t, size_t size>
    simd_t add_saturated_unsigned(simd_t const & a, simd_t const & b, host_simd_tag<size, 32> const & tag) const noexcept
    {
        if constexpr (size == 32)
            return to_packed<simd_t>(_mm256_adds_epu8(to_native(a, tag), to_native(b, tag)));
        else if constexpr (size == 16)
            return to_packed<simd_t>(_mm256_adds_epu16(to_native(a, tag), to_native(b, tag)));
    }

    template <typename simd_t, size_t size>
    simd_t add_saturated_unsigned(simd_t const & a, simd_t const & b, host_simd_tag<size, 16> const & tag) const noexcept
    {
        if constexpr (size == 16)
            return to_packed<simd_t>(_mm_adds_epu8(to_native(a, tag), to_native(b, tag)));
        else if constexpr (size == 8)
            return to_packed<simd_t>(_mm_adds_epu16(to_native(a, tag), to_native(b, tag)));
    }

private:

    template <typename simd_t, size_t size, size_t width>
    constexpr auto const & to_native(simd_t const & packed, host_simd_tag<size, width> const &) const noexcept
    {
        if constexpr (width == 64)
            return reinterpret_cast<__m512i const &>(packed);
        else if constexpr (width == 32)
            return reinterpret_cast<__m256i const &>(packed);
        else
            return reinterpret_cast<__m128i const &>(packed);
    }

    template <typename simd_t>
    constexpr auto const & to_packed(auto const & native) const noexcept
    {
        return reinterpret_cast<simd_t const &>(native);
    }
};

} // namespace detail

inline namespace _cpo {

inline constexpr auto add_saturated = [] (auto const & a, auto const & b)
    -> std::invoke_result_t<detail::add_saturated,
                            decltype(a),
                            decltype(b),
                            detail::host_simd_tag_t<std::common_type_t<std::remove_cvref_t<decltype(a)>,
                                                                       std::remove_cvref_t<decltype(b)>>>>
    requires std::common_with<std::remove_cvref_t<decltype(a)>, std::remove_cvref_t<decltype(b)>>
{
    using simd_t = std::common_type_t<std::remove_cvref_t<decltype(a)>, std::remove_cvref_t<decltype(b)>>;

    return std::invoke(detail::add_saturated{}, a, b, detail::host_simd_tag_t<simd_t>{});
};
} // inline namespace _cpo
} // namespace simd
} // inline namespace v1
} // namespace seqan::pairwise_aligner
