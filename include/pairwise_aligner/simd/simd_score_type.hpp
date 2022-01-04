// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise::simd_score_type.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <array>
#include <seqan3/std/concepts>

#include <seqan3/utility/simd/algorithm.hpp>
#include <seqan3/utility/simd/simd.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

namespace detail
{
#if defined(__AVX512F__)
inline constexpr size_t max_simd_size = 64;
#elif defined(__AVX2__)
inline constexpr size_t max_simd_size = 32;
#else
inline constexpr size_t max_simd_size = 16;
#endif

} // namespace detail

template <std::integral score_t, size_t simd_size = detail::max_simd_size/sizeof(score_t)>
class alignas(detail::max_simd_size) simd_score
{
public:

    inline static constexpr size_t size = simd_size;

    using simd_type = seqan3::simd::simd_type_t<score_t, simd_size>;
    using value_type = score_t;
    using reference = score_t &;
    using const_reference = score_t const &;

private:

    static constexpr bool is_native = simd_size == detail::max_simd_size/sizeof(score_t);

    simd_type values{};

public:

    simd_score() = default;
    simd_score(simd_score const &) = default;
    simd_score(simd_score &&) = default;
    simd_score & operator=(simd_score const &) = default;
    simd_score & operator=(simd_score &&) = default;

    explicit simd_score(score_t const initial_score) noexcept
    {
        if constexpr (is_native) {
            values = seqan3::simd::fill<simd_type>(initial_score);
        } else {
            for (size_t i = 0; i < size; ++i)
                values[i] = initial_score;
        }
    }

    explicit simd_score(simd_type simd_score) noexcept : values{std::move(simd_score)}
    {}

    template <typename ..._score_t>
        requires ((sizeof...(_score_t) == simd_size - 2) && (std::same_as<_score_t, score_t> && ...))
    explicit simd_score(score_t const first_score, score_t const second_score, _score_t const ...remaining_scores)
        noexcept :
        values{{first_score, second_score, remaining_scores...}}
    {}

    template <typename other_score_t>
        requires std::assignable_from<score_t &, other_score_t>
    explicit simd_score(simd_score<other_score_t, simd_size> const & other) noexcept
    {
        for (size_t i = 0; i < size; ++i)
            values[i] = static_cast<score_t>(other[i]);
    }

    constexpr reference operator[](size_t const pos) noexcept
    {
        return values[pos];
    }

    constexpr const_reference operator[](size_t const pos) const noexcept
    {
        return values[pos];
    }

    constexpr simd_score & operator+=(simd_score const & right) noexcept
    {
        values += right.values;
        return *this;
    }

    constexpr simd_score & operator+=(score_t const right_constant) noexcept
    {
        values += simd_score{right_constant}.values;
        return *this;
    }

    constexpr simd_score operator+(simd_score const & right) const noexcept
    {
        simd_score tmp{};
        tmp.values = values + right.values;
        return tmp;
    }

    constexpr simd_score operator+(score_t const right_constant) const noexcept
    {
        simd_score tmp{};
        tmp.values = values + simd_score{right_constant}.values;
        return tmp;
    }

    constexpr simd_score & operator-=(simd_score const & right) noexcept
    {
        values -= right.values;
        return *this;
    }

    constexpr simd_score & operator-=(score_t const right_constant) noexcept
    {
        values -= simd_score{right_constant}.values;
        return *this;
    }

    constexpr simd_score operator-(simd_score const & right) const noexcept
    {
        simd_score tmp{};
        tmp.values = values - right.values;
        return tmp;
    }

    constexpr simd_score operator-(score_t const right_constant) const noexcept
    {
        simd_score tmp{};
        tmp.values = values - simd_score{right_constant}.values;
        return tmp;
    }

    constexpr simd_score & operator*=(simd_score const & right) noexcept
    {
        for (size_t i = 0; i < simd_size; ++i)
            values[i] *= right[i];

        return *this;
    }

    constexpr simd_score & operator*=(score_t const right_constant) noexcept
    {
        return *this *= simd_score{right_constant};
    }

    constexpr simd_score operator*(simd_score const & right) const noexcept
    {
        simd_score tmp{};
        tmp.values = values * right.values;
        return tmp;
    }

    constexpr simd_score operator*(score_t const right_constant) const noexcept
    {
        simd_score tmp{};
        tmp.values = values * simd_score{right_constant}.values;
        return tmp;
    }

    constexpr simd_score operator^(simd_score const & right) const noexcept
    {
        simd_score tmp{};
        tmp.values = values ^ right.values;
        return tmp;
    }

    constexpr friend simd_score max(simd_score const & left, simd_score const & right) noexcept
    {
        simd_score tmp{};
        tmp.values = (left.values < right.values) ? right.values : left.values;
        return tmp;
    }

    template <typename fn_t>
    constexpr friend auto compare(simd_score const & left, simd_score const & right, fn_t && fn) noexcept
        -> std::invoke_result_t<fn_t, simd_score const &, simd_score const &>
    {
        return fn(left, right);
    }

    constexpr auto lt(simd_score const & rhs) const noexcept
    {
        return values < rhs.values;
    }

    constexpr auto le(simd_score const & rhs) const noexcept
    {
        return values <= rhs.values;
    }

    constexpr friend simd_score blend(auto const & mask, simd_score const & left, simd_score const & right) noexcept
    {
        simd_score tmp{};
        tmp.values = mask ? left.values : right.values;
        return tmp;
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
