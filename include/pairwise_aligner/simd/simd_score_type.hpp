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

template <std::integral score_t, size_t simd_size = detail::max_simd_size / sizeof(score_t)>
class alignas(detail::max_simd_size) simd_score
{
private:
    static constexpr size_t native_simd_size = detail::max_simd_size / sizeof(score_t);
    static constexpr size_t native_simd_count = simd_size / native_simd_size;
    static constexpr bool is_native = native_simd_count == 1;
    using native_simd_t = seqan3::simd::simd_type_t<score_t, native_simd_size>;

public:

    inline static constexpr size_t size = simd_size;

    using simd_type = std::array<native_simd_t, native_simd_count>;
    using value_type = score_t;
    using reference = score_t &;
    using const_reference = score_t const &;

private:

    simd_type values{};

public:

    simd_score() = default;
    simd_score(simd_score const &) = default;
    simd_score(simd_score &&) = default;
    simd_score & operator=(simd_score const &) = default;
    simd_score & operator=(simd_score &&) = default;

    explicit simd_score(score_t const initial_score) noexcept
    {
        apply([&] (native_simd_t & native_simd_chunk) {
            native_simd_chunk = seqan3::simd::fill<native_simd_t>(initial_score);
        }, values);
    }

    explicit simd_score(native_simd_t simd_score) noexcept
        requires (is_native)
    : values{std::move(simd_score)}
    {}

    template <typename ..._score_t>
        requires ((sizeof...(_score_t) == simd_size - 2) && (std::same_as<_score_t, score_t> && ...))
    explicit simd_score(score_t const first_score, score_t const second_score, _score_t const ...remaining_scores)
        noexcept :
        values{{first_score, second_score, remaining_scores...}}
    {}

    // cast the other vector in this element
    template <typename other_score_t>
        requires std::assignable_from<score_t &, other_score_t>
    explicit simd_score(simd_score<other_score_t, simd_size> const & other) noexcept
    {
        for (size_t j = 0; j < native_simd_count; ++j)
            for (size_t i = 0; i < native_simd_size; ++i)
                values[j][i] = static_cast<score_t>(other[i + (j * native_simd_size)]);
    }

    constexpr reference operator[](size_t const pos) noexcept
    {
        auto [index, offset] = to_local_position(pos);
        return values[index][offset];
    }

    constexpr const_reference operator[](size_t const pos) const noexcept
    {
        auto [index, offset] = to_local_position(pos);
        return values[index][offset];
    }

    constexpr simd_score & operator+=(simd_score const & right) noexcept
    {
        apply([] (native_simd_t & left, native_simd_t const & right) { left += right; },
              values, right.values);
        return *this;
    }

    constexpr simd_score & operator+=(score_t const right_constant) noexcept
    {
        apply([] (native_simd_t & left, native_simd_t const & right) { left += right; },
              values, simd_score{right_constant}.values);
        return *this;
    }

    constexpr simd_score operator+(simd_score const & right) const noexcept
    {
        simd_score tmp{};
        apply([] (native_simd_t & res, native_simd_t const & left, native_simd_t const & right) { res = left + right; },
              tmp.values, values, right.values);
        return tmp;
    }

    constexpr simd_score operator+(score_t const right_constant) const noexcept
    {
        simd_score tmp{};
        apply([] (native_simd_t & res, native_simd_t const & left, native_simd_t const & right) { res = left + right; },
              tmp.values, values, simd_score{right_constant}.values);
        return tmp;
    }

    constexpr simd_score & operator-=(simd_score const & right) noexcept
    {
        apply([] (native_simd_t & left, native_simd_t const & right) { left -= right; },
              values, right.values);
        return *this;
    }

    constexpr simd_score & operator-=(score_t const right_constant) noexcept
    {
        apply([] (native_simd_t & left, native_simd_t const & right) { left -= right; },
              values, simd_score{right_constant}.values);
        return *this;
    }

    constexpr simd_score operator-(simd_score const & right) const noexcept
    {
        simd_score tmp{};
        apply([] (native_simd_t & res, native_simd_t const & left, native_simd_t const & right) { res = left - right; },
             tmp.values, values, right.values);
        return tmp;
    }

    constexpr simd_score operator-(score_t const right_constant) const noexcept
    {
        simd_score tmp{};
        apply([] (native_simd_t & res, native_simd_t const & left, native_simd_t const & right) { res = left - right; },
              tmp.values, values, simd_score{right_constant}.values);
        return tmp;
    }

    constexpr simd_score & operator*=(simd_score const & right) noexcept
    {
        apply([] (native_simd_t & left, native_simd_t const & right) { left *= right; },
              values, right.values);
        return *this;
    }

    constexpr simd_score & operator*=(score_t const right_constant) noexcept
    {
        apply([] (native_simd_t & left, native_simd_t const & right) { left *= right; },
              values, simd_score{right_constant}.values);
        return *this;
    }

    constexpr simd_score operator*(simd_score const & right) const noexcept
    {
        simd_score tmp{};
        apply([] (native_simd_t & res, native_simd_t const & left, native_simd_t const & right) { res = left * right; },
              tmp.values, values, right.values);
        return tmp;
    }

    constexpr simd_score operator*(score_t const right_constant) const noexcept
    {
        simd_score tmp{};
        apply([] (native_simd_t & res, native_simd_t const & left, native_simd_t const & right) { res = left * right; },
              tmp.values, values, simd_score{right_constant}.values);
        return tmp;
    }

    constexpr simd_score operator^(simd_score const & right) const noexcept
    {
        simd_score tmp{};
        apply([] (native_simd_t & res, native_simd_t const & left, native_simd_t const & right) { res = left ^ right; },
              tmp.values, values, right.values);
        return tmp;
    }

    constexpr friend simd_score max(simd_score const & left, simd_score const & right) noexcept
    {
        simd_score tmp{};
        apply([] (native_simd_t & res, native_simd_t const & left, native_simd_t const & right) {
                res = (left < right) ? right : left;
        }, tmp.values, left.values, right.values);
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
        using result_t = decltype(std::declval<native_simd_t>() < std::declval<native_simd_t>());
        std::array<result_t, native_simd_count> masks{};
        apply([&] (result_t & mask, native_simd_t const & left, native_simd_t const & right) {
                mask = left < right;
        }, masks, values, rhs.values);
        return masks;
    }

    constexpr auto le(simd_score const & rhs) const noexcept
    {
        using result_t = decltype(std::declval<native_simd_t>() <= std::declval<native_simd_t>());
        std::array<result_t, native_simd_count> masks{};
        apply([&] (result_t & mask, native_simd_t const & left, native_simd_t const & right) {
                mask = left <= right;
        }, masks, values, rhs.values);
        return masks;
    }

    constexpr friend simd_score blend(auto const & masks, simd_score const & left, simd_score const & right) noexcept
    {
        simd_score tmp{};
        apply([] (native_simd_t & res, auto const & mask, native_simd_t const & left, native_simd_t const & right) {
                res = mask ? left : right;
        }, tmp.values, masks, left.values, right.values);
        return tmp;
    }

private:

    constexpr auto to_local_position(size_t const position) const noexcept
    {
        return std::pair<size_t, size_t>{position / native_simd_size, position % native_simd_size};
    }

    template <typename fn_t, typename first_simd_vector_t, typename ...remaining_simd_vector_t>
    static constexpr void apply(fn_t && fn, first_simd_vector_t && first, remaining_simd_vector_t && ...remaining) noexcept
    {
        for (size_t i = 0; i < native_simd_count; ++i)
            fn(first[i], remaining[i]...);
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
