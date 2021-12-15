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

#include <seqan3/utility/simd/all.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <std::integral score_t, size_t simd_size = 32/sizeof(score_t)>
class simd_score
{
private:

    static constexpr bool use_intrinsics = sizeof(score_t) <= 1;
    using simd_type = std::conditional_t<use_intrinsics,
                                        seqan3::simd_type_t<score_t, simd_size>,
                                        std::array<score_t, simd_size>>;

public:

    simd_type values{};

    inline static constexpr size_t size = simd_size;

    using value_type = score_t;
    using reference = score_t &;
    using const_reference = score_t const &;

    simd_score() = default;
    simd_score(simd_score const &) = default;
    simd_score(simd_score &&) = default;
    simd_score & operator=(simd_score const &) = default;
    simd_score & operator=(simd_score &&) = default;

    explicit simd_score(score_t const initial_score) noexcept : values{}
    {
        if constexpr (use_intrinsics) {
            values = seqan3::simd::fill<simd_type>(initial_score);
        } else {
            values.fill(initial_score);
        }
    }

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
        // if constexpr (use_intrinsics) {
        //     values = seqan3::simd::upcast<simd_type>(other.values);
        // } else {
        for (size_t i = 0; i < size; ++i)
            values[i] = static_cast<score_t>(other.values[i]);
        // }
    }

    reference operator[](size_t const pos) noexcept
    {
        return values[pos];
    }

    const_reference operator[](size_t const pos) const noexcept
    {
        return values[pos];
    }

    simd_score & operator+=(simd_score const & right) noexcept
    {
        if constexpr (use_intrinsics) {
            values += right.values;
        } else {
            simd_score tmp = *this + right;
            swap(tmp);
        }

        return *this;
    }

    simd_score & operator+=(score_t const right_constant) noexcept
    {
        if constexpr (use_intrinsics) {
            values += simd_type{right_constant};
        } else {
            simd_score tmp = *this + right_constant;
            swap(tmp);
        }

        return *this;
    }

    simd_score operator+(simd_score const & right) const noexcept
    {
        simd_score tmp{};
        if constexpr (use_intrinsics) {
            tmp.values = values + right.values;
        } else {
            for (size_t i = 0; i < simd_size; ++i)
                tmp[i] = values[i] + right[i];
        }

        return tmp;
    }

    simd_score operator+(score_t const right_constant) const noexcept
    {
        simd_score tmp{};
        if constexpr (use_intrinsics) {
            tmp.values = values + simd_type{right_constant};
        } else {
            for (size_t i = 0; i < simd_size; ++i)
                tmp[i] = values[i] + right_constant;
        }

        return tmp;
    }

    simd_score & operator-=(simd_score const & right) noexcept
    {
        if constexpr (use_intrinsics) {
            values -= right.values;
        } else {
            simd_score tmp = *this - right;
            swap(tmp);
        }

        return *this;
    }

    simd_score & operator-=(score_t const right_constant) noexcept
    {
        if constexpr (use_intrinsics) {
            values -= simd_type{right_constant};
        } else {
            simd_score tmp = *this - right_constant;
            swap(tmp);
        }

        return *this;
    }

    simd_score operator-(simd_score const & right) const noexcept
    {
        simd_score tmp{};
        if constexpr (use_intrinsics) {
            tmp.values = values - right.values;
        } else {
            for (size_t i = 0; i < simd_size; ++i)
                tmp[i] = values[i] - right[i];
        }

        return tmp;
    }

    simd_score operator-(score_t const right_constant) const noexcept
    {
        simd_score tmp{};
        if constexpr (use_intrinsics) {
            tmp.values = values - simd_type{right_constant};
        } else {
            for (size_t i = 0; i < simd_size; ++i)
                tmp[i] = values[i] - right_constant;
        }

        return tmp;
    }

    simd_score & operator*=(simd_score const & right) noexcept
    {
        if constexpr (use_intrinsics) {
            values *= right.values;
        } else {
            simd_score tmp = *this * right;
            swap(tmp);
        }

        return *this;
    }

    simd_score & operator*=(score_t const right_constant) noexcept
    {
        if constexpr (use_intrinsics) {
            values *= simd_type{right_constant};
        } else {
            simd_score tmp = *this * right_constant;
            swap(tmp);
        }

        return *this;
    }

    simd_score operator*(simd_score const & right) const noexcept
    {
        simd_score tmp{};
        if constexpr (use_intrinsics) {
            tmp.values = values * right.values;
        } else {
            for (size_t i = 0; i < simd_size; ++i)
                tmp[i] = values[i] * right[i];
        }

        return tmp;
    }

    simd_score operator*(score_t const right_constant) const noexcept
    {
        simd_score tmp{};
        if constexpr (use_intrinsics) {
            tmp.values = values * simd_type{right_constant};
        } else {
            for (size_t i = 0; i < simd_size; ++i)
                tmp[i] = values[i] * right_constant;
        }

        return tmp;
    }

    friend simd_score max(simd_score const & left, simd_score const & right) noexcept
    {
        simd_score tmp{};
        if constexpr (use_intrinsics) {
            tmp.values = left.values < right.values ? right.values : left.values;
        } else {
            for (size_t i = 0; i < simd_size; ++i)
                tmp[i] = std::max(left[i], right[i]);
        }

        return tmp;
    }

    template <typename fn_t>
    friend simd_score compare(simd_score const left, simd_score const right, fn_t && fn) noexcept
    {
        simd_score tmp{};
        if constexpr (use_intrinsics) {
            tmp.values = fn(left.values, right.values);
        } else {
            for (size_t i = 0; i < simd_size; ++i)
                tmp[i] = fn(left[i], right[i]);
        }

        return tmp;
    }

    friend simd_score blend(simd_score const mask, simd_score const left, simd_score const right) noexcept
    {
        simd_score tmp{};
        if constexpr (use_intrinsics) {
            tmp.values = mask.values ? left.values : right.values;
        } else {
            for (size_t i = 0; i < simd_size; ++i)
                tmp[i] = mask[i] ? left[i] : right[i];
        }

        return tmp;
    }

    friend simd_score compare_and_blend(simd_score const & left,
                                        simd_score const & right,
                                        score_t const true_value,
                                        score_t const false_value) noexcept
    {
        simd_score tmp{};
        for (size_t i = 0; i < simd_size; ++i)
            tmp[i] = left[i] == right[i] ? true_value : false_value;

        return tmp;
    }

    friend simd_score compare_and_blend(simd_score const & left,
                                        simd_score const & right,
                                        simd_score const true_value,
                                        simd_score const false_value) noexcept
    {
        simd_score tmp{};
        for (size_t i = 0; i < simd_size; ++i)
            tmp[i] = left[i] == right[i] ? true_value[i] : false_value[i];

        return tmp;
    }

    void swap(simd_score & right) noexcept
    {
        values.swap(right.values);
    }
private:

    // template <typename fn_t>
    // static void apply(fn_t && fn, simd_score const & left, simd_score const & right)
    // {
    //     for (size_t i = 0; i < simd_size; ++i)
    //         fn(left[i], right[i]);
    // }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
