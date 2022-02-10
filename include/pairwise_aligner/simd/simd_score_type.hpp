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
#include <iomanip>

#include <seqan3/utility/simd/algorithm.hpp>
#include <seqan3/utility/simd/simd_traits.hpp>
#include <seqan3/utility/simd/simd.hpp>

#include <pairwise_aligner/simd/simd_base.hpp>
#include <pairwise_aligner/simd/simd_convert.hpp>
#include <pairwise_aligner/simd/simd_mask_type.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
template <std::integral score_t, size_t simd_size>
class alignas(detail::max_simd_size) simd_score : protected detail::simd_convert_base
{
private:
    using base_t = detail::simd_convert_base;
    using native_simd_t = seqan3::simd::simd_type_t<score_t>;

    static constexpr size_t native_simd_size = seqan3::simd_traits<native_simd_t>::length;
    static constexpr size_t native_simd_count = simd_size / native_simd_size;
    static constexpr bool is_native = native_simd_count == 1;

    template <std::integral, size_t>
    friend class simd_score;

public:

    inline static constexpr size_t size = simd_size;
    inline static constexpr size_t count = native_simd_count;

    using simd_type = std::array<native_simd_t, native_simd_count>;
    using mask_type = simd_mask<std::make_unsigned_t<score_t>, simd_size>;
    using native_simd_type = native_simd_t;
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

    constexpr explicit simd_score(score_t const initial_score) noexcept
    {
        apply([&] (native_simd_t & native_simd_chunk) {
            native_simd_chunk = seqan3::simd::fill<native_simd_t>(initial_score);
        }, values);
    }

    constexpr simd_score(mask_type const mask) noexcept
    {
        apply([&] <typename mask_simd_t> (native_simd_t & native_simd_chunk, mask_simd_t const & mask_simd) {
            native_simd_chunk = static_cast<native_simd_t>(mask_simd);
        }, values, mask.values);
    }

    constexpr explicit simd_score(native_simd_t simd_score) noexcept
        requires (is_native)
    : values{std::move(simd_score)}
    {}

    template <typename ...other_score_t>
        requires ((sizeof...(other_score_t) == simd_size) && sizeof...(other_score_t) > 1 &&
                  (std::convertible_to<score_t, other_score_t> && ...))
    constexpr explicit simd_score(other_score_t const ...values)
        noexcept :
        values{{static_cast<score_t>(values)...}}
    {}

    // cast the other vector in this element
    template <typename other_score_t>
        requires (!std::same_as<other_score_t, score_t> && std::assignable_from<score_t &, other_score_t>)
    constexpr explicit simd_score(simd_score<other_score_t, simd_size> const & other) noexcept
    {
        if constexpr (native_simd_count > other.native_simd_count) { // upcast
            base_t::expand_into(values, other.values[0]);
        } else { // downcast
            base_t::merge_into(values[0], other.values);
        }
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

    constexpr void load(value_type const * mem_address) noexcept
    {
        apply([&] (native_simd_t & native_simd_chunk) {
            native_simd_chunk = seqan3::simd::load<native_simd_t>(mem_address);
            mem_address += native_simd_size; // move the pointer to the next element to load.
        }, values);
    }

    constexpr void store(value_type * mem_address) const noexcept
    {
        apply([&] (native_simd_t const & native_simd_chunk) {
            seqan3::simd::store(mem_address, native_simd_chunk);
            mem_address += native_simd_size; // move the pointer to the next element to store.
        }, values);
    }

    constexpr simd_score & operator++() noexcept
    {
        apply([] (native_simd_t & value) { ++value; }, values);
        return *this;
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
        simd_score tmp{right_constant};
        apply([] (native_simd_t & left, native_simd_t const & right) { left += right; },
              tmp.values, values);
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
        simd_score tmp{right_constant};
        apply([] (native_simd_t & left, native_simd_t const & right) { left = right - left; },
              tmp.values, values);
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
        simd_score tmp{right_constant};
        apply([] (native_simd_t & left, native_simd_t const & right) { left *= right; },
              tmp.values, values);
        return tmp;
    }

    constexpr simd_score operator/(simd_score const & right) const noexcept
    {
        simd_score tmp{};
        apply([] (native_simd_t & res, native_simd_t const & left, native_simd_t const & right) { res = left / right; },
              tmp.values, values, right.values);
        return tmp;
    }

    constexpr simd_score operator/(score_t const right_constant) const noexcept
    {
        simd_score tmp{right_constant};
        apply([] (native_simd_t & left, native_simd_t const & right) { left /= right; },
              tmp.values, values);
        return tmp;
    }

    constexpr simd_score operator^(simd_score const & right) const noexcept
    {
        simd_score tmp{};
        apply([] (native_simd_t & res, native_simd_t const & left, native_simd_t const & right) { res = left ^ right; },
              tmp.values, values, right.values);
        return tmp;
    }

    constexpr simd_score operator^(score_t const right_constant) const noexcept
    {
        simd_score tmp{};
        apply([] (native_simd_t & res, native_simd_t const & left, native_simd_t const & right) { res = left ^ right; },
              tmp.values, values, simd_score{right_constant}.values);
        return tmp;
    }

    constexpr simd_score operator|(simd_score const & right) const noexcept
    {
        simd_score tmp{};
        apply([] (native_simd_t & res, native_simd_t const & left, native_simd_t const & right) { res = left | right; },
              tmp.values, values, right.values);
        return tmp;
    }

    constexpr simd_score & operator|=(simd_score const & right) noexcept
    {
        apply([] (native_simd_t & left, native_simd_t const & right) { left |= right; },
              values, right.values);
        return *this;
    }

    constexpr simd_score operator&(simd_score const & right) const noexcept
    {
        simd_score tmp{};
        apply([] (native_simd_t & res, native_simd_t const & left, native_simd_t const & right) { res = left & right; },
              tmp.values, values, right.values);
        return tmp;
    }

    constexpr simd_score operator&(value_type const right_constant) const noexcept
    {
        simd_score tmp{};
        apply([] (native_simd_t & res, native_simd_t const & left, native_simd_t const & right) { res = left & right; },
              tmp.values, values, simd_score{right_constant}.values);
        return tmp;
    }

    constexpr simd_score operator>>(uint32_t const shift) const noexcept
    {
        simd_score tmp{};
        apply([shift] (native_simd_t & res, native_simd_t const & left) { res = left >> shift; },
              tmp.values, values);
        return tmp;
    }

    // constexpr simd_score operator>>(value_type const right_constant) const noexcept
    // {
    //     simd_score tmp{};
    //     apply([] (native_simd_t & res, native_simd_t const & left, native_simd_t const & right) { res = left >> right; },
    //           tmp.values, values, simd_score{right_constant}.values);
    //     return tmp;
    // }

    constexpr simd_score operator<<(uint32_t const shift) const noexcept
    {
        simd_score tmp{};
        apply([shift] (native_simd_t & res, native_simd_t const & left) { res = left << shift; },
              tmp.values, values);
        return tmp;
    }

    // constexpr simd_score operator<<(value_type const right_constant) const noexcept
    // {
    //     simd_score tmp{};
    //     apply([] (native_simd_t & res, native_simd_t const & left, native_simd_t const & right) { res = left << right; },
    //           tmp.values, values, simd_score{right_constant}.values);
    //     return tmp;
    // }

    constexpr friend simd_score max(simd_score const & left, simd_score const & right) noexcept
    {
        simd_score tmp{};
        apply([] (native_simd_t & res, native_simd_t const & left, native_simd_t const & right) {
                res = (left < right) ? right : left;
        }, tmp.values, left.values, right.values);
        return tmp;
    }

    constexpr friend simd_score min(simd_score const & left, simd_score const & right) noexcept
    {
        simd_score tmp{};
        apply([] (native_simd_t & res, native_simd_t const & left, native_simd_t const & right) {
                res = (left < right) ? left : right;
        }, tmp.values, left.values, right.values);
        return tmp;
    }

    template <typename fn_t>
    constexpr friend auto compare(simd_score const & left, simd_score const & right, fn_t && fn) noexcept
        -> std::invoke_result_t<fn_t, simd_score const &, simd_score const &>
    {
        return fn(left, right);
    }

    constexpr mask_type eq(simd_score const & rhs) const noexcept
    {
        mask_type masks{};
        apply([&] <typename result_t> (result_t & mask, native_simd_t const & left, native_simd_t const & right) {
                mask = left == right;
        }, masks.values, values, rhs.values);
        return masks;
    }

    constexpr mask_type lt(simd_score const & rhs) const noexcept
    {
        mask_type masks{};
        apply([&] <typename result_t> (result_t & mask, native_simd_t const & left, native_simd_t const & right) {
                mask = left < right;
        }, masks.values, values, rhs.values);
        return masks;
    }

    constexpr mask_type le(simd_score const & rhs) const noexcept
    {
        mask_type masks{};
        apply([&] <typename result_t> (result_t & mask, native_simd_t const & left, native_simd_t const & right) {
                mask = left <= right;
        }, masks.values, values, rhs.values);
        return masks;
    }

    constexpr friend simd_score blend(mask_type const & masks, simd_score const & left, simd_score const & right) noexcept
    {
        simd_score tmp{};
        apply([] (native_simd_t & res, auto const & mask, native_simd_t const & left, native_simd_t const & right) {
                res = mask ? left : right;
        }, tmp.values, masks.values, left.values, right.values);
        return tmp;
    }

    constexpr friend simd_score abs(simd_score const & value) noexcept
    {
        simd_score tmp;
        apply([&] (native_simd_t & res, native_simd_t const & a) {
            res = (a < 0) ? -a : a;
        }, tmp.values, value.values);
        return tmp;
    }

    constexpr friend simd_score mask_max(simd_score const & source,
                                         mask_type const & mask,
                                         simd_score const & left,
                                         simd_score const & right) noexcept
    {
        simd_score tmp{};
        apply([] <typename mask_t> (native_simd_t & res, native_simd_t const & src, mask_t const & k,
                                    native_simd_t const & a, native_simd_t const & b) {
                res = k ? (a < b ? b : a) : src;
        }, tmp.values, source.values, mask.values, left.values, right.values);
        return tmp;
    }

    constexpr friend simd_score mask_add(simd_score const & source,
                                         mask_type const & mask,
                                         simd_score const & left,
                                         simd_score const & right) noexcept
    {
        simd_score tmp{};
        apply([] <typename mask_t> (native_simd_t & res, native_simd_t const & src, mask_t const & k,
                                    native_simd_t const & a, native_simd_t const & b) {
                res = k ? a + b : src;
        }, tmp.values, source.values, mask.values, left.values, right.values);
        return tmp;
    }

    constexpr friend simd_score mask_subtract(simd_score const & source,
                                              mask_type const & mask,
                                              simd_score const & left,
                                              simd_score const & right) noexcept
    {
        simd_score tmp{};
        apply([] <typename mask_t> (native_simd_t & res, native_simd_t const & src, mask_t const & k,
                                    native_simd_t const & a, native_simd_t const & b) {
                res = k ? a - b : src;
        }, tmp.values, source.values, mask.values, left.values, right.values);
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

    template <typename target_masks_t, typename source_simd_t, size_t ...idx>
    constexpr void downcast(target_masks_t & target,
                            source_simd_t const & source,
                            std::index_sequence<idx...> const &) const noexcept
    {
      constexpr size_t source_mask_size_v = source.native_simd_size;
      constexpr size_t split_factor = source.native_simd_count / native_simd_count;
      constexpr auto index_seq = std::make_index_sequence<source_mask_size_v>();

      ((downcast_impl<(idx * source_mask_size_v) % native_simd_size>(target[idx / split_factor],
                                                                     source.values[idx],
                                                                     index_seq)), ...);
    }

    template <size_t first_idx, typename other_mask_t, size_t ...idx>
    constexpr void downcast_impl(native_simd_t & target,
                                 other_mask_t const & source,
                                 std::index_sequence<idx...> const &) const noexcept
    {
      ((target[first_idx + idx] = static_cast<score_t>(source[idx])), ...);
    }

    template <typename target_masks_t, typename source_simd_t, size_t ...idx>
    constexpr void upcast(target_masks_t & target,
                          source_simd_t const & source,
                          std::index_sequence<idx...> const &) const noexcept
    {
      constexpr size_t split_factor = native_simd_count / source.native_simd_count;
      constexpr auto index_seq = std::make_index_sequence<native_simd_size>();

      ((upcast_impl<(idx * native_simd_size) % source.native_simd_size>(target[idx],
                                                                        source.values[idx / split_factor],
                                                                        index_seq)), ...);
    }

    template <size_t first_idx, typename other_mask_t, size_t ...idx>
    constexpr  void upcast_impl(native_simd_t & target,
                                other_mask_t const & source,
                                std::index_sequence<idx...> const &) const noexcept
    {
      ((target[idx] = static_cast<score_t>(source[first_idx + idx])), ...);
    }
};

template <typename stream_t, typename simd_score_t>
    requires requires {
        std::remove_cvref_t<simd_score_t>::size;
        typename std::remove_cvref_t<simd_score_t>::simd_type;
    }
inline stream_t & operator<<(stream_t & ostream, simd_score_t && simd_score)
{
    size_t width = 4;
    ostream << "\nidx:";
    for (size_t i = 0; i < simd_score.size - 1; ++i)
        ostream << std::setw(width) << (int32_t) i;
    ostream << std::setw(width) << simd_score.size - 1 << "\nval:";
    for (size_t i = 0; i < simd_score.size - 1; ++i)
        ostream << std::setw(width) << (int32_t) simd_score[i]; // << ", ";

    ostream << std::setw(width) << (int32_t) simd_score[simd_score.size - 1] << "\n";
    return ostream;
}

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
