// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise::simd_mask.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <array>
#include <seqan3/std/concepts>

#include <seqan3/utility/simd/algorithm.hpp>
#include <seqan3/utility/simd/simd_traits.hpp>
#include <seqan3/utility/simd/simd.hpp>

#include <pairwise_aligner/simd/simd_base.hpp>
#include <pairwise_aligner/simd/simd_convert.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <std::unsigned_integral score_t, size_t simd_size>
class alignas(detail::max_simd_size) simd_mask : protected detail::simd_convert_base
{
private:
    using base_t = detail::simd_convert_base;
    using native_simd_t = seqan3::simd::simd_type_t<score_t>;
    using native_mask_t = typename seqan3::simd_traits<native_simd_t>::mask_type;

    static constexpr size_t native_simd_size = seqan3::simd_traits<native_simd_t>::length;
    static constexpr size_t native_simd_count = simd_size / native_simd_size;

    template <typename, size_t, template <typename> typename ...>
    friend class simd_score_base;

    template <std::unsigned_integral, size_t>
    friend class simd_mask;

public:

    inline static constexpr size_t size = simd_size;

    using mask_type = std::array<native_mask_t, native_simd_count>;
    using value_type = bool;
    using reference = bool &;
    using const_reference = bool;

    mask_type values{};

public:

    simd_mask() = default;
    simd_mask(simd_mask const &) = default;
    simd_mask(simd_mask &&) = default;
    simd_mask & operator=(simd_mask const &) = default;
    simd_mask & operator=(simd_mask &&) = default;
    ~simd_mask() = default;

    constexpr explicit simd_mask(value_type const initial_value) noexcept
    {
        apply([&] (native_mask_t & native_mask_chunk) {
                native_mask_chunk = (initial_value ? ~native_mask_t{} : native_mask_t{});
        }, values);
    }

    // cast the other vector in this element
    template <typename other_score_t>
        requires (!std::same_as<other_score_t, score_t> && std::assignable_from<score_t &, other_score_t>)
    constexpr explicit simd_mask(simd_mask<other_score_t, simd_size> const & other) noexcept
    {
        if constexpr (native_simd_count < other.native_simd_count) { // downcast: merge
            base_t::merge_into(values[0], other.values);
        } else {  // upcast: expand
            base_t::expand_into(values, other.values[0]);
        }
    }

    constexpr const_reference operator[](size_t const pos) const noexcept
    {
        auto [index, offset] = to_local_position(pos);
        return values[index][offset];
    }

    constexpr simd_mask operator&&(simd_mask tmp) const noexcept
    {
        apply([] (native_mask_t & left, native_mask_t const & right) { left = left && right; },
              tmp.values, values);
        return tmp;
    }

    constexpr simd_mask operator&(simd_mask const & right) const noexcept
    {
        simd_mask tmp{};
        apply([] (native_mask_t & res, native_mask_t const & left, native_mask_t const & right) { res = left & right; },
              tmp.values, values, right.values);
        return tmp;
    }

    constexpr simd_mask & operator&=(simd_mask const & right) noexcept
    {
        apply([] (native_mask_t & left, native_mask_t const & right) { left &= right; },
              values, right.values);
        return *this;
    }

    constexpr simd_mask operator|(simd_mask const & right) const noexcept
    {
        simd_mask tmp{};
        apply([] (native_mask_t & res, native_mask_t const & left, native_mask_t const & right) { res = left | right; },
              tmp.values, values, right.values);
        return tmp;
    }

    constexpr simd_mask & operator|=(simd_mask const & right) noexcept
    {
        apply([] (native_mask_t & left, native_mask_t const & right) { left |= right; },
              values, right.values);
        return *this;
    }

    constexpr simd_mask operator~() const noexcept
    {
        simd_mask tmp{};
        apply([] (native_mask_t & left, native_mask_t const & right) { left = ~right; },
              tmp.values, values);
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
      constexpr size_t source_mask_size_v = source.native_simd_size; // [0]
      constexpr size_t split_factor = source.native_simd_count / native_simd_count;
      constexpr auto index_seq = std::make_index_sequence<source_mask_size_v>();

      ((downcast_impl<(idx * source_mask_size_v) % native_simd_size>(target[idx / split_factor],
                                                                     source.values[idx],
                                                                     index_seq)), ...);
    }

    template <size_t first_idx, typename other_mask_t, size_t ...idx>
    constexpr void downcast_impl(native_mask_t & target,
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
    constexpr  void upcast_impl(native_mask_t & target,
                                other_mask_t const & source,
                                std::index_sequence<idx...> const &) const noexcept
    {
      ((target[idx] = static_cast<score_t>(source[first_idx + idx])), ...);
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
