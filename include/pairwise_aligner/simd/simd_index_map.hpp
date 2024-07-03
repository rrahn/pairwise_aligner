// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise::simd_index_map.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <seqan3/std/concepts>
#include <seqan3/std/ranges>
#include <seqan3/std/type_traits>

#include <seqan3/utility/type_traits/lazy_conditional.hpp>
#include <seqan3/utility/views/slice.hpp>

#include <pairwise_aligner/simd/simd_base.hpp>
// #include <pairwise_aligner/simd/simd_selector_sse4.hpp>
// #include <pairwise_aligner/simd/simd_selector_avx2.hpp>
#include <pairwise_aligner/simd/simd_selector.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <std::integral value_t, std::integral key_t, std::ptrdiff_t size_v>
    requires (sizeof(value_t) == sizeof(key_t))
struct _simd_index_map
{
    class type;
};

template <std::integral value_t, std::integral key_t, std::ptrdiff_t size_v>
    requires (sizeof(value_t) == sizeof(key_t))
using simd_index_map = typename _simd_index_map<value_t, key_t, size_v>::type;

template <std::integral value_t, std::integral key_t, std::ptrdiff_t size_v>
    requires (sizeof(value_t) == sizeof(key_t))
class _simd_index_map<value_t, key_t, size_v>::type
{
    using simd_value_t = simd_score<value_t>;
    using simd_key_t = simd_score<key_t>;

    using select_strategy_t = simd_selector<simd_value_t, simd_key_t, size_v>;

    static constexpr std::ptrdiff_t elements_per_select = select_strategy_t::elements_per_select;
    static constexpr std::ptrdiff_t select_count = (size_v - 1 + elements_per_select) / elements_per_select;

    using address_t = typename select_strategy_t::address_t;
    using container_t = std::array<address_t, select_count>;

    // The container holding the data.
    container_t _data{};

public:
    using value_type = simd_value_t;
    using key_type = simd_key_t;

    type() = default;
    template <std::ranges::forward_range data_t>
    explicit constexpr type(data_t && data) noexcept
    {
        assert(std::ranges::distance(data) == size_v);

        std::ranges::for_each(std::views::iota(0, select_count), [&] (std::ptrdiff_t const cycle) {
            _data[cycle] = select_strategy_t::load(seqan3::views::slice(data,
                                                                        cycle * elements_per_select,
                                                                        (cycle + 1) * elements_per_select));
        });
    }

    constexpr value_type operator[](key_type const & key) const noexcept
    {
        auto selector = select_strategy_t::selector_for(key);

        if constexpr (select_count == 1) {
            return to_value(selector(_data[0]));
        } else if constexpr (select_count == 2) {
            auto a = selector(_data[0]);
            auto b = selector(_data[1]);

            return blend(key.lt(key_type{static_cast<key_t>(elements_per_select)}), to_value(a), to_value(b));
        } else {
            value_type result{};
            for (std::ptrdiff_t cycle = 0; cycle < select_count; ++cycle) {
                result = blend(key_type{static_cast<key_t>(cycle * elements_per_select)}.le(key),
                               to_value(selector(_data[cycle])),
                               result);
            }
            return result;
        }
    }

private:

    template <typename native_simd_t>
    value_type to_value(native_simd_t && native) const noexcept
    {
        return reinterpret_cast<value_type &&>(native);
    }

    template <typename native_simd_t>
    value_type const & to_value(native_simd_t const & native) const noexcept
    {
        return reinterpret_cast<value_type const &>(native);
    }

};

} // v1
} // namespace seqan::pairwise_aligner
