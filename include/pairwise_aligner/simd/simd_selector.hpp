// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise::simd_selector_avx512.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <array>
#include <seqan3/std/concepts>
#include <seqan3/std/type_traits>

#include <pairwise_aligner/simd/concept.hpp>
#include <pairwise_aligner/simd/simd_base.hpp>
// #include <pairwise_aligner/simd/simd_selector_impl_sse4.hpp>
#include <pairwise_aligner/simd/simd_selector_impl_avx2.hpp>
#include <pairwise_aligner/simd/simd_selector_impl_avx512.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

namespace detail {
// ----------------------------------------------------------------------------
// Default
// ----------------------------------------------------------------------------
template <typename simd_offset_t, size_t element_count, size_t simd_lane_width, size_t simd_register_width>
struct simd_selector<simd_offset_t, selector_tag<element_count, simd_lane_width, simd_register_width>>
{
    // Note: element count can not be larger than maximal value range of operand type.
    static constexpr bool in_lane_shuffle = false;
    static constexpr size_t max_operand_count = element_count;

    template <typename simd_t> // expected to be simd value type
    using address_t = std::array<std::array<typename simd_t::value_type, element_count>, 1>;

    simd_offset_t const & offsets;

    template <typename value_t>
    constexpr simd_offset_t operator()(std::array<std::array<value_t, element_count>, 1> const & data) const noexcept {
        simd_offset_t tmp{};
        for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(simd_offset_t::size_v); ++i) {
            tmp[i] = data[0][offsets[i]];
        }
        return tmp;
    }
};

} // namespace detail
template <typename simd_value_t, typename simd_offset_t, size_t size_v>
struct _simd_selector
{
    class type;
};

template <simd::simd_type simd_value_t, simd::simd_type simd_offset_t, size_t size_v>
using simd_selector = typename _simd_selector<simd_value_t, simd_offset_t, size_v>::type;

template <typename simd_value_t, typename simd_offset_t, size_t size_v>
class _simd_selector<simd_value_t, simd_offset_t, size_v>::type
{
private:

    using offset_t = typename simd_offset_t::value_type;

    static_assert(size_v - 1 <= std::numeric_limits<offset_t>::max(),
                  "The offset type is too small to access all elements in the range.");

    using selector_impl_t = detail::simd_selector<simd_offset_t,
                                                  detail::selector_tag<size_v, // what is size_v: 231 -> matrix size
                                                                       sizeof(offset_t) * 8, // int32_t -> 32bit
                                                                       detail::max_simd_size * 8>>; // 256 bit

public:

    using address_t = typename selector_impl_t::template address_t<simd_value_t>;
    static constexpr std::ptrdiff_t elements_per_select = selector_impl_t::max_operand_count;

    template <std::ranges::forward_range data_slice_t>
    static constexpr address_t load(data_slice_t && data_slice) noexcept
    {
        assert(std::ranges::distance(data_slice) <= elements_per_select);

        using address_value_t = typename address_t::value_type;

        address_t address{};
        for (size_t idx = 0; idx < address.size(); ++idx) {
            address_value_t & values = address[idx];

            auto data_it = std::ranges::next(std::ranges::begin(data_slice),
                                             idx * values.size(),
                                             std::ranges::end(data_slice));
            auto data_end = std::ranges::next(std::ranges::begin(data_slice),
                                              (idx + 1) * values.size(),
                                              std::ranges::end(data_slice));

            for (size_t data_idx = 0; data_it != data_end; ++data_it, ++data_idx)
                values[data_idx] = static_cast<typename address_value_t::value_type>(*data_it);

            // Handle in lane shuffle by copying the data of the first lane in all following lanes.
            if constexpr (selector_impl_t::in_lane_shuffle) {
                for (size_t simd_idx = selector_impl_t::max_operand_count; simd_idx < simd_value_t::size_v; ++simd_idx) {
                    values[simd_idx] = values[simd_idx % selector_impl_t::max_operand_count];
                }
            }
        }

        return address;
    }

    static constexpr selector_impl_t selector_for(simd_offset_t const & offset) noexcept
    {
        return selector_impl_t{offset};
    }
};

} // v1
} // namespace seqan::pairwise_aligner
