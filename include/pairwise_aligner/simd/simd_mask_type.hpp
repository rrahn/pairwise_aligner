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

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <std::unsigned_integral score_t, size_t simd_size>
class alignas(detail::max_simd_size) simd_mask
{
private:
    using native_simd_t = seqan3::simd::simd_type_t<score_t>;
    using native_mask_t = typename seqan3::simd_traits<native_simd_t>::mask_type;

    static constexpr size_t native_simd_size = seqan3::simd_traits<native_simd_t>::length;
    static constexpr size_t native_simd_count = simd_size / native_simd_size;

    template <std::integral, size_t>
    friend class simd_score;

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

// template <typename stream_t, typename simd_mask_t>
//     requires requires {
//         std::remove_cvref_t<simd_mask_t>::size;
//         typename std::remove_cvref_t<simd_mask_t>::simd_type;
//     }
// inline stream_t & operator<<(stream_t & ostream, simd_mask_t && simd_mask)
// {
//     ostream << "<";
//     for (size_t i = 0; i < simd_mask.size - 1; ++i)
//         ostream << (int32_t) simd_mask[i] << ", ";

//     ostream << (int32_t) simd_mask[simd_mask.size - 1] << ">";
//     return ostream;
// }

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
