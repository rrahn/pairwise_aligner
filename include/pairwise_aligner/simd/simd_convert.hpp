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

#include <pairwise_aligner/simd/simd_base.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace detail {

class simd_convert_base
{
public:

    template <typename target_simd_t, typename source_simd_t, size_t source_count>
    constexpr void merge_into(target_simd_t & target,
                              std::array<source_simd_t, source_count> const & source_array) const noexcept
    {
        merge_into_impl(target, source_array, std::make_index_sequence<seqan3::simd_traits<target_simd_t>::length>());
    }

    template <typename target_simd_t, size_t target_count, typename source_simd_t>
    constexpr void expand_into(std::array<target_simd_t, target_count> & target_array,
                               source_simd_t const & source) const noexcept
    {
        expand_into_impl(target_array, source, std::make_index_sequence<seqan3::simd_traits<source_simd_t>::length>());
    }

    private:

    template <typename target_simd_t, typename source_simd_t, size_t source_count, size_t ...vec_idx>
    constexpr void merge_into_impl(target_simd_t & target,
                                    std::array<source_simd_t, source_count> const & source_array,
                                    std::index_sequence<vec_idx...> const &) const noexcept
    {
        using scalar_t = typename seqan3::simd_traits<target_simd_t>::scalar_type;
        constexpr size_t simd_size_v = seqan3::simd_traits<source_simd_t>::length;
        target = target_simd_t{
            (static_cast<scalar_t>(std::get<vec_idx / simd_size_v>(source_array)[vec_idx % simd_size_v]))...
        };
    }

    template <typename target_simd_t, size_t target_count, typename source_simd_t, size_t ...idx>
    constexpr void expand_into_impl(std::array<target_simd_t, target_count> & target_array,
                                    source_simd_t const & source,
                                    std::index_sequence<idx...> const &) const noexcept
    {
        using scalar_t = typename seqan3::simd_traits<target_simd_t>::scalar_type;
        constexpr size_t simd_size_v = seqan3::simd_traits<target_simd_t>::length;

        ((target_array[idx / simd_size_v][idx % simd_size_v] = static_cast<scalar_t>(source[idx])), ...);
    }
};

} // namespace detail
} // v1
} // namespace seqan::pairwise_aligner
