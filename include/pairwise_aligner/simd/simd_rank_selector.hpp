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

#include <seqan3/utility/type_traits/lazy_conditional.hpp>

#include <pairwise_aligner/simd/simd_base.hpp>
#include <pairwise_aligner/simd/simd_rank_selector_impl_sse4.hpp>
#include <pairwise_aligner/simd/simd_rank_selector_impl_avx2.hpp>
#include <pairwise_aligner/simd/simd_rank_selector_impl_avx512.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

namespace detail {

template <typename key_t>
    requires (key_t::size_v != detail::max_simd_size)
struct simd_rank_selector_default
{
protected:
    using rank_map_t = std::vector<key_t, seqan3::aligned_allocator<key_t, alignof(key_t)>>;

    template <std::ranges::random_access_range ranks_t>
        requires (std::same_as<std::remove_cvref_t<ranks_t>, rank_map_t>)
    static rank_map_t initialise_rank_map(ranks_t && ranks) noexcept
    {
        return ranks;
    }

    static key_t select_rank_for(rank_map_t const & rank_map, key_t const & keys) noexcept
    {
        key_t ranks{};
        for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(key_t::size_v); ++i) {
            std::ptrdiff_t index = keys[i] / key_t::size_v;
            std::ptrdiff_t position = keys[i] % key_t::size_v;
            ranks[i] = rank_map[index][position];
        }
        return ranks;
    }
};

template <typename index_t>
using eight_bit_rank_selector_t = seqan3::detail::lazy_conditional_t<
                                    detail::max_simd_size == 64,
                                    seqan3::detail::lazy<simd_rank_selector_impl_avx512, index_t>,
                                    seqan3::detail::lazy_conditional_t<
                                        detail::max_simd_size == 32,
                                        seqan3::detail::lazy<simd_rank_selector_impl_avx2, index_t>,
                                        seqan3::detail::lazy<simd_rank_selector_impl_sse4, index_t>
                                    >
                                >;

template <typename index_t>
using simd_rank_selector_t =
        seqan3::detail::lazy_conditional_t<index_t::size_v == detail::max_simd_size,
                            seqan3::detail::lazy<eight_bit_rank_selector_t, index_t>,
                            seqan3::detail::lazy<simd_rank_selector_default, index_t>
                        >;

} // namespace detail

} // v1
} // namespace seqan::pairwise_aligner
