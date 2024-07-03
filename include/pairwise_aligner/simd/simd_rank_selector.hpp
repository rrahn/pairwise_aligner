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

template <typename index_t>
using simd_rank_selector_t = seqan3::detail::lazy_conditional_t<
    detail::max_simd_size == 64,
    seqan3::detail::lazy<simd_rank_selector_impl_avx512, index_t>,
    seqan3::detail::lazy_conditional_t<
        detail::max_simd_size == 32,
        seqan3::detail::lazy<simd_rank_selector_impl_avx2, index_t>,
        seqan3::detail::lazy<simd_rank_selector_impl_sse4, index_t>>>;

} // namespace detail

} // v1
} // namespace seqan::pairwise_aligner
