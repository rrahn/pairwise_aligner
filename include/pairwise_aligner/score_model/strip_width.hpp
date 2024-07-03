// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::interleaved_score_profile.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <seqan3/std/type_traits>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <size_t width_v>
using strip_width_t = std::integral_constant<size_t, width_v>;

template <size_t width_v>
inline constexpr strip_width_t<width_v> strip_width;

} // inline namespace v1
} // namespace seqan::pairwise_aligner
