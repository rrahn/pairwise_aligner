// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides pack algorithms.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <seqan3/utility/type_pack/traits.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename predicate_t, typename ...args_t>
using find_if_t = seqan3::pack_traits::at<seqan3::pack_traits::find_if<predicate_t, args_t...>, args_t...>;

} // inline namespace v1
} // namespace seqan::pairwise_aligner
