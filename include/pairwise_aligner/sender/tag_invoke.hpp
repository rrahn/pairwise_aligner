// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides .
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <unifex/tag_invoke.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

using unifex::tag_t;
using unifex::tag_invoke_result_t;
using unifex::tag_invoke;

} // inline namespace v1
} // namespace seqan::pairwise_aligner
