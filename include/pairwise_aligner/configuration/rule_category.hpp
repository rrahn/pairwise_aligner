// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::cfg::detail::rule_category.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <cstdint>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace cfg::detail
{

enum struct rule_category : uint8_t
{
    score_model = 0,
    gap_model = 1,
    size = 2
};

} // namespace cfg::detail
} // inline namespace v1
} // namespace seqan::pairwise_aligner
