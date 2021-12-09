// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_trailing_gaps.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

namespace seqan::pairwise_aligner
{
inline namespace v1
{
enum struct dp_trailing_gaps
{
    regular,
    free
};

struct trailing_gap_setting
{
    dp_trailing_gaps column{dp_trailing_gaps::regular};
    dp_trailing_gaps row{dp_trailing_gaps::regular};
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
