// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::end_gap, seqan::pairwise_aligner::leading_end_gap, and
 *        seqan::pairwise_aligner::trailing_end_gap.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

namespace seqan::pairwise_aligner {
inline namespace v1 {
namespace cfg {

enum struct end_gap : bool
{
    penalised,
    free
};

struct leading_end_gap
{
    end_gap first_column{end_gap::penalised};
    end_gap first_row{end_gap::penalised};
};

struct trailing_end_gap
{
    end_gap last_column{end_gap::penalised};
    end_gap last_row{end_gap::penalised};
};

} // namespace cfg
} // inline namespace v1
} // namespace seqan::pairwise_aligner
