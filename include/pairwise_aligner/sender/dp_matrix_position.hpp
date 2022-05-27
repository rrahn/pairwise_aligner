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

namespace seqan::align
{
inline namespace v1
{

struct dp_matrix_position {
    size_t row_index;
    size_t column_index;
};

} // inline namespace v1
} // namespace seqan::align
