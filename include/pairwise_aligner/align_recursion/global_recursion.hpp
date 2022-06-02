// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides initialisation switch.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once


namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename row_initialiser_t, typename column_initialiser_t>
struct dp_matrix_initialiser {
    row_initialiser_t row_initialiser;
    column_initialiser_t column_initialiser;
};

} // inline namespace v1
} // namespace seqan::pairwise_aligner
