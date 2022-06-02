// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides score scout.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

namespace seqan::align
{
inline namespace v1
{

template <typename score_matrix_t, bool last_row, bool last_column>
struct global_score_scout {

    score_matrix_t _score_matrix;

    score_scout(score_matrix_t score_matrix) : _score_matrix{std::forward<score_matrix_t>(score_matrix_t)}
    {}

private:

    template <typename CPO, typename ...args_t>
    constexpr friend void tag_invoke(CPO cpo, score_scout & me, args_t && ...args) noexcept {
        cpo(me._score_matrix, std::forward<args_t>(args)...);
    }

    // its sole task is to track the optimal score
    // this depends a bit on the way how the score is stored.
    // so to make this independent the sore matrix offers the interface for
        // last column, last row
    // must connect to the score_matrix
    // is set by the aligner that we use.
};
} // inline namespace v1
} // namespace seqan::align
