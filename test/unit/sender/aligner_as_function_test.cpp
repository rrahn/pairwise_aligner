// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string_view>

#include <pairwise_aligner/score_model/score_model_unitary.hpp>
#include <pairwise_aligner/configuration/configure_aligner.hpp>
#include <pairwise_aligner/configuration/gap_model_affine.hpp>
#include <pairwise_aligner/configuration/method_global.hpp>
#include <pairwise_aligner/configuration/score_model_unitary.hpp>

#include <unifex/just.hpp>
#include <unifex/sync_wait.hpp>
#include <unifex/then.hpp>

TEST(basic_sender_test, call_configured_aligner)
{
    namespace pa = seqan::pairwise_aligner;

    auto aligner = pa::cfg::configure_aligner(
        pa::cfg::method_global(
            pa::cfg::gap_model_affine(
                pa::cfg::score_model_unitary(4, -5),
                -10, -1
            ),
            pa::cfg::leading_end_gap{.first_column = pa::cfg::end_gap::penalised, .first_row = pa::cfg::end_gap::free},
            pa::cfg::trailing_end_gap{.last_column = pa::cfg::end_gap::penalised, .last_row = pa::cfg::end_gap::free}
        )
    );

    std::string_view seq1{"ACGTGACTGACACTACGACT"};
    std::string_view seq2{"ACGTGACTGACACTACGACT"};

    auto s = unifex::just(seq1, seq2) |
             unifex::then([&] (auto const & seq1, auto const & seq2) { return aligner.compute(seq1, seq2); });

    auto alignment = unifex::sync_wait(s);

    EXPECT_EQ(alignment->score(), 80);
}
