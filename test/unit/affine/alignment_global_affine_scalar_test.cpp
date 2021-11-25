// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string_view>

#include <pairwise_aligner/affine/affine_dp_algorithm.hpp>

TEST(affine_test, all_match)
{
    seqan::pairwise_aligner::pairwise_aligner_affine<> aligner{};

    std::string_view seq1{"ACGTGACTGACACTACGACT"};
    std::string_view seq2{"ACGTGACTGACACTACGACT"};

    EXPECT_EQ((aligner.compute(seq1, seq2)), 80);
}

TEST(affine_test, all_mismatch)
{
    seqan::pairwise_aligner::pairwise_aligner_affine<> aligner{};

    std::string_view seq1{"AAAAAAAAAA"};
    std::string_view seq2{"TTTTTTTTTT"};

    EXPECT_EQ((aligner.compute(seq1, seq2)), -40);
}
