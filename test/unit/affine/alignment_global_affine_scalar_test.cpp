// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string_view>

#include <pairwise_aligner/affine/affine_dp_algorithm.hpp>
#include <pairwise_aligner/affine/gap_model_affine.hpp>
#include <pairwise_aligner/affine/initialisation_strategy_affine.hpp>

TEST(affine_test, all_match)
{
    namespace pa = seqan::pairwise_aligner;
    pa::gap_model_affine<int32_t> gap_model{-10, -1};
    pa::initialisation_strategy_affine init{gap_model,
                                            pa::dp_initialisation_rule::regular,
                                            pa::dp_initialisation_rule::regular};
    pa::pairwise_aligner_affine aligner{gap_model, init};

    std::string_view seq1{"ACGTGACTGACACTACGACT"};
    std::string_view seq2{"ACGTGACTGACACTACGACT"};

    EXPECT_EQ((aligner.compute(seq1, seq2)), 80);
}

TEST(affine_test, all_mismatch)
{
    namespace pa = seqan::pairwise_aligner;

    pa::gap_model_affine<int32_t> gap_model{-10, -1};
    pa::initialisation_strategy_affine init{gap_model,
                                            pa::dp_initialisation_rule::regular,
                                            pa::dp_initialisation_rule::regular};
    pa::pairwise_aligner_affine aligner{gap_model, init};

    std::string_view seq1{"AAAAAAAAAA"};
    std::string_view seq2{"TTTTTTTTTT"};

    EXPECT_EQ((aligner.compute(seq1, seq2)), -40);
}
