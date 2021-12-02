// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string_view>

#include <pairwise_aligner/configuration/configure_aligner.hpp>
#include <pairwise_aligner/configuration/gap_model_affine.hpp>
#include <pairwise_aligner/configuration/score_model_unitary_simd.hpp>

TEST(affine_test, all_match)
{
    namespace pa = seqan::pairwise_aligner;

    auto aligner = pa::cfg::configure_aligner(
        pa::cfg::gap_model_affine(
            pa::cfg::score_model_unitary_simd(static_cast<int16_t>(4), static_cast<int16_t>(-5)),
            -10, -1
        )
    );

    std::string_view seq1{"ACGTGACTGACACTACGACT"};
    std::string_view seq2{"ACGTGACTGACACTACGACT"};

    std::vector<std::string_view> seq1_collection{};
    seq1_collection.resize(16, seq1);

    std::vector<std::string_view> seq2_collection{};
    seq2_collection.resize(16, seq2);

    auto res = aligner.compute(seq1_collection, seq2_collection);
    for (auto i = 0; i < 16; ++i)
        EXPECT_EQ(res[i], 80);
}

TEST(affine_test, all_mismatch)
{
    namespace pa = seqan::pairwise_aligner;

    auto aligner = pa::cfg::configure_aligner(
        pa::cfg::gap_model_affine(
            pa::cfg::score_model_unitary_simd(static_cast<int16_t>(4), static_cast<int16_t>(-5)),
            -10, -1
        )
    );

    std::string_view seq1{"AAAAAAAAAA"};
    std::string_view seq2{"TTTTTTTTTT"};

    std::vector<std::string_view> seq1_collection{};
    seq1_collection.resize(16, seq1);

    std::vector<std::string_view> seq2_collection{};
    seq2_collection.resize(16, seq2);

    auto res = aligner.compute(seq1_collection, seq2_collection);
    for (auto i = 0; i < 16; ++i)
        EXPECT_EQ(res[i], -40);
}
