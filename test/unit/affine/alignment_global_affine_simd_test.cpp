// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string_view>

#include <pairwise_aligner/affine/affine_dp_algorithm.hpp>
#include <pairwise_aligner/affine/affine_gap_model.hpp>
#include <pairwise_aligner/affine/affine_initialisation_strategy.hpp>
#include <pairwise_aligner/simd_score_type.hpp>
#include <pairwise_aligner/interface/interface_one_to_one_bulk.hpp>
#include <pairwise_aligner/pairwise_aligner.hpp>
#include <pairwise_aligner/score_model/score_model_unitary.hpp>

TEST(affine_test, all_match)
{
    namespace pa = seqan::pairwise_aligner;

    using score_t = int16_t;
    using simd_score_t = pa::simd_score<score_t>;
    pa::score_model_unitary<simd_score_t> score_model{simd_score_t{4}, simd_score_t{-5}};

    pa::affine_gap_model<score_t> gap_model{score_t{-10}, score_t{-1}};
    pa::affine_initialisation_strategy init{gap_model,
                                            pa::dp_initialisation_rule::regular,
                                            pa::dp_initialisation_rule::regular};

    using dp_vector_column_t = pa::simd_intermediate_dp_vector<pa::affine_cell<simd_score_t, pa::dp_vector_order::column>>;
    using dp_vector_row_t = pa::simd_intermediate_dp_vector<pa::affine_cell<simd_score_t, pa::dp_vector_order::row>>;
    using dp_algorithm_t = decltype(pa::pairwise_aligner_affine{score_model, gap_model, init});
    using aligner_t = pa::interface_one_to_one_bulk<dp_algorithm_t, dp_vector_column_t, dp_vector_row_t>;

    aligner_t aligner{score_model, gap_model, init};

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

    using score_t = int16_t;
    using simd_score_t = pa::simd_score<score_t>;
    pa::score_model_unitary<simd_score_t> score_model{simd_score_t{4}, simd_score_t{-5}};

    pa::affine_gap_model<score_t> gap_model{score_t{-10}, score_t{-1}};
    pa::affine_initialisation_strategy init{gap_model,
                                            pa::dp_initialisation_rule::regular,
                                            pa::dp_initialisation_rule::regular};

    using dp_vector_column_t = pa::simd_intermediate_dp_vector<pa::affine_cell<simd_score_t, pa::dp_vector_order::column>>;
    using dp_vector_row_t = pa::simd_intermediate_dp_vector<pa::affine_cell<simd_score_t, pa::dp_vector_order::row>>;
    using dp_algorithm_t = decltype(pa::pairwise_aligner_affine{score_model, gap_model, init});
    using aligner_t = pa::interface_one_to_one_bulk<dp_algorithm_t, dp_vector_column_t, dp_vector_row_t>;

    aligner_t aligner{score_model, gap_model, init};

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
