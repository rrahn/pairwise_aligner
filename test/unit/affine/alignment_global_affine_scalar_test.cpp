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
#include <pairwise_aligner/interface/interface_one_to_one_single.hpp>
#include <pairwise_aligner/score_model/score_model_unitary.hpp>
#include <pairwise_aligner/pairwise_aligner.hpp>

TEST(affine_test, all_match)
{
    namespace pa = seqan::pairwise_aligner;

    pa::score_model_unitary<int32_t> score_model{4, -5};
    using score_t = typename decltype(score_model)::score_type;
    pa::gap_model_affine<score_t> gap_model{-10, -1};
    pa::initialisation_strategy_affine init{gap_model,
                                            pa::dp_initialisation_rule::regular,
                                            pa::dp_initialisation_rule::regular};

    using dp_vector_column_t = pa::intermediate_dp_vector<pa::affine_cell<score_t, pa::dp_vector_order::column>>;
    using dp_vector_row_t = pa::intermediate_dp_vector<pa::affine_cell<score_t, pa::dp_vector_order::row>>;
    using dp_algorithm_t = decltype(pa::pairwise_aligner_affine{score_model, gap_model, init});
    using aligner_t = pa::interface_one_to_one_single<dp_algorithm_t, dp_vector_column_t, dp_vector_row_t>;

    aligner_t aligner{score_model, gap_model, init};

    std::string_view seq1{"ACGTGACTGACACTACGACT"};
    std::string_view seq2{"ACGTGACTGACACTACGACT"};

    EXPECT_EQ((aligner.compute(seq1, seq2)), 80);
}

TEST(affine_test, all_mismatch)
{
    namespace pa = seqan::pairwise_aligner;

    pa::score_model_unitary<int32_t> score_model{4, -5};
    using score_t = typename decltype(score_model)::score_type;
    pa::gap_model_affine<score_t> gap_model{-10, -1};
    pa::initialisation_strategy_affine init{gap_model,
                                            pa::dp_initialisation_rule::regular,
                                            pa::dp_initialisation_rule::regular};

    using score_t = int32_t;
    using dp_vector_column_t = pa::intermediate_dp_vector<pa::affine_cell<score_t, pa::dp_vector_order::column>>;
    using dp_vector_row_t = pa::intermediate_dp_vector<pa::affine_cell<score_t, pa::dp_vector_order::row>>;
    using dp_algorithm_t = decltype(pa::pairwise_aligner_affine{score_model, gap_model, init});
    using aligner_t = pa::interface_one_to_one_single<dp_algorithm_t, dp_vector_column_t, dp_vector_row_t>;

    aligner_t aligner{score_model, gap_model, init};

    std::string_view seq1{"AAAAAAAAAA"};
    std::string_view seq2{"TTTTTTTTTT"};

    EXPECT_EQ((aligner.compute(seq1, seq2)), -40);
}
