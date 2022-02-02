// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string_view>

#include <pairwise_aligner/configuration/method_local.hpp>
#include <pairwise_aligner/configuration/gap_model_affine.hpp>
#include <pairwise_aligner/configuration/score_model_matrix.hpp>
#include <pairwise_aligner/configuration/score_model_unitary.hpp>
#include <pairwise_aligner/score_model/substitution_matrix.hpp>

#include "alignment_scalar_test_template.hpp"

using namespace std::literals;

namespace local::affine::scalar {

// ----------------------------------------------------------------------------
// Base aligner configuration
// ----------------------------------------------------------------------------

namespace aligner = seqan::pairwise_aligner;

inline constexpr auto base_config =
    aligner::cfg::method_local(
        aligner::cfg::gap_model_affine(-10, -1)
    );

// ----------------------------------------------------------------------------
// Test definitions
// ----------------------------------------------------------------------------

DEFINE_TEST_VALUES(same_sequence,
    .configurator = aligner::cfg::score_model_unitary(base_config, 4, -5),
    .sequence1 = "TTTTTTTTTTTTTTT"sv,
    .sequence2 = "TTTTTTTTTTTTTTT"sv,
    .expected_score = 60
)

DEFINE_TEST_VALUES(local_infix_sequence,
    .configurator = aligner::cfg::score_model_unitary(base_config, 4, -5),
    .sequence1 = "CTGCACACGACGGGGGG"sv,
    .sequence2 = "TTTTTCTGACACACT"sv,
    .expected_score = 21
)

DEFINE_TEST_VALUES(unequal_sequence,
    .configurator = aligner::cfg::score_model_unitary(base_config, 4, -5),
    .sequence1 = "AAAAAAAAAA"sv,
    .sequence2 = "TTTTTTTTTT"sv,
    .expected_score = 0
)

DEFINE_TEST_VALUES(single_match,
    .configurator = aligner::cfg::score_model_unitary(base_config, 4, -5),
    .sequence1 = "AAAAATAAAA"sv,
    .sequence2 = "CCCCCTCCCC"sv,
    .expected_score = 4
)

DEFINE_TEST_VALUES(infix_with_one_mismatch,
    .configurator = aligner::cfg::score_model_unitary(base_config, 4, -5),
    .sequence1 = "AACCGGTTTAACCGGTT"sv,
    .sequence2 = "ACGTCTACGTA"sv,
    .expected_score = 11
)

DEFINE_TEST_VALUES(prefix_with_suffix,
    .configurator = aligner::cfg::score_model_unitary(base_config, 4, -5),
    .sequence1 = "AAAAAATCCCCCC"sv,
    .sequence2 = "CCCCCCTAAAAAA"sv,
    .expected_score = 24
)

DEFINE_TEST_VALUES(infix_starting_in_first_row,
    .configurator = aligner::cfg::score_model_unitary(
                        aligner::cfg::method_local(
                            aligner::cfg::gap_model_affine(-1, -1)
                        ),
                        2, -1
                    ),
    .sequence1 = "ATAAGCGTCTCG"sv,
    .sequence2 = "TCATAGAGTTGC"sv,
    .expected_score = 9
)

DEFINE_TEST_VALUES(matrix_infix,
    .configurator = aligner::cfg::score_model_matrix(base_config, aligner::blosum62_extended),
    .sequence1 = "ALIGATOR"sv,
    .sequence2 = "GALORA"sv,
    .expected_score = 13
)

DEFINE_TEST_VALUES(matrix_empty_sequence,
    .configurator = aligner::cfg::score_model_matrix(base_config, aligner::blosum62_extended),
    .sequence1 = "ALIGATOR"sv,
    .sequence2 = ""sv,
    .expected_score = std::numeric_limits<int32_t>::lowest()
)

using test_types =
    ::testing::Types<
        alignment::test::fixture<&same_sequence>,
        alignment::test::fixture<&local_infix_sequence>,
        alignment::test::fixture<&unequal_sequence>,
        alignment::test::fixture<&single_match>,
        alignment::test::fixture<&infix_with_one_mismatch>,
        alignment::test::fixture<&prefix_with_suffix>,
        alignment::test::fixture<&infix_starting_in_first_row>,
        alignment::test::fixture<&matrix_infix>,
        alignment::test::fixture<&matrix_empty_sequence>
    >;

} // namespace local::affine::scalar

INSTANTIATE_TYPED_TEST_SUITE_P(test,
                               test_suite,
                               local::affine::scalar::test_types,);
