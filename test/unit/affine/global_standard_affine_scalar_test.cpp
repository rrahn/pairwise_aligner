// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string_view>

#include <pairwise_aligner/configuration/method_global.hpp>
#include <pairwise_aligner/configuration/gap_model_affine.hpp>
#include <pairwise_aligner/configuration/score_model_matrix.hpp>
#include <pairwise_aligner/configuration/score_model_unitary.hpp>

#include <pairwise_aligner/score_model/substitution_matrix.hpp>

#include "alignment_scalar_test_template.hpp"

using namespace std::literals;

namespace global::standard::affine::scalar {

namespace aligner = seqan::pairwise_aligner;

inline constexpr auto base_config =
    aligner::cfg::method_global(
        aligner::cfg::gap_model_affine(-10, -1),
        aligner::cfg::leading_end_gap{}, aligner::cfg::trailing_end_gap{}
    );

// ----------------------------------------------------------------------------
// Unitary substitution model
// ----------------------------------------------------------------------------

DEFINE_TEST_VALUES(unitary_same_sequence,
    .configurator = aligner::cfg::score_model_unitary(base_config, 4, -5),
    .sequence1 = "ACGTGACTGACACTACGACT"sv,
    .sequence2 = "ACGTGACTGACACTACGACT"sv,
    .expected_score = 80
)

DEFINE_TEST_VALUES(unitary_unequal_sequence,
    .configurator = aligner::cfg::score_model_unitary(base_config, 4, -5),
    .sequence1 = "AAAAAAAAAA"sv,
    .sequence2 = "TTTTTTTTTT"sv,
    .expected_score = -40
)

DEFINE_TEST_VALUES(unitary_related_sequence,
    .configurator = aligner::cfg::score_model_unitary(base_config, 4, -5),
    .sequence1 = "AGATCGACTAGCGAGCTACGAGCTAGC"sv,
    .sequence2 = "AGACGATCGACGAGCGACTACGTACGA"sv,
    .expected_score = 26
)

using unitary_test_types =
    ::testing::Types<
        pairwise_aligner::test::fixture<&unitary_same_sequence>,
        pairwise_aligner::test::fixture<&unitary_unequal_sequence>,
        pairwise_aligner::test::fixture<&unitary_related_sequence>
    >;

// ----------------------------------------------------------------------------
// Matrix substitution model
// ----------------------------------------------------------------------------

DEFINE_TEST_VALUES(matrix_same_sequence,
    .configurator = aligner::cfg::score_model_matrix(base_config, aligner::blosum62_standard),
    .sequence1 = "ACKLMNPQRRTVWYNMPQHIK"sv,
    .sequence2 = "ACKLMNPQRRTVWYNMPQHIK"sv,
    .expected_score = 122
)

DEFINE_TEST_VALUES(matrix_unequal_sequence,
    .configurator = aligner::cfg::score_model_matrix(base_config, aligner::blosum62_standard),
    .sequence1 = "AAAAAAAAAAAAAAAAAAAA"sv,
    .sequence2 = "CDEFGHIKLMNPQRSTVWLY"sv,
    .expected_score = -21
)

DEFINE_TEST_VALUES(matrix_sequence_different_size,
    .configurator = aligner::cfg::score_model_matrix(base_config, aligner::blosum62_standard),
    .sequence1 = "FNQSAEYPDISHCGVMQLKWRATLGT"sv,
    .sequence2 = "EIKSDVLLHRWSMKNPGNILMIDVGMQVAESYFAT"sv,
    .expected_score = -26
)

DEFINE_TEST_VALUES(matrix_short_sequence_different_size,
    .configurator = aligner::cfg::score_model_matrix(base_config, aligner::blosum62_standard),
    .sequence1 = "RKFCYMD"sv,
    .sequence2 = "GAYQW"sv,
    .expected_score = -11
)

using matrix_test_types =
    ::testing::Types<
        pairwise_aligner::test::fixture<&matrix_same_sequence>,
        pairwise_aligner::test::fixture<&matrix_unequal_sequence>,
        pairwise_aligner::test::fixture<&matrix_sequence_different_size>,
        pairwise_aligner::test::fixture<&matrix_short_sequence_different_size>
    >;

} // namespace global::standard::affine::scalar

INSTANTIATE_TYPED_TEST_SUITE_P(unitary_test,
                               test_suite,
                               global::standard::affine::scalar::unitary_test_types,);

INSTANTIATE_TYPED_TEST_SUITE_P(matrix_test,
                               test_suite,
                               global::standard::affine::scalar::matrix_test_types,);
