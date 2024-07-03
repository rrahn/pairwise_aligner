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

#include <pairwise_aligner/score_model/substitution_matrix.hpp>

#include "alignment_scalar_test_template.hpp"

using namespace std::literals;

namespace global::standard::affine::scalar::matrix {

namespace aligner = seqan::pairwise_aligner;

inline constexpr auto base_config =
    aligner::cfg::method_global(
        aligner::cfg::gap_model_affine(-10, -1),
        aligner::cfg::leading_end_gap{}, aligner::cfg::trailing_end_gap{}
    );

DEFINE_TEST_VALUES(same_sequence,
    .configurator = aligner::cfg::score_model_matrix(base_config, aligner::blosum62_standard<>),
    .sequence1 = "ACKLMNPQRRTVWYNMPQHIK"sv,
    .sequence2 = "ACKLMNPQRRTVWYNMPQHIK"sv,
    .expected_score = 122
)

DEFINE_TEST_VALUES(unequal_sequence,
    .configurator = aligner::cfg::score_model_matrix(base_config, aligner::blosum62_standard<>),
    .sequence1 = "AAAAAAAAAAAAAAAAAAAA"sv,
    .sequence2 = "CDEFGHIKLMNPQRSTVWLY"sv,
    .expected_score = -21
)

DEFINE_TEST_VALUES(sequence_different_size,
    .configurator = aligner::cfg::score_model_matrix(base_config, aligner::blosum62_standard<>),
    .sequence1 = "FNQSAEYPDISHCGVMQLKWRATLGT"sv,
    .sequence2 = "EIKSDVLLHRWSMKNPGNILMIDVGMQVAESYFAT"sv,
    .expected_score = -26
)

DEFINE_TEST_VALUES(short_sequence_different_size,
    .configurator = aligner::cfg::score_model_matrix(base_config, aligner::blosum62_standard<>),
    .sequence1 = "RKFCYMD"sv,
    .sequence2 = "GAYQW"sv,
    .expected_score = -11
)

using test_types =
    ::testing::Types<
        pairwise_aligner::test::fixture<&same_sequence>,
        pairwise_aligner::test::fixture<&unequal_sequence>,
        pairwise_aligner::test::fixture<&sequence_different_size>,
        pairwise_aligner::test::fixture<&short_sequence_different_size>
    >;

} // namespace global::standard::affine::scalar::matrix

INSTANTIATE_TYPED_TEST_SUITE_P(test,
                               test_suite,
                               global::standard::affine::scalar::matrix::test_types,);
