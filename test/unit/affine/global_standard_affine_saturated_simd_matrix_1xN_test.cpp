// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <pairwise_aligner/configuration/gap_model_affine.hpp>
#include <pairwise_aligner/configuration/method_global.hpp>
#include <pairwise_aligner/configuration/score_model_matrix_simd_saturated_1xN.hpp>

#include <pairwise_aligner/score_model/substitution_matrix.hpp>

#include "alignment_simd_test_template.hpp"

namespace global::standard::affine::saturated_simd::matrix {

namespace aligner = seqan::pairwise_aligner;

inline constexpr size_t sequence_count = aligner::simd_score<int8_t>::size_v;

inline constexpr auto base_config =
    aligner::cfg::method_global(
        aligner::cfg::gap_model_affine(-10, -1),
        aligner::cfg::leading_end_gap{}, aligner::cfg::trailing_end_gap{}
    );

// ----------------------------------------------------------------------------
// Equal size
// ----------------------------------------------------------------------------

DEFINE_TEST_VALUES(equal_size_64,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_matrix_simd_saturated_1xN,
    .substitution_scores = alignment::test::simd::matrix_model{aligner::blosum62_standard<int64_t>},
    .sequence_generation_param{sequence_count, 93, 93},
    .one_vs_many = std::true_type{},
)

DEFINE_TEST_VALUES(equal_size_32,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_matrix_simd_saturated_1xN,
    .substitution_scores = alignment::test::simd::matrix_model{aligner::blosum62_standard<int32_t>},
    .sequence_generation_param{sequence_count, 210, 210},
    .one_vs_many = std::true_type{},
)

DEFINE_TEST_VALUES(equal_size_16,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_matrix_simd_saturated_1xN,
    .substitution_scores = alignment::test::simd::matrix_model{aligner::blosum62_standard<int16_t>},
    .sequence_generation_param{sequence_count, 150, 150},
    .one_vs_many = std::true_type{},
)

DEFINE_TEST_VALUES(equal_size_8,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_matrix_simd_saturated_1xN,
    .substitution_scores = alignment::test::simd::matrix_model{aligner::blosum62_standard<int8_t>},
    .sequence_generation_param{sequence_count, 25, 25},
    .one_vs_many = std::true_type{},
)

using equal_size_types =
    ::testing::Types<
        pairwise_aligner::test::fixture<&equal_size_64>,
        pairwise_aligner::test::fixture<&equal_size_32>,
        pairwise_aligner::test::fixture<&equal_size_16>,
        pairwise_aligner::test::fixture<&equal_size_8>
    >;
// ----------------------------------------------------------------------------
// Variable size
// ----------------------------------------------------------------------------

DEFINE_TEST_VALUES(variable_size_64,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_matrix_simd_saturated_1xN,
    .substitution_scores = alignment::test::simd::matrix_model{aligner::blosum62_standard<int64_t>},
    .sequence_generation_param{sequence_count, 75, 93},
    .one_vs_many = std::true_type{},
)

DEFINE_TEST_VALUES(variable_size_32,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_matrix_simd_saturated_1xN,
    .substitution_scores = alignment::test::simd::matrix_model{aligner::blosum62_standard<int32_t>},
    .sequence_generation_param{sequence_count, 11, 200},
    .one_vs_many = std::true_type{},
)

DEFINE_TEST_VALUES(variable_size_16,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_matrix_simd_saturated_1xN,
    .substitution_scores = alignment::test::simd::matrix_model{aligner::blosum62_standard<int16_t>},
    .sequence_generation_param{sequence_count, 133, 136},
    .one_vs_many = std::true_type{},
)

DEFINE_TEST_VALUES(variable_size_8,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_matrix_simd_saturated_1xN,
    .substitution_scores = alignment::test::simd::matrix_model{aligner::blosum62_standard<int8_t>},
    .sequence_generation_param{sequence_count, 10, 15},
    .one_vs_many = std::true_type{},
)

using variable_size_types =
    ::testing::Types<
        pairwise_aligner::test::fixture<&variable_size_64>,
        pairwise_aligner::test::fixture<&variable_size_32>,
        pairwise_aligner::test::fixture<&variable_size_16>,
        pairwise_aligner::test::fixture<&variable_size_8>
    >;
} // global::affine::saturated_simd

INSTANTIATE_TYPED_TEST_SUITE_P(equal_size_test,
                               test_suite,
                               global::standard::affine::saturated_simd::matrix::equal_size_types,);

INSTANTIATE_TYPED_TEST_SUITE_P(variable_size_test,
                               test_suite,
                               global::standard::affine::saturated_simd::matrix::variable_size_types,);
