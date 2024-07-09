// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <pairwise_aligner/configuration/gap_model_affine.hpp>
#include <pairwise_aligner/configuration/method_local.hpp>
#include <pairwise_aligner/configuration/score_model_unitary_simd_saturated.hpp>

#include "alignment_simd_test_template.hpp"

namespace local::affine::saturated_simd {

namespace aligner = seqan::pairwise_aligner;

inline constexpr size_t sequence_count = aligner::simd_score<int8_t>::size_v;

inline constexpr auto base_config =
    aligner::cfg::method_local(
        aligner::cfg::gap_model_affine(-10, -1)
    );

// ----------------------------------------------------------------------------
// Equal size
// ----------------------------------------------------------------------------

DEFINE_TEST_VALUES(equal_size_64,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_unitary_simd_saturated,
    .substitution_scores = alignment::test::simd::unitary_model<int64_t>{4, -5},
    .sequence_generation_param{sequence_count, 93, 93}
)

DEFINE_TEST_VALUES(equal_size_32,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_unitary_simd_saturated,
    .substitution_scores = alignment::test::simd::unitary_model<int32_t>{4, -5},
    .sequence_generation_param{sequence_count, 210, 210},
)

DEFINE_TEST_VALUES(equal_size_16,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_unitary_simd_saturated,
    .substitution_scores = alignment::test::simd::unitary_model<int16_t>{4, -5},
    .sequence_generation_param{sequence_count, 150, 150}
)

DEFINE_TEST_VALUES(sequence_size_1000_equal_32,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_unitary_simd_saturated,
    .substitution_scores = alignment::test::simd::unitary_model<int32_t>{4, -5},
    .sequence_generation_param{sequence_count, 1000, 1000},
)

using equal_size_types =
    ::testing::Types<
        pairwise_aligner::test::fixture<&equal_size_64>,
        pairwise_aligner::test::fixture<&equal_size_32>,
        pairwise_aligner::test::fixture<&equal_size_16>,
        pairwise_aligner::test::fixture<&sequence_size_1000_equal_32>
    >;
// ----------------------------------------------------------------------------
// Variable size
// ----------------------------------------------------------------------------

DEFINE_TEST_VALUES(variable_size_64,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_unitary_simd_saturated,
    .substitution_scores = alignment::test::simd::unitary_model<int64_t>{4, -5},
    .sequence_generation_param{sequence_count, 75, 93}
)

DEFINE_TEST_VALUES(variable_size_32,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_unitary_simd_saturated,
    .substitution_scores = alignment::test::simd::unitary_model<int32_t>{4, -5},
    .sequence_generation_param{sequence_count, 11, 200},
)

DEFINE_TEST_VALUES(variable_size_16,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_unitary_simd_saturated,
    .substitution_scores = alignment::test::simd::unitary_model<int16_t>{4, -5},
    .sequence_generation_param{sequence_count, 133, 136}
)

DEFINE_TEST_VALUES(sequence_size_1000_variable_32,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_unitary_simd_saturated,
    .substitution_scores = alignment::test::simd::unitary_model<int32_t>{4, -5},
    .sequence_generation_param{sequence_count, 900, 1100},
)

using variable_size_types =
    ::testing::Types<
        pairwise_aligner::test::fixture<&variable_size_64>,
        pairwise_aligner::test::fixture<&variable_size_32>,
        pairwise_aligner::test::fixture<&variable_size_16>,
        pairwise_aligner::test::fixture<&sequence_size_1000_variable_32>
    >;
} // global::affine::saturated_simd

INSTANTIATE_TYPED_TEST_SUITE_P(equal_size_test,
                               test_suite,
                               local::affine::saturated_simd::equal_size_types,);

INSTANTIATE_TYPED_TEST_SUITE_P(variable_size_test,
                               test_suite,
                               local::affine::saturated_simd::variable_size_types,);
