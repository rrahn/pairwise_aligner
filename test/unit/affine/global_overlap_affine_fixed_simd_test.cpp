// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <pairwise_aligner/configuration/gap_model_affine.hpp>
#include <pairwise_aligner/configuration/method_global.hpp>
#include <pairwise_aligner/configuration/score_model_unitary_simd.hpp>

#include "alignment_simd_test_template.hpp"

namespace global::overlap::affine::fixed_simd {

namespace aligner = seqan::pairwise_aligner;

inline constexpr auto base_config =
    aligner::cfg::method_global(
        aligner::cfg::gap_model_affine(-10, -1),
        aligner::cfg::leading_end_gap{aligner::cfg::end_gap::free, aligner::cfg::end_gap::free},
        aligner::cfg::trailing_end_gap{aligner::cfg::end_gap::free, aligner::cfg::end_gap::free}
    );

// ----------------------------------------------------------------------------
// Equal size
// ----------------------------------------------------------------------------

DEFINE_TEST_VALUES(equal_size_64,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_unitary_simd,
    .substitution_scores = alignment::test::simd::unitary_model<int64_t>{4, -5},
    .sequence_generation_param{aligner::simd_score<int64_t>::size_v, 93, 93}
)

DEFINE_TEST_VALUES(equal_size_32,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_unitary_simd,
    .substitution_scores = alignment::test::simd::unitary_model<int32_t>{4, -5},
    .sequence_generation_param{aligner::simd_score<int32_t>::size_v, 210, 210}
)

DEFINE_TEST_VALUES(equal_size_16,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_unitary_simd,
    .substitution_scores = alignment::test::simd::unitary_model<int16_t>{4, -5},
    .sequence_generation_param{aligner::simd_score<int16_t>::size_v, 150, 150}
)

DEFINE_TEST_VALUES(equal_size_8,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_unitary_simd,
    .substitution_scores = alignment::test::simd::unitary_model<int8_t>{4, -5},
    .sequence_generation_param{aligner::simd_score<int8_t>::size_v, 25, 25}
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
    .score_configurator = aligner::cfg::score_model_unitary_simd,
    .substitution_scores = alignment::test::simd::unitary_model<int64_t>{4, -5},
    .sequence_generation_param{aligner::simd_score<int64_t>::size_v, 75, 93}
)

DEFINE_TEST_VALUES(variable_size_32,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_unitary_simd,
    .substitution_scores = alignment::test::simd::unitary_model<int32_t>{4, -5},
    .sequence_generation_param{aligner::simd_score<int32_t>::size_v, 11, 200}
)

DEFINE_TEST_VALUES(variable_size_16,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_unitary_simd,
    .substitution_scores = alignment::test::simd::unitary_model<int16_t>{4, -5},
    .sequence_generation_param{aligner::simd_score<int16_t>::size_v, 133, 136}
)

DEFINE_TEST_VALUES(variable_size_8,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_unitary_simd,
    .substitution_scores = alignment::test::simd::unitary_model<int8_t>{4, -5},
    .sequence_generation_param{aligner::simd_score<int8_t>::size_v, 10, 15}
)

using variable_size_types =
    ::testing::Types<
        pairwise_aligner::test::fixture<&variable_size_64>,
        pairwise_aligner::test::fixture<&variable_size_32>,
        pairwise_aligner::test::fixture<&variable_size_16>,
        pairwise_aligner::test::fixture<&variable_size_8>
    >;
} // global::affine::fixed_simd

INSTANTIATE_TYPED_TEST_SUITE_P(equal_size_test,
                               test_suite,
                               global::overlap::affine::fixed_simd::equal_size_types,);

INSTANTIATE_TYPED_TEST_SUITE_P(variable_size_test,
                               test_suite,
                               global::overlap::affine::fixed_simd::variable_size_types,);
