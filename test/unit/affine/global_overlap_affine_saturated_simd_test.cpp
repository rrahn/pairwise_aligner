// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <pairwise_aligner/configuration/gap_model_affine.hpp>
#include <pairwise_aligner/configuration/method_global.hpp>
#include <pairwise_aligner/configuration/score_model_unitary_simd_saturated.hpp>

#include "alignment_simd_test_template.hpp"

namespace global::overlap::affine::saturated_simd {

namespace aligner = seqan::pairwise_aligner;

inline constexpr size_t sequence_count = aligner::simd_score<int8_t>::size;

inline constexpr auto base_config =
    aligner::cfg::method_global(
        aligner::cfg::gap_model_affine(-10, -1),
        aligner::cfg::leading_end_gap{aligner::cfg::end_gap::free, aligner::cfg::end_gap::free},
        aligner::cfg::trailing_end_gap{aligner::cfg::end_gap::free, aligner::cfg::end_gap::free}
    );

// ----------------------------------------------------------------------------
// Equal size
// ----------------------------------------------------------------------------

DEFINE_TEST_VALUES(equal_size_64, int64_t,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_unitary_simd_saturated,
    .substitution_scores{4, -5},
    .sequence_generation_param{sequence_count, 93, 93}
)

DEFINE_TEST_VALUES(equal_size_32, int32_t,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_unitary_simd_saturated,
    .substitution_scores{4, -5},
    .sequence_generation_param{sequence_count, 210, 210},
)

DEFINE_TEST_VALUES(equal_size_16, int16_t,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_unitary_simd_saturated,
    .substitution_scores{4, -5},
    .sequence_generation_param{sequence_count, 150, 150}
)

DEFINE_TEST_VALUES(equal_size_8, int8_t,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_unitary_simd_saturated,
    .substitution_scores{4, -5},
    .sequence_generation_param{sequence_count, 25, 25},
)

DEFINE_TEST_VALUES(sequence_size_1000_equal_32, int32_t,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_unitary_simd_saturated,
    .substitution_scores{4, -5},
    .sequence_generation_param{sequence_count, 1000, 1000},
)

using equal_size_types =
    ::testing::Types<
        pairwise_aligner::test::fixture<&equal_size_64>,
        pairwise_aligner::test::fixture<&equal_size_32>,
        pairwise_aligner::test::fixture<&equal_size_16>,
        pairwise_aligner::test::fixture<&equal_size_8>,
        pairwise_aligner::test::fixture<&sequence_size_1000_equal_32>
    >;
// ----------------------------------------------------------------------------
// Variable size
// ----------------------------------------------------------------------------

DEFINE_TEST_VALUES(variable_size_64, int64_t,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_unitary_simd_saturated,
    .substitution_scores{4, -5},
    .sequence_generation_param{sequence_count, 75, 93}
)

DEFINE_TEST_VALUES(variable_size_32, int32_t,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_unitary_simd_saturated,
    .substitution_scores{4, -5},
    .sequence_generation_param{sequence_count, 11, 200},
)

DEFINE_TEST_VALUES(variable_size_16, int16_t,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_unitary_simd_saturated,
    .substitution_scores{4, -5},
    .sequence_generation_param{sequence_count, 133, 136}
)

DEFINE_TEST_VALUES(variable_size_8, int8_t,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_unitary_simd_saturated,
    .substitution_scores{4, -5},
    .sequence_generation_param{sequence_count, 10, 15},
)

DEFINE_TEST_VALUES(sequence_size_1000_variable_32, int32_t,
    .base_configurator = base_config,
    .score_configurator = aligner::cfg::score_model_unitary_simd_saturated,
    .substitution_scores{4, -5},
    .sequence_generation_param{sequence_count, 900, 1100},
)

using variable_size_types =
    ::testing::Types<
        pairwise_aligner::test::fixture<&variable_size_64>,
        pairwise_aligner::test::fixture<&variable_size_32>,
        pairwise_aligner::test::fixture<&variable_size_16>,
        pairwise_aligner::test::fixture<&variable_size_8>,
        pairwise_aligner::test::fixture<&sequence_size_1000_variable_32>
    >;
} // global::affine::saturated_simd

INSTANTIATE_TYPED_TEST_SUITE_P(equal_size_test,
                               test_suite,
                               global::overlap::affine::saturated_simd::equal_size_types,);

INSTANTIATE_TYPED_TEST_SUITE_P(variable_size_test,
                               test_suite,
                               global::overlap::affine::saturated_simd::variable_size_types,);
