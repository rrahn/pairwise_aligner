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
#include <pairwise_aligner/configuration/score_model_unitary.hpp>

#include "alignment_scalar_test_template.hpp"

using namespace std::literals;

namespace global::semi_first::affine::scalar {

namespace aligner = seqan::pairwise_aligner;

inline constexpr auto base_config =
    aligner::cfg::method_global(
        aligner::cfg::gap_model_affine(-10, -1),
        aligner::cfg::leading_end_gap{.first_column = aligner::cfg::end_gap::free },
        aligner::cfg::trailing_end_gap{.last_column = aligner::cfg::end_gap::free }
    );

// ----------------------------------------------------------------------------
// Equal size
// ----------------------------------------------------------------------------

DEFINE_TEST_VALUES(infix_match,
    .configurator = aligner::cfg::score_model_unitary(base_config, 4, -5),
    .sequence1 = "TTTTTACGTATGTCCCCC"sv,
    .sequence2 = "ACGTAAAACGT"sv,
    .expected_score = 10
)

DEFINE_TEST_VALUES(prefix_match,
    .configurator = aligner::cfg::score_model_unitary(base_config, 4, -5),
    .sequence1 = "AAGACTACGGGGGG"sv,
    .sequence2 = "ACGAGGCTAC"sv,
    .expected_score = 11
)

DEFINE_TEST_VALUES(suffix_match,
    .configurator = aligner::cfg::score_model_unitary(base_config, 4, -5),
    .sequence1 = "GGGGGGAAGACTAC"sv,
    .sequence2 = "ACGAGGCTAC"sv,
    .expected_score = 11
)

using test_types =
    ::testing::Types<
        alignment::test::fixture<&infix_match>,
        alignment::test::fixture<&prefix_match>,
        alignment::test::fixture<&suffix_match>
    >;

} // namespace global::semi_first::affine::scalar

INSTANTIATE_TYPED_TEST_SUITE_P(test,
                               test_suite,
                               global::semi_first::affine::scalar::test_types,);
