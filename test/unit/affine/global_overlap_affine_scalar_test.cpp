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

namespace global::overlap::affine::scalar {

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

DEFINE_TEST_VALUES(all_match,
    .configurator = aligner::cfg::score_model_unitary(base_config, 4, -5),
    .sequence1 = "ACGTGACTGACACTACGACT"sv,
    .sequence2 = "ACGTGACTGACACTACGACT"sv,
    .expected_score = 80
)

DEFINE_TEST_VALUES(overlap_left_match,
    .configurator = aligner::cfg::score_model_unitary(base_config, 4, -5),
    .sequence1 = "TTTTTCTGACACACT"sv,
    .sequence2 = "CTGCACACGACGGGGGG"sv,
    .expected_score = 16
)

DEFINE_TEST_VALUES(overlap_right_match,
    .configurator = aligner::cfg::score_model_unitary(base_config, 4, -5),
    .sequence1 = "CTGCACACGACGGGGGG"sv,
    .sequence2 = "TTTTTCTGACACACT"sv,
    .expected_score = 16
)

using test_types =
    ::testing::Types<
        alignment::test::fixture<&all_match>,
        alignment::test::fixture<&overlap_left_match>,
        alignment::test::fixture<&overlap_right_match>
    >;

} // namespace global::overlap::affine::scalar

INSTANTIATE_TYPED_TEST_SUITE_P(test,
                               test_suite,
                               global::overlap::affine::scalar::test_types,);
