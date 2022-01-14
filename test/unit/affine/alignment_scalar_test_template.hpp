// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <pairwise_aligner/configuration/configure_aligner.hpp>

#include "fixture_base.hpp"

// ----------------------------------------------------------------------------
// Helper macro to define test values
// ----------------------------------------------------------------------------

#define DEFINE_TEST_VALUES(name, ...)  \
static auto name = [](){ return alignment::test::scalar::values{__VA_ARGS__}; }();

namespace alignment::test::scalar {

// ----------------------------------------------------------------------------
// Values for test
// ----------------------------------------------------------------------------

template <typename aligner_configurator_t, typename sequence1_t, typename sequence2_t, typename score_t>
struct values
{
    aligner_configurator_t configurator;
    sequence1_t sequence1;
    sequence2_t sequence2;
    score_t expected_score;
};

// ----------------------------------------------------------------------------
// Define test
// ----------------------------------------------------------------------------

template <typename fixture_t>
struct test : public fixture_t
{};

} // namespec alignment::test::scalar

template <typename fixture_t>
using test_suite = alignment::test::scalar::test<fixture_t>;

TYPED_TEST_SUITE_P(test_suite);

// ----------------------------------------------------------------------------
// Define test cases
// ----------------------------------------------------------------------------

TYPED_TEST_P(test_suite, score)
{
    auto aligner = seqan::pairwise_aligner::cfg::configure_aligner(this->GetParam().configurator);

    auto result = aligner.compute(this->GetParam().sequence1, this->GetParam().sequence2);
    EXPECT_EQ(result.score(), this->GetParam().expected_score);
}

// ----------------------------------------------------------------------------
// Register test cases
// ----------------------------------------------------------------------------

REGISTER_TYPED_TEST_SUITE_P(test_suite, score);
