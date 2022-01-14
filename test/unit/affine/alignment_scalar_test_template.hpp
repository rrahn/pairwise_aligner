// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <pairwise_aligner/configuration/configure_aligner.hpp>

namespace aligner = seqan::pairwise_aligner;

// ----------------------------------------------------------------------------
// Helper macro to define test values
// ----------------------------------------------------------------------------

#define DEFINE_TEST_VALUES(name, ...)  \
static auto name = [](){ return pairwise_aligner_fixture_values{__VA_ARGS__}; }();

// ----------------------------------------------------------------------------
// Values for test fixture
// ----------------------------------------------------------------------------

template <typename aligner_configurator_t, typename sequence1_t, typename sequence2_t, typename score_t>
struct pairwise_aligner_fixture_values
{
    aligner_configurator_t configurator;
    sequence1_t sequence1;
    sequence2_t sequence2;
    score_t expected_score;
};

// ----------------------------------------------------------------------------
// Fixture wrapper to access static test values
// ----------------------------------------------------------------------------

template <auto _fixture>
struct pairwise_aligner_fixture : public ::testing::Test
{
    // Method in same naming scheme as used by Google Test
    auto GetParam() -> decltype(pairwise_aligner_fixture_values{*_fixture}) const &
    {
        return *_fixture;
    }
};

// ----------------------------------------------------------------------------
// Define test fixture
// ----------------------------------------------------------------------------

template <typename fixture_t>
struct pairwise_aligner_test : public fixture_t
{};

TYPED_TEST_SUITE_P(pairwise_aligner_test);

// ----------------------------------------------------------------------------
// Define test cases
// ----------------------------------------------------------------------------

TYPED_TEST_P(pairwise_aligner_test, score)
{
    auto aligner = aligner::cfg::configure_aligner(this->GetParam().configurator);

    auto result = aligner.compute(this->GetParam().sequence1, this->GetParam().sequence2);
    EXPECT_EQ(result.score(), this->GetParam().expected_score);
}

// ----------------------------------------------------------------------------
// Register test cases
// ----------------------------------------------------------------------------

REGISTER_TYPED_TEST_SUITE_P(pairwise_aligner_test, score);
