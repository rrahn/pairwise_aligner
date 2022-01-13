// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/std/algorithm>
#include <seqan3/std/concepts>
#include <random>
#include <string>
#include <vector>

#include <pairwise_aligner/configuration/configure_aligner.hpp>
#include <pairwise_aligner/configuration/score_model_unitary.hpp>

namespace aligner = seqan::pairwise_aligner;

// ----------------------------------------------------------------------------
// Helper macro to define test values
// ----------------------------------------------------------------------------

#define DEFINE_TEST_VALUES(name, score_t, ...)  \
static auto name = [](){ return pairwise_aligner_fixture_values{.score_type = score_t{}__VA_OPT__(,) __VA_ARGS__}; }();

// ----------------------------------------------------------------------------
// Values for test fixture
// ----------------------------------------------------------------------------

template <typename score_t, typename base_configurator_t, typename score_configurator_t>
struct pairwise_aligner_fixture_values
{
    score_t score_type;
    base_configurator_t base_configurator;
    score_configurator_t score_configurator;
    std::pair<int32_t, int32_t> substitution_scores;
    std::tuple<std::size_t, std::size_t, std::size_t> sequence_generation_param;
    unsigned seed{42};
};

// ----------------------------------------------------------------------------
// Fixture wrapper to access static test values
// ----------------------------------------------------------------------------

template <auto _fixture>
struct pairwise_aligner_fixture : public ::testing::Test
{
    using scalar_score_type = decltype(_fixture->score_type);

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
{
    std::vector<std::string> sequence_collection1{};
    std::vector<std::string> sequence_collection2{};

    void SetUp() override
    {
        std::mt19937 random_engine{this->GetParam().seed};
        generate_sequences(random_engine);
    }

private:

    template <typename random_engine_t>
    void generate_sequences(random_engine_t & random_engine)
    {
        auto [count, min_size, max_size] = this->GetParam().sequence_generation_param;
        sequence_collection1.resize(count);
        sequence_collection2.resize(count);

        std::uniform_int_distribution<size_t> sequence_size_distribution{min_size, max_size};
        constexpr std::array dna_symbol_table{'A', 'C', 'G', 'T'};

        auto generate_sequence = [&] (auto & sequence) {
            sequence.resize(sequence_size_distribution(random_engine));
            std::ranges::generate(sequence, [&] () { return dna_symbol_table[std::rand() % dna_symbol_table.size()]; });
        };

        for (size_t i = 0; i < sequence_collection1.size(); ++i) {
            generate_sequence(sequence_collection1[i]);
            generate_sequence(sequence_collection2[i]);
        }
    }
};

TYPED_TEST_SUITE_P(pairwise_aligner_test);

// ----------------------------------------------------------------------------
// Define test cases
// ----------------------------------------------------------------------------

TYPED_TEST_P(pairwise_aligner_test, score)
{
    auto [match_score, mismatch_score] = this->GetParam().substitution_scores;

    auto scalar_aligner = aligner::cfg::configure_aligner(
        aligner::cfg::score_model_unitary(this->GetParam().base_configurator, match_score, mismatch_score)
    );

    auto simd_aligner = aligner::cfg::configure_aligner(
        this->GetParam().score_configurator(this->GetParam().base_configurator,
                                            static_cast<TestFixture::scalar_score_type>(match_score),
                                            static_cast<TestFixture::scalar_score_type>(mismatch_score))
    );

    auto simd_results = simd_aligner.compute(this->sequence_collection1, this->sequence_collection2);
    for (auto result : simd_results) {
        EXPECT_EQ(result.score(), (scalar_aligner.compute(result.sequence1(), result.sequence2()).score()));
    }
}

// ----------------------------------------------------------------------------
// Register test cases
// ----------------------------------------------------------------------------

REGISTER_TYPED_TEST_SUITE_P(pairwise_aligner_test, score);