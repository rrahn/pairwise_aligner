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

#include "fixture_base.hpp"

// ----------------------------------------------------------------------------
// Helper macro to define test values
// ----------------------------------------------------------------------------

#define DEFINE_TEST_VALUES(name, score_t, ...)  \
static auto name = [](){ return alignment::test::simd::values{.score_v = score_t{}__VA_OPT__(,) __VA_ARGS__}; }();

namespace alignment::test::simd {

// ----------------------------------------------------------------------------
// Values for test fixture
// ----------------------------------------------------------------------------

template <typename score_t, typename base_configurator_t, typename score_configurator_t>
struct values
{
    using score_type =  score_t;

    score_t score_v;
    base_configurator_t base_configurator;
    score_configurator_t score_configurator;
    std::pair<int32_t, int32_t> substitution_scores;
    std::tuple<std::size_t, std::size_t, std::size_t> sequence_generation_param;
    unsigned seed{42};
};

// ----------------------------------------------------------------------------
// Define test fixture
// ----------------------------------------------------------------------------

template <typename fixture_t>
struct test : public fixture_t
{
    using test_values_type = typename fixture_t::test_values_type;
    using scalar_score_type = typename test_values_type::score_type;

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
        constexpr std::array dna_symbol_table{'A', 'C', 'G', 'T'};

        auto [count, min_size, max_size] = this->GetParam().sequence_generation_param;

        sequence_collection1.resize(count);
        sequence_collection2.resize(count);

        std::uniform_int_distribution<size_t> sequence_size_distribution{min_size, max_size};
        std::uniform_int_distribution<size_t> symbol_distribution{0, dna_symbol_table.size() - 1};

        auto generate_sequence = [&] (auto & sequence) {
            sequence.resize(sequence_size_distribution(random_engine));
            std::ranges::generate(sequence, [&] () { return dna_symbol_table[symbol_distribution(random_engine)]; });
        };

        for (size_t i = 0; i < sequence_collection1.size(); ++i) {
            generate_sequence(sequence_collection1[i]);
            generate_sequence(sequence_collection2[i]);
        }
    }
};

} // namespace alignment::test::simd

template <typename fixture_t>
using test_suite = alignment::test::simd::test<fixture_t>;

TYPED_TEST_SUITE_P(test_suite);

// ----------------------------------------------------------------------------
// Define test cases
// ----------------------------------------------------------------------------

TYPED_TEST_P(test_suite, score)
{
    auto [match_score, mismatch_score] = this->GetParam().substitution_scores;

    auto scalar_aligner = seqan::pairwise_aligner::cfg::configure_aligner(
        seqan::pairwise_aligner::cfg::score_model_unitary(this->GetParam().base_configurator,
                                                          match_score,
                                                          mismatch_score)
    );

    auto simd_aligner = seqan::pairwise_aligner::cfg::configure_aligner(
        this->GetParam().score_configurator(this->GetParam().base_configurator,
                                            static_cast<TestFixture::scalar_score_type>(match_score),
                                            static_cast<TestFixture::scalar_score_type>(mismatch_score))
    );

    auto simd_results = simd_aligner.compute(this->sequence_collection1, this->sequence_collection2);
    size_t index = 0;
    for (auto result : simd_results) {
        EXPECT_EQ(result.score(), (scalar_aligner.compute(result.sequence1(), result.sequence2()).score()))
            << "index: " << index << "\n"
            << "s1: (" << result.sequence1().size() << ") " << result.sequence1() << "\n"
            << "s2: (" << result.sequence2().size() << ") " << result.sequence2() << "\n";
        ++index;
    }
}

// ----------------------------------------------------------------------------
// Register test cases
// ----------------------------------------------------------------------------

REGISTER_TYPED_TEST_SUITE_P(test_suite, score);
