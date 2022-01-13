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
#include <string_view>
#include <string>
#include <variant>
#include <vector>

#include <pairwise_aligner/configuration/configure_aligner.hpp>
#include <pairwise_aligner/configuration/gap_model_affine.hpp>
#include <pairwise_aligner/configuration/method_global.hpp>
#include <pairwise_aligner/configuration/score_model_unitary_simd.hpp>
#include <pairwise_aligner/configuration/score_model_unitary.hpp>

namespace aligner = seqan::pairwise_aligner;

template <typename score_t>
struct alignment_test_values
{
    std::tuple<std::size_t, std::size_t, std::size_t> sequence_generation_param{};
    std::pair<int32_t, int32_t> substitution_scores{};
    std::pair<int32_t, int32_t> gap_scores{};
    aligner::initialisation_rule leading_gaps{};
    aligner::trailing_gap_setting trailing_gaps{};
    unsigned seed{42};

};

struct alignment_test_value_wrapper
{
    using test_value_types = std::variant<alignment_test_values<int8_t>,
                                          alignment_test_values<int16_t>,
                                          alignment_test_values<int32_t>,
                                          alignment_test_values<int64_t>>;

    test_value_types values{};

    auto sequence_generation_param() const noexcept
    {
        return std::visit([] (auto value) { return value.sequence_generation_param; }, values);
    }

    auto substitution_scores() const noexcept
    {
        return std::visit([] (auto value) { return value.substitution_scores; }, values);
    }

    auto gap_scores() const noexcept
    {
        return std::visit([] (auto value) { return value.gap_scores; }, values);
    }

    auto leading_gaps() const noexcept
    {
        return std::visit([] (auto value) { return value.leading_gaps; }, values);
    }

    auto trailing_gaps() const noexcept
    {
        return std::visit([] (auto value) { return value.trailing_gaps; }, values);
    }

    unsigned seed() const noexcept
    {
        return std::visit([] (auto value) { return value.seed; }, values);
    }
};

struct simd_alignment_test : public ::testing::TestWithParam<alignment_test_value_wrapper>
{
    std::vector<std::string> sequence_collection1{};
    std::vector<std::string> sequence_collection2{};

    void SetUp() override
    {
        std::mt19937 random_engine{GetParam().seed()};
        generate_sequences(random_engine);
    }

    template <typename fn_t>
    void run_test(fn_t && fn) const
    {
        std::visit([&] (auto const & value) { std::invoke(fn, value); } , GetParam().values);
    }

private:

    template <typename random_engine_t>
    void generate_sequences(random_engine_t & random_engine)
    {
        auto [count, min_size, max_size] = GetParam().sequence_generation_param();
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

TEST_P(simd_alignment_test, compute_score)
{
    auto [match_score, mismatch_score] = GetParam().substitution_scores();
    auto [gap_open_score, gap_extension_score] = GetParam().gap_scores();

    auto aligner_base = aligner::cfg::method_global(
        aligner::cfg::gap_model_affine(gap_open_score, gap_extension_score),
        GetParam().leading_gaps(), GetParam().trailing_gaps()
    );

    auto single_aligner = aligner::cfg::configure_aligner(
        aligner::cfg::score_model_unitary(aligner_base, match_score, mismatch_score)
    );

    run_test([&] <typename score_t> (alignment_test_values<score_t> const &) {

        auto simd_aligner = aligner::cfg::configure_aligner(
            aligner::cfg::score_model_unitary_simd(aligner_base,
                                                   static_cast<score_t>(match_score),
                                                   static_cast<score_t>(mismatch_score))
        );

        auto simd_results = simd_aligner.compute(sequence_collection1, sequence_collection2);
        for (auto result : simd_results) {
            EXPECT_EQ(result.score(), (single_aligner.compute(result.sequence1(), result.sequence2()).score()));
        }
    });
}

// ----------------------------------------------------------------------------
// Test cases
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Regular alignment

INSTANTIATE_TEST_SUITE_P(regular_equal_size_64, simd_alignment_test, testing::Values(alignment_test_values<int64_t>
{
    .sequence_generation_param{aligner::simd_score<int64_t>::size, 93, 93},
    .substitution_scores{4, -5},
    .gap_scores{-10, -1},
    .leading_gaps{},
    .trailing_gaps{}
}));

INSTANTIATE_TEST_SUITE_P(regular_equal_size_32, simd_alignment_test, testing::Values(alignment_test_values<int32_t>
{
    .sequence_generation_param{aligner::simd_score<int32_t>::size, 210, 210},
    .substitution_scores{4, -5},
    .gap_scores{-10, -1},
    .leading_gaps{},
    .trailing_gaps{}
}));

INSTANTIATE_TEST_SUITE_P(regular_equal_size_16, simd_alignment_test, testing::Values(alignment_test_values<int16_t>
{
    .sequence_generation_param{aligner::simd_score<int16_t>::size, 150, 150},
    .substitution_scores{4, -5},
    .gap_scores{-10, -1},
    .leading_gaps{},
    .trailing_gaps{}
}));

INSTANTIATE_TEST_SUITE_P(regular_equal_size_8, simd_alignment_test, testing::Values(alignment_test_values<int8_t>
{
    .sequence_generation_param{aligner::simd_score<int8_t>::size, 25, 25},
    .substitution_scores{4, -5},
    .gap_scores{-10, -1},
    .leading_gaps{},
    .trailing_gaps{}
}));

// ----------------------------------------------------------------------------
// Semi-global with free end-gaps in sequence 1

INSTANTIATE_TEST_SUITE_P(semi_seq1_equal_size_32, simd_alignment_test, testing::Values(alignment_test_values<int32_t>
{
    .sequence_generation_param{aligner::simd_score<int32_t>::size, 333, 333},
    .substitution_scores{4, -5},
    .gap_scores{-10, -1},
    .leading_gaps{.column_initialisation = aligner::dp_initialisation_rule::zero},
    .trailing_gaps{.column = aligner::dp_trailing_gaps::free}
}));

// ----------------------------------------------------------------------------
// Semi-global with free end-gaps in sequence 2

INSTANTIATE_TEST_SUITE_P(semi_seq2_equal_size_16, simd_alignment_test, testing::Values(alignment_test_values<int16_t>
{
    .sequence_generation_param{aligner::simd_score<int16_t>::size, 70, 70},
    .substitution_scores{4, -5},
    .gap_scores{-10, -1},
    .leading_gaps{.row_initialisation = aligner::dp_initialisation_rule::zero},
    .trailing_gaps{.row = aligner::dp_trailing_gaps::free}
}));

// ----------------------------------------------------------------------------
// Overlap

INSTANTIATE_TEST_SUITE_P(overlap_equal_size_16, simd_alignment_test, testing::Values(alignment_test_values<int16_t>
{
    .sequence_generation_param{aligner::simd_score<int16_t>::size, 12, 12},
    .substitution_scores{4, -5},
    .gap_scores{-10, -1},
    .leading_gaps{.column_initialisation = aligner::dp_initialisation_rule::zero,
                  .row_initialisation = aligner::dp_initialisation_rule::zero},
    .trailing_gaps{.column = aligner::dp_trailing_gaps::free, .row = aligner::dp_trailing_gaps::free}
}));
