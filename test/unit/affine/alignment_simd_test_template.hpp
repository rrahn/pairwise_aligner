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
#include <seqan3/std/ranges>
#include <string_view>
#include <string>
#include <vector>

#include <seqan3/core/debug_stream.hpp>

#include <pairwise_aligner/configuration/configure_aligner.hpp>
#include <pairwise_aligner/configuration/score_model_matrix.hpp>
#include <pairwise_aligner/configuration/score_model_unitary.hpp>

#include "../fixture/fixture_base.hpp"

// ----------------------------------------------------------------------------
// Helper macro to define test values
// ----------------------------------------------------------------------------

#define DEFINE_TEST_VALUES(name, ...)  \
static auto name = [](){ return alignment::test::simd::values{__VA_ARGS__}; }();

namespace alignment::test::simd {

using namespace std::literals;
// ----------------------------------------------------------------------------
// Values for test fixture
// ----------------------------------------------------------------------------

template <typename score_t>
struct unitary_model{
    using score_type = score_t;

    score_t match_score{};
    score_t mismatch_score{};

    constexpr auto symbol_list() const noexcept
    {
        return "ACGT"sv;
    }

    constexpr auto compare_configurator() const noexcept
    {
        return seqan::pairwise_aligner::cfg::score_model_unitary;
    }
};

template <typename matrix_t>
struct matrix_model{
    using score_type = typename std::tuple_element_t<1, typename matrix_t::value_type>::value_type;
    matrix_t matrix{};

    constexpr auto symbol_list() const noexcept
    {
        return matrix | std::views::elements<0>;
    }

    constexpr auto compare_configurator() const noexcept
    {
        return seqan::pairwise_aligner::cfg::score_model_matrix;
    }
};

template <typename base_configurator_t, typename score_configurator_t, typename score_model_t, typename bool_t = std::false_type>
struct values
{
    using score_type = typename score_model_t::score_type;

    base_configurator_t base_configurator;
    score_configurator_t score_configurator;
    score_model_t substitution_scores;
    std::tuple<std::size_t, std::size_t, std::size_t> sequence_generation_param;
    bool_t one_vs_many{};
    unsigned seed{42};

    constexpr auto symbol_list() const noexcept
    {
        return substitution_scores.symbol_list();
    }

    constexpr auto compare_configurator() const noexcept
    {
        return substitution_scores.compare_configurator();
    }
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

    auto const & sequence1() const noexcept
    {
        if constexpr (decltype(std::declval<fixture_t>().GetParam().one_vs_many)::value == true) {
            return sequence_collection1.front();
        } else {
            return sequence_collection1;
        }
    }

    auto const & sequence2() const noexcept
    {
        return sequence_collection2;
    }

    template <typename score_config_t, typename score_t>
    auto configure_aligner_for(score_config_t const & config, unitary_model<score_t> const & score_model) const noexcept
    {
        auto [match_score, mismatch_score] = score_model;

        return seqan::pairwise_aligner::cfg::configure_aligner(
                    config(this->GetParam().base_configurator,
                           static_cast<score_t>(match_score),
                           static_cast<score_t>(mismatch_score)));
    }

    template <typename score_config_t, typename matrix_t>
    auto configure_aligner_for(score_config_t const & config, matrix_model<matrix_t> const & score_model) const noexcept
    {
        return seqan::pairwise_aligner::cfg::configure_aligner(config(this->GetParam().base_configurator,
                                                                      score_model.matrix));
    }

private:

    template <typename random_engine_t>
    void generate_sequences(random_engine_t & random_engine)
    {
        auto symbol_table = this->GetParam().symbol_list();

        auto [count, min_size, max_size] = this->GetParam().sequence_generation_param;

        sequence_collection1.resize(count);
        sequence_collection2.resize(count);

        std::uniform_int_distribution<size_t> sequence_size_distribution{min_size, max_size};
        std::uniform_int_distribution<size_t> symbol_distribution{0, symbol_table.size() - 1};

        auto generate_sequence = [&] (auto & sequence) {
            sequence.resize(sequence_size_distribution(random_engine));
            std::ranges::generate(sequence, [&] () { return symbol_table[symbol_distribution(random_engine)]; });
        };

        for (size_t i = 0; i < sequence_collection1.size(); ++i) {
            generate_sequence(sequence_collection1[i]);
            seqan3::debug_stream << "seq1 = " << sequence_collection1[i] << "\n";
            generate_sequence(sequence_collection2[i]);
            seqan3::debug_stream << "seq2 = " << sequence_collection2[i] << "\n";
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
    // we can only configure inside of the name
    // auto [match_score, mismatch_score] = this->GetParam().substitution_scores;

    auto scalar_aligner = this->configure_aligner_for(this->GetParam().compare_configurator(),
                                                      this->GetParam().substitution_scores);
    // seqan::pairwise_aligner::cfg::configure_aligner(
    //     seqan::pairwise_aligner::cfg::score_model_unitary(this->GetParam().base_configurator,
    //                                                       static_cast<int32_t>(match_score),
    //                                                       static_cast<int32_t>(mismatch_score))
    // );

    auto simd_aligner = this->configure_aligner_for(this->GetParam().score_configurator,
                                                    this->GetParam().substitution_scores);
    // seqan::pairwise_aligner::cfg::configure_aligner(
    //     this->GetParam().score_configurator(this->GetParam().base_configurator,
    //                                         static_cast<TestFixture::scalar_score_type>(match_score),
    //                                         static_cast<TestFixture::scalar_score_type>(mismatch_score))
    // );

    auto simd_results = simd_aligner.compute(this->sequence1(), this->sequence2());
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
