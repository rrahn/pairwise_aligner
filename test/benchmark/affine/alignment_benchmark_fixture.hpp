// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <memory>
#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/performance/units.hpp>
#include <seqan3/utility/views/zip.hpp>

#include <pairwise_aligner/configuration/configure_aligner.hpp>

// ----------------------------------------------------------------------------
// Helper macro to define benchmark values
// ----------------------------------------------------------------------------

#define DEFINE_BENCHMARK_VALUES(name, ...)  \
static auto name = [](){ return aligner::benchmark::values{__VA_ARGS__}; }();


#define ALIGNER_BENCHMARK(name, type)  \
BENCHMARK_TEMPLATE_F(test, name##_##type, aligner::benchmark::fixture<&type>)(::benchmark::State & state) \
{ this->run(state); }

namespace aligner::benchmark
{

// ----------------------------------------------------------------------------
// Templatized fixture values.
// ----------------------------------------------------------------------------

template <typename configurator_t, typename seqan_configurator_t>
struct values
{
    configurator_t configurator;
    seqan_configurator_t seqan_configurator;

    size_t sequence_size_mean;
    size_t sequence_size_variance;
    size_t sequence_count;
};

// ----------------------------------------------------------------------------
// Helper struct to capture fixture values as template argument
// ----------------------------------------------------------------------------

template <auto _fixture>
struct fixture : public ::benchmark::Fixture
{
    auto GetParam() const noexcept
        -> decltype(*_fixture) const &
    {
        return *_fixture;
    }
};

// ----------------------------------------------------------------------------
// Benchmark fixture to generate values for each benchmark.
// ----------------------------------------------------------------------------

template <typename fixture_t>
class test : public fixture_t
{
public:
    using collection_t = std::vector<std::string>;

    std::unique_ptr<collection_t> _seq1_collection{};
    std::unique_ptr<collection_t> _seq2_collection{};

    using fixture_t::GetParam;

    void SetUp([[maybe_unused]] ::benchmark::State & state) override {

        if (!_seq1_collection) {

            auto seq_collection_tmp =
                seqan3::test::generate_sequence_pairs<seqan3::dna4>(GetParam().sequence_size_mean,
                                                                    GetParam().sequence_count,
                                                                    GetParam().sequence_size_variance);

            collection_t seq1_collection{};
            collection_t seq2_collection{};
            for (auto const & seq_tmp : seq_collection_tmp)
            {
                auto char_seq1 = seq_tmp.first | seqan3::views::to_char;
                auto char_seq2 = seq_tmp.second | seqan3::views::to_char;

                seq1_collection.emplace_back(std::ranges::begin(char_seq1), std::ranges::end(char_seq1));
                seq2_collection.emplace_back(std::ranges::begin(char_seq2), std::ranges::end(char_seq2));
            }

            _seq1_collection = std::make_unique<collection_t>(std::move(seq1_collection));
            _seq2_collection = std::make_unique<collection_t>(std::move(seq2_collection));
        }

        if (!_seq1_collection || !_seq2_collection)
            throw std::runtime_error{"Setting up the benchmark failed!"};
    }

    auto & sequence1() noexcept { return *_seq1_collection; }

    auto & sequence2() noexcept { return *_seq2_collection; }

public:
    void run(::benchmark::State & state) noexcept {

        auto aligner = seqan::pairwise_aligner::cfg::configure_aligner(GetParam().configurator);

        int32_t score{};

        for (auto _ : state)
            for (auto const & res : aligner.compute(sequence1(), sequence2()))
                score += res.score();

        state.counters["score"] = score;
        state.counters["cells"] = seqan3::test::pairwise_cell_updates(seqan3::views::zip(sequence1(),
                                                                                         sequence2()),
                                                                    GetParam().seqan_configurator);
        state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
    }
};

} // namespace aligner::benchmark

template <typename values_t>
using test = aligner::benchmark::test<values_t>;
