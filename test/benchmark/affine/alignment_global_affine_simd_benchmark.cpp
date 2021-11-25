// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <seqan3/std/ranges>
#include <string_view>

#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/performance/units.hpp>
#include <seqan3/core/configuration/configuration.hpp>

#include <pairwise_aligner/affine/affine_dp_algorithm.hpp>
#include <pairwise_aligner/affine/gap_model_affine.hpp>
#include <pairwise_aligner/affine/initialisation_strategy_affine.hpp>
#include <pairwise_aligner/dp_initialisation_rule.hpp>

void alignment_global_affine_bulk_scalar(benchmark::State & state)
{
    size_t sequence_length = 150;
    auto seq_collection_tmp = seqan3::test::generate_sequence_pairs<seqan3::dna4>(sequence_length, 16, 0);

    std::vector<std::string> seq1_collection{};

    for (auto const & seq_tmp : seq_collection_tmp)
    {
        auto char_seq = seq_tmp.first | seqan3::views::to_char;
        seq1_collection.emplace_back(std::ranges::begin(char_seq), std::ranges::end(char_seq));
    }

    std::vector<std::string> seq2_collection{};
    for (auto const &  seq_tmp : seq_collection_tmp)
    {
        auto char_seq = seq_tmp.second | seqan3::views::to_char;
        seq2_collection.emplace_back(std::ranges::begin(char_seq), std::ranges::end(char_seq));
    }

    namespace pa = seqan::pairwise_aligner;
    pa::gap_model_affine<int32_t> gap_model{-10, -1};
    pa::initialisation_strategy_affine init{gap_model,
                                            pa::dp_initialisation_rule::regular,
                                            pa::dp_initialisation_rule::regular};

    pa::pairwise_aligner_affine aligner{gap_model, init};
    int32_t score{};

    for (auto _ : state)
        for (size_t i = 0; i < seq_collection_tmp.size(); ++i)
            score += aligner.compute(seq1_collection[i], seq2_collection[i]);

    state.counters["score"] = score;
    state.counters["cells"] = seqan3::test::pairwise_cell_updates(seq_collection_tmp,
                                                                  seqan3::configuration{} |
                                                                  seqan3::align_cfg::method_global{});
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
}

void alignment_global_affine_bulk_simd(benchmark::State & state)
{
    size_t sequence_length = 150;
    auto seq_collection_tmp = seqan3::test::generate_sequence_pairs<seqan3::dna4>(sequence_length, 16, 0);

    std::vector<std::string> seq1_collection{};

    for (auto const & seq_tmp : seq_collection_tmp)
    {
        auto char_seq = seq_tmp.first | seqan3::views::to_char;
        seq1_collection.emplace_back(std::ranges::begin(char_seq), std::ranges::end(char_seq));
    }

    std::vector<std::string> seq2_collection{};
    for (auto const &  seq_tmp : seq_collection_tmp)
    {
        auto char_seq = seq_tmp.second | seqan3::views::to_char;
        seq2_collection.emplace_back(std::ranges::begin(char_seq), std::ranges::end(char_seq));
    }

    namespace pa = seqan::pairwise_aligner;
    pa::gap_model_affine<int32_t> gap_model{-10, -1};
    pa::initialisation_strategy_affine init{gap_model,
                                            pa::dp_initialisation_rule::regular,
                                            pa::dp_initialisation_rule::regular};

    pa::pairwise_aligner_affine aligner{gap_model, init};
    int32_t score{};

    for (auto _ : state)
        score += aligner.compute_simd(seq1_collection, seq2_collection)[0];

    state.counters["score"] = score;
    state.counters["cells"] = seqan3::test::pairwise_cell_updates(seq_collection_tmp,
                                                                  seqan3::configuration{} |
                                                                  seqan3::align_cfg::method_global{});
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
}

BENCHMARK(alignment_global_affine_bulk_scalar);
BENCHMARK(alignment_global_affine_bulk_simd);
