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

void alignment_global_affine(benchmark::State & state)
{
    size_t sequence_length = 500;
    auto seq1_tmp = seqan3::test::generate_sequence<seqan3::dna4>(sequence_length, 0, 0) | seqan3::views::to_char;
    auto seq2_tmp = seqan3::test::generate_sequence<seqan3::dna4>(sequence_length, 0, 1) | seqan3::views::to_char;

    std::string seq1{std::ranges::begin(seq1_tmp), std::ranges::end(seq1_tmp)};
    std::string seq2{std::ranges::begin(seq2_tmp), std::ranges::end(seq2_tmp)};

    namespace pa = seqan::pairwise_aligner;
    pa::gap_model_affine<int32_t> gap_model{-10, -1};
    pa::initialisation_strategy_affine init{gap_model,
                                            pa::dp_initialisation_rule::regular,
                                            pa::dp_initialisation_rule::regular};

    pa::pairwise_aligner_affine aligner{gap_model, init};
    int32_t score{};

    for (auto _ : state)
        score += aligner.compute(seq1, seq2);

    state.counters["score"] = score;
    state.counters["cells"] = seqan3::test::pairwise_cell_updates(std::views::single(std::tie(seq1, seq2)),
                                                                  seqan3::configuration{} |
                                                                  seqan3::align_cfg::method_global{});
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
}

BENCHMARK(alignment_global_affine);
