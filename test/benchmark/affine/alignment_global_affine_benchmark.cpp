// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <ranges>
#include <string_view>

#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/performance/units.hpp>
#include <seqan3/core/configuration/configuration.hpp>

#include <pairwise_aligner/configuration/configure_aligner.hpp>
#include <pairwise_aligner/configuration/score_model_unitary.hpp>
#include <pairwise_aligner/configuration/gap_model_affine.hpp>

#include "alignment_benchmark_fixture.hpp"

void alignment_global_affine(benchmark::State & state)
{
    size_t sequence_length = aligner::benchmark::sequence_size;

    auto generated_seq1 = seqan3::test::generate_sequence<seqan3::dna4>(sequence_length, 0, 0);
    auto generated_seq2 = seqan3::test::generate_sequence<seqan3::dna4>(sequence_length, 0, 1);

    auto char_seq1_view = generated_seq1 | seqan3::views::to_char;
    auto char_seq2_view = generated_seq2 | seqan3::views::to_char;

    std::string seq1{std::ranges::begin(char_seq1_view), std::ranges::end(char_seq1_view)};
    std::string seq2{std::ranges::begin(char_seq2_view), std::ranges::end(char_seq2_view)};

    namespace pa = seqan::pairwise_aligner;

    auto aligner = pa::cfg::configure_aligner(
        pa::cfg::gap_model_affine(
            pa::cfg::score_model_unitary(4, -5),
            -10, -1
        )
    );
    int32_t score{};

    for (auto _ : state)
        score += aligner.compute(seq1, seq2).score();

    state.counters["score"] = score;
    state.counters["cells"] = seqan3::test::pairwise_cell_updates(std::views::single(std::tie(seq1, seq2)),
                                                                  seqan3::configuration{} |
                                                                  seqan3::align_cfg::method_global{});
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
}

BENCHMARK(alignment_global_affine);

BENCHMARK_MAIN();
