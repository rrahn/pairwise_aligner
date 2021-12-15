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

#include <pairwise_aligner/configuration/configure_aligner_saturated.hpp>
#include <pairwise_aligner/configuration/score_model_unitary_simd.hpp>
#include <pairwise_aligner/configuration/gap_model_affine.hpp>

void alignment_global_affine_bulk_simd_saturated(benchmark::State & state)
{
    size_t sequence_length = 150;
    auto seq_collection_tmp = seqan3::test::generate_sequence_pairs<seqan3::dna4>(sequence_length, 32, 0);

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

    auto aligner = pa::cfg::configure_aligner_saturated(
        pa::cfg::gap_model_affine(
            pa::cfg::score_model_unitary_simd(static_cast<int16_t>(4), static_cast<int16_t>(-5)),
            -10, -1
        )
    );

    int32_t score{};

    for (auto _ : state)
        for (auto const & res : aligner.compute(seq1_collection, seq2_collection))
            score += res.score();

    state.counters["score"] = score;
    state.counters["cells"] = seqan3::test::pairwise_cell_updates(seq_collection_tmp,
                                                                  seqan3::configuration{} |
                                                                  seqan3::align_cfg::method_global{});
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
}

BENCHMARK(alignment_global_affine_bulk_simd_saturated);
