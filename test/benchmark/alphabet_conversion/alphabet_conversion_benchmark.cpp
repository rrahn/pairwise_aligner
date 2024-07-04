// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <algorithm>
#include <random>
#include <ranges>
#include <string_view>
#include <tuple>

#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/aminoacid/aa20.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/utility/simd/all.hpp>

#include <pairwise_aligner/simd/simd_score_type.hpp>
#include <pairwise_aligner/alphabet_conversion/alphabet_rank_map_scalar.hpp>
#include <pairwise_aligner/alphabet_conversion/alphabet_rank_map_simd.hpp>

inline constexpr size_t sequence_count = seqan::pairwise_aligner::detail::max_simd_size;
using rank_t = int8_t;

inline auto generate(size_t const sequence_size, std::string_view symbol_list)
{
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<char> symbol_distribution(0, symbol_list.size() - 1);

    std::vector<std::string> sequences;
    sequences.resize(sequence_count);

    std::ranges::for_each(sequences, [&] (std::string & sequence) {
        sequence.resize(sequence_size);
        std::ranges::generate(sequence, [&] () -> char { return symbol_distribution(gen); });
    });

    return sequences;
}

template <typename ...args_t>
void alphabet_conversion_scalar(benchmark::State& state, args_t&&... args) {
    std::string_view symbol_list = get<0>(std::make_tuple(std::move(args)...));
    size_t const sequence_size = state.range(0);

    // Generate sequences
    auto sequence_collection = generate(sequence_size, symbol_list);

    // Prepare memory for conversion.
    std::vector<std::vector<rank_t>> ranks{};
    ranks.resize(sequence_collection.size());
    std::ranges::for_each(ranks, [&] (std::vector<rank_t> & rank_vec) { rank_vec.resize(sequence_size); });

    // Initialise rank map.
    seqan::pairwise_aligner::alphabet_rank_map_scalar rank_map{symbol_list};

    for (auto _ : state) {
        for (size_t i = 0; i < sequence_count; ++i) {
            std::ranges::copy(sequence_collection[i] | std::views::transform([&] (auto const & symbol) {
                return rank_map[symbol];
            }), ranks[i].begin());
        }
    }

    // Output bytes per second.
    state.SetBytesProcessed(int64_t(state.iterations()) * int64_t(state.range(0)) * int64_t(sequence_count));
}

template <typename ...args_t>
void alphabet_conversion_simd(benchmark::State& state, args_t&&... args) {
    std::string_view symbol_list = get<0>(std::make_tuple(std::move(args)...));
    size_t const sequence_size = state.range(0);

    // Generate sequences
    auto sequence_collection = generate(sequence_size, symbol_list);

    // Prepare memory for conversion.
    using simd_t = seqan::pairwise_aligner::simd_score<int8_t>;
    using native_simd_t = typename simd_t::native_simd_type;

    auto to_simd = [&] (auto & _sequence_collection) {

        std::vector<simd_t, seqan3::aligned_allocator<simd_t, alignof(simd_t)>> sequences{};
        sequences.reserve(sequence_size);

        auto simd_view = _sequence_collection | seqan3::views::to_simd<native_simd_t>();

        for (auto && simd_vector_chunk : simd_view) {
            for (auto && simd_vector : simd_vector_chunk) {
                sequences.emplace_back(std::move(simd_vector));
            }
        }
        return sequences;
    };

    auto sequences = to_simd(sequence_collection);
    std::vector<simd_t, seqan3::aligned_allocator<simd_t, alignof(simd_t)>> ranks{};
    ranks.resize(sequence_size);

    // Initialise rank map.
    seqan::pairwise_aligner::alphabet_rank_map_simd<simd_t> rank_map{symbol_list};

    for (auto _ : state) {
        std::ranges::copy(sequences | std::views::transform([&] (auto const & symbol) {
            return rank_map[symbol];
        }), sequences.begin());
    }

    // Output bytes per second.
    state.SetBytesProcessed(int64_t(state.iterations()) * int64_t(state.range(0)) * int64_t(sequence_count));
}

template <typename ...args_t>
void alphabet_conversion_seqan3(benchmark::State& state, args_t&&... args) {
    auto tuple_args = std::make_tuple(std::move(args)...);
    std::string_view symbol_list = get<0>(tuple_args);

    using alphabet_t = std::tuple_element_t<1, decltype(tuple_args)>;

    size_t const sequence_size = state.range(0);

    // Generate sequences
    auto sequence_collection = generate(sequence_size, symbol_list);

    // Prepare memory for conversion.
    std::vector<std::vector<alphabet_t>> ranks{};
    ranks.resize(sequence_collection.size());
    std::ranges::for_each(ranks, [&] (std::vector<alphabet_t> & rank_vec) { rank_vec.resize(sequence_size); });

    // Initialise rank map.
    seqan::pairwise_aligner::alphabet_rank_map_scalar rank_map{symbol_list};

    for (auto _ : state) {
        for (size_t i = 0; i < sequence_count; ++i) {
            std::ranges::copy(sequence_collection[i] | seqan3::views::char_to<alphabet_t>, ranks[i].begin());
        }
    }

    // Output bytes per second.
    state.SetBytesProcessed(int64_t(state.iterations()) * int64_t(state.range(0)) * int64_t(sequence_count));
}


inline constexpr std::string_view dna5_symbols{"ACGTN"};
inline constexpr std::string_view aa20_symbols{"ACDEFGHIKLMNPQRSTVWY"};
inline constexpr std::string_view aa27_symbols{"ABCDEFGHIJKLMNOPQRSTUVWXYZ@"};
inline constexpr std::string_view printable_symbols{" !\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"};

BENCHMARK_CAPTURE(alphabet_conversion_scalar, dna5, dna5_symbols)->Arg(100)->Arg(250)->Arg(500)->Arg(1000);
BENCHMARK_CAPTURE(alphabet_conversion_scalar, aa20, aa20_symbols)->Arg(100)->Arg(250)->Arg(500)->Arg(1000);
BENCHMARK_CAPTURE(alphabet_conversion_scalar, aa27, aa27_symbols)->Arg(100)->Arg(250)->Arg(500)->Arg(1000);
BENCHMARK_CAPTURE(alphabet_conversion_scalar, printable_char, printable_symbols)->Arg(100)->Arg(250)->Arg(500)->Arg(1000);

BENCHMARK_CAPTURE(alphabet_conversion_simd, dna5, dna5_symbols)->Arg(100)->Arg(250)->Arg(500)->Arg(1000);
BENCHMARK_CAPTURE(alphabet_conversion_simd, aa20, aa20_symbols)->Arg(100)->Arg(250)->Arg(500)->Arg(1000);
BENCHMARK_CAPTURE(alphabet_conversion_simd, aa27, aa27_symbols)->Arg(100)->Arg(250)->Arg(500)->Arg(1000);
BENCHMARK_CAPTURE(alphabet_conversion_simd, printable_char, printable_symbols)->Arg(100)->Arg(250)->Arg(500)->Arg(1000);

BENCHMARK_CAPTURE(alphabet_conversion_seqan3, dna5, dna5_symbols, seqan3::dna5{})->Arg(100)->Arg(250)->Arg(500)->Arg(1000);
BENCHMARK_CAPTURE(alphabet_conversion_seqan3, aa20, aa20_symbols, seqan3::aa20{})->Arg(100)->Arg(250)->Arg(500)->Arg(1000);
BENCHMARK_CAPTURE(alphabet_conversion_seqan3, aa27, aa27_symbols, seqan3::aa27{})->Arg(100)->Arg(250)->Arg(500)->Arg(1000);
