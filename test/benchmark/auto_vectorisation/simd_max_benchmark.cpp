// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <seqan3/std/ranges>
#include <string_view>
#include <seqan3/utility/simd/all.hpp>

#include <pairwise_aligner/simd/simd_score_type.hpp>

inline constexpr size_t size = seqan::pairwise_aligner::detail::max_simd_size;

template <typename score_t>
void simd_max(benchmark::State& state) {
    namespace pa = seqan::pairwise_aligner;

    using simd_type = pa::simd_score<score_t, size>;
    simd_type a{};
    simd_type b{};

    for (size_t i = 0; i < size; ++i)
    {
        a[i] = (std::rand() % (sizeof(score_t) << 3));
        b[i] = (std::rand() % (sizeof(score_t) << 3));
    }

    simd_type c{};
    for (auto _ : state) {
        c = max(a, b);
    }

    int32_t score{};
    for (size_t i = 0; i < size; ++i)
        score += c[i];

    state.counters["score"] = score;
}

// C++11 or newer, you can use the BENCHMARK macro with template parameters:
BENCHMARK_TEMPLATE(simd_max, int32_t);
BENCHMARK_TEMPLATE(simd_max, int16_t);
BENCHMARK_TEMPLATE(simd_max, int8_t);
