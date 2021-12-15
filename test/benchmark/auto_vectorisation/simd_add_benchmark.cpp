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

#include <pairwise_aligner/simd_score_type.hpp>

// inline constexpr size_t size = 32;

template <typename score_t>
void simd_add_auto(benchmark::State& state) {
    namespace pa = seqan::pairwise_aligner;

    using simd_type = pa::simd_score<score_t>;
    constexpr size_t size = simd_type::size;
    simd_type a{};
    simd_type b{};

    for (size_t i = 0; i < size; ++i)
    {
        a[i] = (std::rand() % (sizeof(score_t) << 3));
        b[i] = (std::rand() % (sizeof(score_t) << 3));
    }

    simd_type c{};
    for (auto _ : state) {
        c = a + b;
    }

    int32_t score{};
    for (size_t i = 0; i < size; ++i)
        score += c[i];

    state.counters["score"] = score;
}

template <typename score_t>
void simd_add_intrinsics(benchmark::State& state) {
    namespace pa = seqan::pairwise_aligner;

    constexpr size_t size = 32;
    using simd_type = seqan3::simd_type_t<score_t, size>;
    simd_type a{};
    simd_type b{};

    for (size_t i = 0; i < size; ++i)
    {
        a[i] = (std::rand() % (sizeof(score_t) << 3));
        b[i] = (std::rand() % (sizeof(score_t) << 3));
    }

    simd_type c{};
    for (auto _ : state) {
        c = a + b;
    }

    int32_t score{};
    for (size_t i = 0; i < size; ++i)
        score += c[i];

    state.counters["score"] = score;
}

// C++11 or newer, you can use the BENCHMARK macro with template parameters:
BENCHMARK_TEMPLATE(simd_add_auto, int32_t);
BENCHMARK_TEMPLATE(simd_add_auto, int16_t);
BENCHMARK_TEMPLATE(simd_add_auto, int8_t);

BENCHMARK_TEMPLATE(simd_add_intrinsics, int32_t);
BENCHMARK_TEMPLATE(simd_add_intrinsics, int16_t);
BENCHMARK_TEMPLATE(simd_add_intrinsics, int8_t);
