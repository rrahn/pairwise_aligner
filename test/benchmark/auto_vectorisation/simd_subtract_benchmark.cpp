// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <ranges>
#include <string_view>
#include <seqan3/utility/simd/all.hpp>

#include <pairwise_aligner/simd/simd_score_type.hpp>

inline constexpr size_t size = seqan::pairwise_aligner::detail::max_simd_size;

template <typename score_t>
void simd_subtract(benchmark::State& state) {
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
        c = a - b;
    }

    int32_t score{};
    for (size_t i = 0; i < size; ++i)
        score += c[i];

    state.counters["score"] = score;
}

template <typename score_t>
void simd_subtract_self(benchmark::State& state) {
    namespace pa = seqan::pairwise_aligner;

    using simd_type = pa::simd_score<score_t, size>;
    simd_type a{};
    simd_type b{};

    for (size_t i = 0; i < size; ++i)
    {
        a[i] = (std::rand() % (sizeof(score_t) << 3));
        b[i] = (std::rand() % (sizeof(score_t) << 3));
    }

    for (auto _ : state) {
        a -= b;
    }

    int32_t score{};
    for (size_t i = 0; i < size; ++i)
        score += a[i];

    state.counters["score"] = score;
}

template <typename score_t>
void simd_subtract_constant(benchmark::State& state) {
    namespace pa = seqan::pairwise_aligner;

    using simd_type = seqan3::simd_type_t<score_t, size>;
    simd_type a{};
    score_t constant = (std::rand() % (sizeof(score_t) << 3));

    for (size_t i = 0; i < size; ++i)
        a[i] = (std::rand() % (sizeof(score_t) << 3));

    simd_type c{};
    for (auto _ : state) {
        c = a - constant;
    }

    int32_t score{};
    for (size_t i = 0; i < size; ++i)
        score += c[i];

    state.counters["score"] = score;
}

template <typename score_t>
void simd_subtract_constant_self(benchmark::State& state) {
    namespace pa = seqan::pairwise_aligner;

    using simd_type = seqan3::simd_type_t<score_t, size>;
    simd_type a{};
    score_t constant = (std::rand() % (sizeof(score_t) << 3));

    for (size_t i = 0; i < size; ++i)
        a[i] = (std::rand() % (sizeof(score_t) << 3));

    for (auto _ : state) {
        a -= constant;
    }

    int32_t score{};
    for (size_t i = 0; i < size; ++i)
        score += a[i];

    state.counters["score"] = score;
}

// C++11 or newer, you can use the BENCHMARK macro with template parameters:
BENCHMARK_TEMPLATE(simd_subtract, int32_t);
BENCHMARK_TEMPLATE(simd_subtract, int16_t);
BENCHMARK_TEMPLATE(simd_subtract, int8_t);

BENCHMARK_TEMPLATE(simd_subtract_self, int32_t);
BENCHMARK_TEMPLATE(simd_subtract_self, int16_t);
BENCHMARK_TEMPLATE(simd_subtract_self, int8_t);

BENCHMARK_TEMPLATE(simd_subtract_constant, int32_t);
BENCHMARK_TEMPLATE(simd_subtract_constant, int16_t);
BENCHMARK_TEMPLATE(simd_subtract_constant, int8_t);

BENCHMARK_TEMPLATE(simd_subtract_constant_self, int32_t);
BENCHMARK_TEMPLATE(simd_subtract_constant_self, int16_t);
BENCHMARK_TEMPLATE(simd_subtract_constant_self, int8_t);
