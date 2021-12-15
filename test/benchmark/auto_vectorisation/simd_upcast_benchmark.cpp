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

template <typename lo_score_t, typename hi_score_t>
void simd_upcast_auto(benchmark::State& state) {
    namespace pa = seqan::pairwise_aligner;

    using lo_simd_t = pa::simd_score<lo_score_t>;
    constexpr size_t size = lo_simd_t::size;
    using hi_simd_t = pa::simd_score<hi_score_t, size>;

    lo_simd_t a{};
    for (size_t i = 0; i < size; ++i)
    {
        a[i] = (std::rand() % (sizeof(lo_score_t) << 3));
    }

    hi_simd_t c{};
    for (auto _ : state) {
        c = hi_simd_t{a};
    }

    int32_t score{};
    for (size_t i = 0; i < size; ++i)
        score += c[i];

    state.counters["score"] = score;
}

// template <typename score_t>
// void simd_upcast_intrinsics(benchmark::State& state) {
//     namespace pa = seqan::pairwise_aligner;

//     constexpr size_t size = 32;
//     using simd_type = seqan3::simd_type_t<score_t, size>;
//     simd_type a{};
//     simd_type b{};

//     for (size_t i = 0; i < size; ++i)
//     {
//         a[i] = (std::rand() % (sizeof(score_t) << 3));
//         b[i] = (std::rand() % (sizeof(score_t) << 3));
//     }

//     simd_type c{};
//     for (auto _ : state) {
//         c = a + b;
//     }

//     int32_t score{};
//     for (size_t i = 0; i < size; ++i)
//         score += c[i];

//     state.counters["score"] = score;
// }

// C++11 or newer, you can use the BENCHMARK macro with template parameters:
BENCHMARK_TEMPLATE(simd_upcast_auto, int8_t,  int16_t);
BENCHMARK_TEMPLATE(simd_upcast_auto, int8_t,  int32_t);
BENCHMARK_TEMPLATE(simd_upcast_auto, int16_t, int32_t);
BENCHMARK_TEMPLATE(simd_upcast_auto, uint8_t,  uint16_t);
BENCHMARK_TEMPLATE(simd_upcast_auto, uint8_t,  uint32_t);
BENCHMARK_TEMPLATE(simd_upcast_auto, uint16_t, uint32_t);

// BENCHMARK_TEMPLATE(simd_add_intrinsics, int32_t);
// BENCHMARK_TEMPLATE(simd_add_intrinsics, int16_t);
// BENCHMARK_TEMPLATE(simd_add_intrinsics, int8_t);
