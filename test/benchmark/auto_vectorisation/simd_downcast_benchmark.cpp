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

// inline constexpr size_t size = 32;

template <typename hi_score_t, typename lo_score_t>
void simd_downcast_auto(benchmark::State& state) {
    namespace pa = seqan::pairwise_aligner;

    using lo_simd_t = pa::simd_score<lo_score_t>;
    constexpr size_t size = lo_simd_t::size;
    using hi_simd_t = pa::simd_score<hi_score_t, size>;

    hi_simd_t a{};
    for (size_t i = 0; i < size; ++i)
    {
        a[i] = (std::rand() % (sizeof(lo_score_t) << 3));
    }

    lo_simd_t c{};
    for (auto _ : state) {
        c = lo_simd_t{a};
    }

    int32_t score{};
    for (size_t i = 0; i < size; ++i)
        score += c[i];

    state.counters["score"] = score;
}

BENCHMARK_TEMPLATE(simd_downcast_auto, int16_t, int8_t);
BENCHMARK_TEMPLATE(simd_downcast_auto, int32_t, int8_t);
BENCHMARK_TEMPLATE(simd_downcast_auto, int32_t, int16_t);
BENCHMARK_TEMPLATE(simd_downcast_auto, uint16_t, uint8_t);
BENCHMARK_TEMPLATE(simd_downcast_auto, uint32_t, uint8_t);
BENCHMARK_TEMPLATE(simd_downcast_auto, uint32_t, uint16_t);
